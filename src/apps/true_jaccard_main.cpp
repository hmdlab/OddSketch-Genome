#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

struct PairInfo {
  int pair_id{};
  string file1;
  string file2;
  long long mutation_count{};
  long long genome_length{};
};

static inline string read_fasta_concat(const string &path) {
  ifstream ifs(path);
  if (!ifs) throw runtime_error("Cannot open FASTA: " + path);
  string line, seq;
  seq.reserve(1<<20);
  while (getline(ifs, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') continue;
    if (!line.empty() && line.back()=='\r') line.pop_back();
    seq += line;
  }
  return seq;
}

static inline size_t kmer_count_possible(size_t L, size_t k) {
  if (L < k) return 0; return L - k + 1;
}

static inline char comp_base(char c) {
  switch (c) {
    case 'A': case 'a': return 'T';
    case 'C': case 'c': return 'G';
    case 'G': case 'g': return 'C';
    case 'T': case 't': return 'A';
    default: return 'N';
  }
}

static inline string revcomp(const string &s) {
  string rc; rc.resize(s.size());
  for (size_t i = 0, n = s.size(); i < n; ++i) {
    rc[n - 1 - i] = comp_base(s[i]);
  }
  return rc;
}

static inline string canonical_kmer(const string &s) {
  string rc = revcomp(s);
  return (rc < s) ? rc : s;
}

struct JaccardOut {
  int pair_id{};
  long long mutation_count{};
  long long genome_length{};
  double mutation_rate{};
  double jaccard_true{};
  size_t kmers1_count{};
  size_t kmers2_count{};
  size_t inter{};
  size_t uni{};
};

// very lightweight config parser: find integer after key occurrences
static optional<int> find_kmer_from_config(const string &json, const string &ns_prefix) {
  // If ns_prefix provided (e.g., "\"true_jaccard\""), search within it first.
  size_t pos = 0;
  if (!ns_prefix.empty()) {
    pos = json.find(ns_prefix);
    if (pos == string::npos) pos = 0;
  }
  const string key = "\"kmerlen\"";
  size_t kpos = json.find(key, pos);
  if (kpos == string::npos) return nullopt;
  size_t colon = json.find(':', kpos);
  if (colon == string::npos) return nullopt;
  size_t beg = json.find_first_of("0123456789", colon);
  if (beg == string::npos) return nullopt;
  size_t end = json.find_first_not_of("0123456789", beg);
  try {
    int v = stoi(json.substr(beg, end - beg));
    return v;
  } catch (...) { return nullopt; }
}

static int parse_int_arg(const vector<string> &args, const string &name, int defval) {
  string key = "--" + name + "=";
  for (auto &a : args) if (a.rfind(key, 0) == 0) {
    try { return stoi(a.substr(key.size())); } catch (...) { return defval; }
  }
  return defval;
}

static string parse_str_arg(const vector<string> &args, const string &name, const string &defval) {
  string key = "--" + name + "=";
  for (auto &a : args) if (a.rfind(key, 0) == 0) return a.substr(key.size());
  return defval;
}

int main(int argc, char** argv) {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  // Defaults
  int kmer = 64;
  // sequential execution only
  string pair_info = "data/test_genomes/pair_info.txt"; // relative to working dir
  string out_path = "data/test_genomes/jaccard_true_results.txt";
  string cfg_path = ""; // optional

  vector<string> args; args.reserve(argc);
  for (int i=1;i<argc;++i) args.emplace_back(argv[i]);
  // first read CLI to know config path
  cfg_path = parse_str_arg(args, "config", cfg_path);

  // read config if available to set defaults (will be overridden by explicit CLI)
  auto apply_config = [&](const string &path){
    ifstream ifs(path);
    if (!ifs) return;
    string json((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
    if (auto v = find_kmer_from_config(json, "\"true_jaccard\"")) kmer = *v;
    else if (auto v2 = find_kmer_from_config(json, "")) kmer = *v2;
  };

  if (!cfg_path.empty()) {
    apply_config(cfg_path);
  } else {
    // try default config locations relative to working dir
    const vector<string> cands = {"pipeline_config.json", string("..")+"/pipeline_config.json", string("../..")+"/test/pipeline_config.json"};
    for (auto &c : cands) { apply_config(c); }
  }

  // now parse remaining CLI overrides
  kmer = parse_int_arg(args, "kmer", kmer);
  pair_info = parse_str_arg(args, "pair-info", pair_info);
  out_path = parse_str_arg(args, "out", out_path);

  // cfg already applied above

  ifstream pifs(pair_info);
  if (!pifs) {
    cerr << "[true_jaccard] pair_info not found: " << pair_info << "\n";
    return 2;
  }
  string header;
  getline(pifs, header); // skip header
  vector<PairInfo> pairs;
  string line;
  while (getline(pifs, line)) {
    if (line.empty()) continue;
    vector<string> cols; cols.reserve(5);
    size_t prev = 0; size_t pos;
    while ((pos = line.find('\t', prev)) != string::npos) { cols.emplace_back(line.substr(prev, pos - prev)); prev = pos + 1; }
    cols.emplace_back(line.substr(prev));
    if (cols.size() != 5) continue;
    PairInfo pi;
    pi.pair_id = stoi(cols[0]);
    pi.file1 = cols[1];
    pi.file2 = cols[2];
    pi.mutation_count = stoll(cols[3]);
    pi.genome_length = stoll(cols[4]);
    pairs.push_back(std::move(pi));
  }

  vector<JaccardOut> outs(pairs.size());
  const long long total = (long long)pairs.size();
  long long completed = 0;

  for (long long i = 0; i < total; ++i) {
    const auto &p = pairs[i];
    try {
      string s1 = read_fasta_concat(p.file1);
      string s2 = read_fasta_concat(p.file2);

      const size_t L1 = s1.size();
      const size_t L2 = s2.size();
      const int K = kmer;
      size_t n1 = kmer_count_possible(L1, K);
      size_t n2 = kmer_count_possible(L2, K);

      unordered_set<string> set1; set1.reserve(n1*1.3);
      unordered_set<string> set2; set2.reserve(n2*1.3);

      for (size_t j = 0; j < n1; ++j) set1.emplace(canonical_kmer(s1.substr(j, K)));
      for (size_t j = 0; j < n2; ++j) set2.emplace(canonical_kmer(s2.substr(j, K)));

      size_t c1 = set1.size();
      size_t c2 = set2.size();

      // intersection via iterating smaller set
      const unordered_set<string> *small = &set1, *large = &set2;
      if (set2.size() < set1.size()) { small = &set2; large = &set1; }
      size_t inter = 0;
      for (const auto &kmer : *small) if (large->find(kmer) != large->end()) ++inter;
      size_t uni = c1 + c2 - inter;
      double jac = (uni == 0) ? 0.0 : (double)inter / (double)uni;

      JaccardOut o;
      o.pair_id = p.pair_id;
      o.mutation_count = p.mutation_count;
      o.genome_length = p.genome_length;
      o.mutation_rate = (p.genome_length > 0) ? (double)p.mutation_count / (double)p.genome_length : 0.0;
      o.jaccard_true = jac;
      o.kmers1_count = c1;
      o.kmers2_count = c2;
      o.inter = inter;
      o.uni = uni;
      outs[i] = std::move(o);
      ++completed;
      if (completed % 20 == 0) {
        cerr << "[true_jaccard] progress " << completed << "/" << total << "\n";
      }
    } catch (const exception &e) {
      // leave default zeros; but report
      cerr << "[true_jaccard] error on pair " << p.pair_id << ": " << e.what() << "\n";
    }
  }

  ofstream ofs(out_path);
  if (!ofs) { cerr << "[true_jaccard] cannot open out: " << out_path << "\n"; return 3; }
  ofs << "pair_id\tmutation_count\tgenome_length\tmutation_rate\tjaccard_true\tkmers1_count\tkmers2_count\tintersection\tunion\n";
  for (const auto &o : outs) {
    ofs << o.pair_id << '\t' << o.mutation_count << '\t' << o.genome_length << '\t'
        << fixed << setprecision(8) << o.mutation_rate << '\t'
        << setprecision(10) << o.jaccard_true << '\t'
        << o.kmers1_count << '\t' << o.kmers2_count << '\t' << o.inter << '\t' << o.uni << "\n";
  }
  cerr << "[true_jaccard] done. pairs=" << outs.size() << ", kmer=" << kmer << "\n";
  return 0;
}
