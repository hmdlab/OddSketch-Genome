#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

struct K128 {
  uint64_t hi{0}, lo{0};
  bool operator<(K128 const& o) const noexcept { return (hi<o.hi) || (hi==o.hi && lo<o.lo); }
  bool operator==(K128 const& o) const noexcept { return hi==o.hi && lo==o.lo; }
};

static inline bool is_atcg(char c){
  switch(c){case 'A':case 'C':case 'G':case 'T':case 'a':case 'c':case 'g':case 't': return true;} return false;
}

static inline uint64_t two_bit(char c){
  switch(c){
    case 'A': case 'a': return 0u;
    case 'C': case 'c': return 1u;
    case 'G': case 'g': return 2u;
    case 'T': case 't': return 3u;
    default: return 0u;
  }
}

static inline char comp_base(char c){
  switch(c){
    case 'A': case 'a': return 'T';
    case 'C': case 'c': return 'G';
    case 'G': case 'g': return 'C';
    case 'T': case 't': return 'A';
    default: return 'N';
  }
}

static inline void push2(K128 &v, uint64_t two){
  uint64_t new_hi = (v.hi<<2) | (v.lo>>62);
  uint64_t new_lo = (v.lo<<2) | (two & 0x3u);
  v.hi = new_hi; v.lo = new_lo;
}

static K128 encode_kmer64(const string &s){
  K128 v{};
  for(char c: s){ push2(v, two_bit(c)); }
  return v;
}

static K128 canonical_kmer64(const string &s){
  // build revcomp string
  string rc; rc.resize(s.size());
  for(size_t i=0,n=s.size(); i<n; ++i) rc[n-1-i] = comp_base(s[i]);
  K128 a = encode_kmer64(s);
  K128 b = encode_kmer64(rc);
  return (b<a)?b:a;
}

static string read_fasta_concat(const string &path){
  ifstream ifs(path);
  if(!ifs) throw runtime_error("cannot open fasta: "+path);
  string line, seq; seq.reserve(1<<20);
  while(getline(ifs,line)){
    if(line.empty()) continue;
    if(line[0]=='>') continue;
    if(!line.empty() && line.back()=='\r') line.pop_back();
    seq += line;
  }
  return seq;
}

static size_t build_index_for_fasta(const string &fasta, int k){
  namespace fs = std::filesystem;
  string bin = fasta + ".k" + to_string(k) + ".bin";
  string idx = fasta + ".k" + to_string(k) + ".idx";
  // Reuse only if both index files exist AND are newer than the FASTA
  try{
    if(fs::exists(bin) && fs::exists(idx)){
      auto t_fast = fs::last_write_time(fasta);
      auto t_bin  = fs::last_write_time(bin);
      auto t_idx  = fs::last_write_time(idx);
      if(t_bin >= t_fast && t_idx >= t_fast){
        ifstream idxin(idx); size_t n=0; if(idxin>>n) return n; // reuse
      }
    }
  }catch(...){ /* fallthrough to rebuild */ }
  string seq = read_fasta_concat(fasta);
  vector<K128> km;
  if(seq.size() >= (size_t)k){
    km.reserve(seq.size()-k+1);
    for(size_t i=0;i+ (size_t)k <= seq.size(); ++i){
      const char *p = seq.data()+i;
      bool ok=true; for(int j=0;j<k;j++){ if(!is_atcg(p[j])){ ok=false; break; } }
      if(!ok) continue;
      km.push_back(canonical_kmer64(string(p, k)));
    }
    sort(km.begin(), km.end());
    km.erase(unique(km.begin(), km.end()), km.end());
  }
  ofstream bout(bin, ios::binary);
  for(const auto &v: km){
    uint64_t hi=v.hi, lo=v.lo;
    // big endian
    auto wr64=[&](uint64_t x){ unsigned char b[8]; for(int i=7;i>=0;--i){ b[7-i] = (unsigned char)((x>>(i*8)) & 0xFF); } bout.write((char*)b,8); };
    wr64(hi); wr64(lo);
  }
  ofstream idxout(idx); idxout<<km.size();
  return km.size();
}

static void load_index(const string &bin, vector<K128> &out){
  ifstream in(bin, ios::binary);
  out.clear();
  K128 v; unsigned char b[16];
  while(in.read((char*)b,16)){
    uint64_t hi=0, lo=0;
    for(int i=0;i<8;++i) hi = (hi<<8) | b[i];
    for(int i=8;i<16;++i) lo = (lo<<8) | b[i];
    v.hi=hi; v.lo=lo; out.push_back(v);
  }
}

//

static tuple<size_t,size_t,size_t> intersect_count(const vector<K128>&A, const vector<K128>&B){
  size_t i=0,j=0, inter=0; size_t na=A.size(), nb=B.size();
  while(i<na && j<nb){
    if(A[i]==B[j]){ ++inter; ++i; ++j; }
    else if(A[i]<B[j]) ++i; else ++j;
  }
  return {inter, na, nb};
}

static void usage(){
  cerr<<"true_index_pairs usage:\n"
         "  preprocess --list list.txt [--k 64]\n"
         "  pairs --qlist queries.list --dblist db.list --out-pairs true_pairs.tsv --out-nn true_nn.tsv [--k 64]\n";
}

static bool should_report_progress(size_t done, size_t total, size_t step_items = 50){
  if(total == 0 || done == 0) return false;
  return done >= total || done % step_items == 0;
}

struct IndexedGenome {
  string path;
  filesystem::path filename;
  vector<K128> kmers;
};

int main(int argc, char**argv){
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  if(argc<2){ usage(); return 1; }
  string mode = argv[1];
  int k=64; string list, qlist, dblist, outpairs, outnn;
  for(int i=2;i<argc;++i){ string a=argv[i];
    if(a=="--k" && i+1<argc){ k = stoi(argv[++i]); }
    else if(a=="--list" && i+1<argc){ list = argv[++i]; }
    else if(a=="--qlist" && i+1<argc){ qlist = argv[++i]; }
    else if(a=="--dblist" && i+1<argc){ dblist = argv[++i]; }
    else if(a=="--out-pairs" && i+1<argc){ outpairs = argv[++i]; }
    else if(a=="--out-nn" && i+1<argc){ outnn = argv[++i]; }
  }
  try{
    if(mode=="preprocess"){
      if(list.empty()){ usage(); return 2; }
      ifstream lf(list); if(!lf){ cerr<<"cannot open list: "<<list<<"\n"; return 3; }
      vector<string> paths;
      string p;
      while(getline(lf,p)){ if(!p.empty()) paths.push_back(p); }
      const size_t total = paths.size();
      size_t n=0;
      for(const auto &path : paths){
        build_index_for_fasta(path, k);
        ++n;
        if(should_report_progress(n, total)){
          cerr<<"[preprocess] done="<<n<<"/"<<total<<"\n";
        }
      }
      cerr<<"[preprocess] completed files="<<n<<"/"<<total<<" k="<<k<<"\n";
      return 0;
    }else if(mode=="pairs"){
      if(qlist.empty()||dblist.empty()||outpairs.empty()||outnn.empty()){ usage(); return 2; }
      ifstream qf(qlist), df(dblist);
      if(!qf||!df){ cerr<<"cannot open qlist/dblist\n"; return 3; }
      vector<string> qs, db; string s;
      while(getline(qf,s)) if(!s.empty()) qs.push_back(s);
      while(getline(df,s)) if(!s.empty()) db.push_back(s);
      // ensure indices exist
      for(auto &p: qs) build_index_for_fasta(p, k);
      for(auto &p: db) build_index_for_fasta(p, k);

      cerr<<"[pairs] loading DB indices: "<<db.size()<<" genomes\n";
      vector<IndexedGenome> db_indices;
      db_indices.reserve(db.size());
      for(size_t di=0; di<db.size(); ++di){
        IndexedGenome item;
        item.path = db[di];
        item.filename = filesystem::path(db[di]).filename();
        load_index(db[di]+".k"+to_string(k)+".bin", item.kmers);
        db_indices.push_back(std::move(item));
        if(should_report_progress(di+1, db.size())){
          cerr<<"[pairs] loaded DB indices "<<(di+1)<<"/"<<db.size()<<"\n";
        }
      }

      ofstream op(outpairs); op<<"query\tdb\tinter\tn1\tn2\tjaccard_true\n";
      unordered_map<string, pair<double,string>> best; best.reserve(qs.size()*2);
      const size_t total_pairs = qs.size() * db_indices.size();
      size_t processed_pairs = 0;
      cerr<<"[pairs] start comparisons: queries="<<qs.size()
          <<" db="<<db_indices.size()
          <<" total_pairs="<<total_pairs<<"\n";
      for(size_t qi=0; qi<qs.size(); ++qi){
        string q = qs[qi];
        vector<K128> qa; load_index(q+".k"+to_string(k)+".bin", qa);
        filesystem::path q_filename = filesystem::path(q).filename();
        for(size_t di=0; di<db_indices.size(); ++di){
          const auto &d = db_indices[di];
          if(d.filename == q_filename){ continue; }
          auto [inter, n1, n2] = intersect_count(qa, d.kmers);
          size_t uni = n1 + n2 - inter; double jac = (uni? (double)inter / (double)uni : 0.0);
          op<<q<<"\t"<<d.path<<"\t"<<inter<<"\t"<<n1<<"\t"<<n2<<"\t"<<fixed<<setprecision(10)<<jac<<"\n";
          auto it = best.find(q);
          if(it==best.end() || jac > it->second.first){ best[q] = {jac, d.path}; }
          ++processed_pairs;
        }
        if(should_report_progress(qi+1, qs.size())){
          cerr<<"[pairs] processed queries "<<(qi+1)<<"/"<<qs.size()
              <<", "<<processed_pairs<<"/"<<total_pairs<<" pairs\n";
        }
      }
      ofstream on(outnn); on<<"query\tnn_true\tjaccard_true\n";
      for(auto &kv: best){ on<<kv.first<<"\t"<<kv.second.second<<"\t"<<fixed<<setprecision(10)<<kv.second.first<<"\n"; }
      cerr<<"[pairs] wrote "<<outpairs<<" and "<<outnn<<"\n"; return 0;
    }else{
      usage(); return 1;
    }
  }catch(const exception &e){ cerr<<"error: "<<e.what()<<"\n"; return 10; }
}
