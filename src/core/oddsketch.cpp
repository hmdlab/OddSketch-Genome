#include "oddsketch/oddsketch_cli.hpp"
#include "third_party/xxhash.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr uint64_t kHashSeed = 0x9E3779B97F4A7C15ULL;
constexpr uint64_t kEmptyBucket = std::numeric_limits<uint64_t>::max();

enum class PosMode { Value = 0, Mix = 1, Stripe = 2 };

// CLI / 実験スクリプトから与えられる実行時パラメータ。
struct OddsketchOptions {
    size_t sketch_size = 8192; // 64 の倍数
    size_t kmer = 64;
    bool canonical = true;
    double j0 = 0.75;
    PosMode pos_mode = PosMode::Value;
    size_t threads = 1;
};

struct CliArgs {
    std::string mode;
    std::string qlist_path;
    std::string dblist_path;
    std::string pairlist_path;
    bool explicit_all_to_all = false;
    bool explicit_bipartite = false;
    OddsketchOptions options;
};

// OPH で得た稠密 MinHash と、そのとき実際に使ったバケット数。
struct MinhashResult {
    std::vector<uint64_t> values;
    size_t num_buckets = 0;
};

// sketch 本体に加え、後段の距離推定で使うメタ情報も返す。
struct SketchBuildResult {
    std::vector<uint64_t> words;
    size_t num_buckets = 0;
};

// on-disk 形式の簡易ヘッダ。旧形式との互換のため version を持つ。
struct SketchHeader {
    char magic[4];
    uint32_t version;
    uint64_t nbits;
    uint64_t kbuckets;
    double j0;
};

struct Sketch {
    std::vector<uint64_t> words;
    uint64_t bit_size = 0;
    uint64_t k_buckets = 0;
    double j0 = 0.75;
};

struct SketchPair {
    std::string left;
    std::string right;
};

struct PairlistGroup {
    std::string left;
    std::vector<size_t> pair_indices;
};

struct UnivHash {
    uint64_t a;
    uint64_t b;
    size_t k;
    unsigned shift;

    explicit UnivHash(size_t k_, uint64_t seed = kHashSeed) : a(0), b(0), k(k_), shift(64) {
        a = 0x2545F4914F6CDD1DULL ^ seed;
        if ((a & 1ULL) == 0) {
            a ^= 1ULL;
        }
        b = 0x9E3779B185EBCA87ULL + (seed << 1);
        unsigned m = 0;
        while ((1ULL << m) < k) {
            ++m;
        }
        shift = 64 - m;
    }

    size_t operator()(uint64_t x) const {
        return static_cast<size_t>((a * x + b) >> shift);
    }
};

uint64_t hash_kmer(std::string_view kmer) {
    return XXH64(kmer.data(), kmer.length(), 0);
}

char comp_base(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }
}

uint64_t hash_kmer_canonical(const char* kmer, size_t k, std::string& rcbuf) {
    rcbuf.resize(k);
    for (size_t i = 0; i < k; ++i) {
        rcbuf[k - 1 - i] = comp_base(kmer[i]);
    }

    std::string_view fwd(kmer, k);
    std::string_view rev(rcbuf.data(), k);
    const std::string_view canon = (rev < fwd) ? rev : fwd;
    return XXH64(canon.data(), canon.size(), 0);
}

size_t floor_pow2(size_t x) {
    if (x == 0) {
        return 1;
    }

    size_t p = 1;
    while ((p << 1) <= x) {
        p <<= 1;
    }
    return p;
}

size_t compute_num_buckets(const OddsketchOptions& options) {
    // OddSketch 論文の L = n / (4 * (1 - J0)) に基づいて OPH バケット数を決める。
    // 実装では均等性と扱いやすさのため 2 の冪に丸める。
    const double denom = 4.0 * (1.0 - options.j0);
    const size_t desired = (denom > 0.0)
        ? static_cast<size_t>(std::max(1.0, std::round(static_cast<double>(options.sketch_size) / denom)))
        : options.sketch_size;
    return std::max<size_t>(1, floor_pow2(desired));
}

size_t normalize_thread_count(size_t requested_threads, size_t task_count) {
    const size_t safe_threads = std::max<size_t>(1, requested_threads);
    return task_count == 0 ? 1 : std::min(safe_threads, task_count);
}

template <typename Func>
void run_in_parallel(size_t task_count, size_t requested_threads, Func&& func) {
    const size_t worker_count = normalize_thread_count(requested_threads, task_count);
    if (task_count == 0) {
        return;
    }
    if (worker_count == 1) {
        for (size_t i = 0; i < task_count; ++i) {
            func(i);
        }
        return;
    }

    std::atomic<size_t> next_index{0};
    std::atomic<bool> failed{false};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    auto worker = [&]() {
        while (!failed.load(std::memory_order_relaxed)) {
            const size_t idx = next_index.fetch_add(1, std::memory_order_relaxed);
            if (idx >= task_count) {
                break;
            }
            try {
                func(idx);
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) {
                        first_error = std::current_exception();
                    }
                }
                failed.store(true, std::memory_order_relaxed);
                break;
            }
        }
    };

    for (size_t i = 0; i < worker_count; ++i) {
        workers.emplace_back(worker);
    }
    for (auto& thread : workers) {
        thread.join();
    }

    if (first_error) {
        std::rethrow_exception(first_error);
    }
}

uint64_t pack_pair(uint32_t i, uint32_t attempt) {
    return (uint64_t(i) << 32) ^ uint64_t(attempt);
}

std::vector<uint64_t> densify_optimal(const std::vector<uint64_t>& buckets) {
    const size_t k = buckets.size();
    std::vector<uint64_t> out(k);
    if (k == 0) {
        return out;
    }

    const bool any = std::any_of(buckets.begin(), buckets.end(), [](uint64_t v) {
        return v != kEmptyBucket;
    });
    if (!any) {
        std::fill(out.begin(), out.end(), 0ULL);
        return out;
    }

    // 空バケットは Shrivastava 2017 の optimal densification に従って埋める。
    const UnivHash huniv(k);
    for (uint32_t i = 0; i < static_cast<uint32_t>(k); ++i) {
        if (buckets[i] != kEmptyBucket) {
            out[i] = buckets[i];
            continue;
        }

        uint32_t attempt = 1;
        while (true) {
            const size_t next = huniv(pack_pair(i, attempt));
            if (buckets[next] != kEmptyBucket) {
                out[i] = buckets[next];
                break;
            }
            ++attempt;
        }
    }
    return out;
}

MinhashResult get_minhash_one_permutation(const std::string& seq, const OddsketchOptions& options) {
    const size_t num_buckets = compute_num_buckets(options);
    std::vector<uint64_t> bucket_min_hash(num_buckets, kEmptyBucket);

    const size_t num_kmers = (seq.size() >= options.kmer) ? (seq.size() - options.kmer + 1) : 0;
    unsigned m = 0;
    while ((1ULL << m) < num_buckets) {
        ++m;
    }
    const unsigned shift = 64 - m;

    // OPH の各バケットに対して最小 hash 値だけを保持する。
    if (options.canonical) {
        std::string rcbuf;
        for (size_t i = 0; i < num_kmers; ++i) {
            const char* kptr = &seq[i];
            const uint64_t hv = hash_kmer_canonical(kptr, options.kmer, rcbuf);
            const size_t bucket_id = static_cast<size_t>(hv >> shift);
            if (hv < bucket_min_hash[bucket_id]) {
                bucket_min_hash[bucket_id] = hv;
            }
        }
    } else {
        for (size_t i = 0; i < num_kmers; ++i) {
            const std::string_view kmer(&seq[i], options.kmer);
            const uint64_t hv = hash_kmer(kmer);
            const size_t bucket_id = static_cast<size_t>(hv >> shift);
            if (hv < bucket_min_hash[bucket_id]) {
                bucket_min_hash[bucket_id] = hv;
            }
        }
    }

    return {densify_optimal(bucket_min_hash), num_buckets};
}

void flip_bit(std::vector<uint64_t>& words, size_t idx) {
    const size_t w = idx / 64;
    const size_t b = idx % 64;
    words[w] ^= (uint64_t(1) << b);
}

size_t map_position(uint32_t idx, uint64_t hv, size_t nbits, size_t kbuckets, PosMode pos_mode) {
    if (pos_mode == PosMode::Value) {
        return static_cast<size_t>(hv % nbits);
    }

    if (pos_mode == PosMode::Stripe) {
        // bucket ごとに担当領域を割り当て、局所的に bit を立てる。
        const size_t stride = (kbuckets > 0) ? (nbits / kbuckets) : 0;
        if (stride >= 2) {
            const size_t base = static_cast<size_t>(idx) * stride;
            const size_t offset = static_cast<size_t>(hv % stride);
            const size_t pos = base + offset;
            return (pos < nbits) ? pos : (pos % nbits);
        }
    }

    // Mix は bucket index と hash 値を再混合して位置相関を弱める。
    uint64_t buf[2];
    buf[0] = hv;
    buf[1] = static_cast<uint64_t>(idx);
    const uint64_t mixed = XXH64(reinterpret_cast<const void*>(buf), sizeof(buf), kHashSeed);
    return static_cast<size_t>(mixed % nbits);
}

SketchBuildResult make_odd_sketch_from_fasta(const std::string& fname, const OddsketchOptions& options) {
    std::ifstream ifs(fname);
    if (!ifs) {
        throw std::runtime_error("Cannot open " + fname);
    }

    std::string line;
    std::string seq;
    // 現在は FASTA を単純に 1 本の塩基列として連結して扱う。
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '>') {
            continue;
        }
        seq += line;
    }

    if (options.sketch_size == 0 || (options.sketch_size % 64) != 0) {
        throw std::runtime_error("Invalid sketch size (must be a positive multiple of 64)");
    }

    const MinhashResult minhash = get_minhash_one_permutation(seq, options);
    std::vector<uint64_t> words(options.sketch_size / 64, 0);
    for (size_t i = 0; i < minhash.values.size(); ++i) {
        const uint64_t hv = minhash.values[i];
        const size_t pos = map_position(
            static_cast<uint32_t>(i),
            hv,
            options.sketch_size,
            minhash.values.size(),
            options.pos_mode
        );
        // OddSketch の本体は「該当 bit を反転する」操作。
        flip_bit(words, pos);
    }

    return {std::move(words), minhash.num_buckets};
}

bool has_valid_sketch_header(const SketchHeader& h) {
    return h.magic[0] == 'O' &&
           h.magic[1] == 'D' &&
           h.magic[2] == 'S' &&
           h.magic[3] == 'K' &&
           h.version >= 1;
}

Sketch load_legacy_sketch(std::ifstream& ifs, const std::string& fname) {
    // 旧形式はヘッダ無しで word 配列だけが並んでいる前提で読む。
    ifs.clear();
    ifs.seekg(0, std::ios::end);
    const std::streampos end = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if (end <= 0 || (end % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
        throw std::runtime_error("Invalid sketch file size: " + fname);
    }

    const size_t num_words = static_cast<size_t>(end / static_cast<std::streampos>(sizeof(uint64_t)));
    Sketch s;
    s.words.resize(num_words);
    ifs.read(reinterpret_cast<char*>(s.words.data()), num_words * sizeof(uint64_t));
    if (!ifs) {
        throw std::runtime_error("Error reading sketch file: " + fname);
    }
    s.bit_size = static_cast<uint64_t>(num_words * 64);
    return s;
}

void write_sketch_binary(const std::vector<uint64_t>& words,
                         const std::string& outfname,
                         uint64_t nbits,
                         uint64_t kbuckets,
                         double j0val) {
    std::ofstream ofs(outfname, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open output file: " + outfname);
    }

    SketchHeader h{};
    h.magic[0] = 'O';
    h.magic[1] = 'D';
    h.magic[2] = 'S';
    h.magic[3] = 'K';
    h.version = 1;
    h.nbits = nbits;
    h.kbuckets = kbuckets;
    h.j0 = j0val;

    ofs.write(reinterpret_cast<const char*>(&h), sizeof(h));
    if (!ofs) {
        throw std::runtime_error("Error writing header: " + outfname);
    }

    ofs.write(reinterpret_cast<const char*>(words.data()), words.size() * sizeof(uint64_t));
    if (!ofs) {
        throw std::runtime_error("Error writing body: " + outfname);
    }
}

Sketch load_sketch(const std::string& fname) {
    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Cannot open sketch file: " + fname);
    }

    SketchHeader h{};
    ifs.read(reinterpret_cast<char*>(&h), sizeof(h));
    // 先頭が新形式ヘッダでなければ、旧形式 reader にフォールバックする。
    if (!ifs || !has_valid_sketch_header(h)) {
        return load_legacy_sketch(ifs, fname);
    }

    ifs.seekg(0, std::ios::end);
    const std::streampos end = ifs.tellg();
    const std::streampos body = end - static_cast<std::streampos>(sizeof(SketchHeader));
    if (body < 0 || (body % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
        throw std::runtime_error("Invalid sketch body size: " + fname);
    }

    const size_t num_words = static_cast<size_t>(body / static_cast<std::streampos>(sizeof(uint64_t)));
    Sketch s;
    s.words.resize(num_words);
    ifs.seekg(sizeof(SketchHeader), std::ios::beg);
    ifs.read(reinterpret_cast<char*>(s.words.data()), num_words * sizeof(uint64_t));
    if (!ifs) {
        throw std::runtime_error("Error reading sketch body: " + fname);
    }

    s.bit_size = h.nbits ? h.nbits : static_cast<uint64_t>(s.words.size() * 64);
    s.k_buckets = h.kbuckets;
    s.j0 = (h.j0 > 0.0 && h.j0 < 1.0) ? h.j0 : 0.75;
    return s;
}

std::vector<Sketch> load_all_sketches(const std::vector<std::string>& paths) {
    std::vector<Sketch> sketches;
    sketches.reserve(paths.size());
    for (const auto& path : paths) {
        sketches.push_back(load_sketch(path));
    }
    return sketches;
}

double jaccard_distance(const Sketch& a, const Sketch& b, const OddsketchOptions& options) {
    if (a.words.size() != b.words.size()) {
        throw std::runtime_error("Sketch size mismatch");
    }

    // XOR した bit 差分数から論文の推定式を評価する。
    uint64_t popcnt = 0;
    for (size_t w = 0; w < a.words.size(); ++w) {
        popcnt += __builtin_popcountll(a.words[w] ^ b.words[w]);
    }

    const double n = static_cast<double>(a.words.size() * 64);
    const uint64_t k_header = (a.k_buckets > 0 && b.k_buckets > 0)
        ? std::min(a.k_buckets, b.k_buckets)
        : (a.k_buckets > 0 ? a.k_buckets : b.k_buckets);

    double k = 0.0;
    // 生成時の k がヘッダにあればそれを優先し、無ければ j0 から再推定する。
    if (k_header > 0) {
        k = static_cast<double>(k_header);
    } else {
        const double denom = 4.0 * (1.0 - options.j0);
        const size_t k_est = (denom > 0.0)
            ? floor_pow2(static_cast<size_t>(std::max(1.0, std::round(n / denom))))
            : static_cast<size_t>(n);
        k = static_cast<double>(k_est);
    }

    const double x = 1.0 - (2.0 * static_cast<double>(popcnt)) / n;
    const double term = std::log(std::max(1e-12, x));
    double jacc = 1.0 + (n / (4.0 * k)) * term;
    jacc = std::clamp(jacc, 0.0, 1.0);
    return jacc;
}

std::vector<std::string> read_nonempty_lines(std::istream& input) {
    std::vector<std::string> paths;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) {
            paths.push_back(line);
        }
    }
    return paths;
}

std::vector<std::string> read_paths_from_stdin() {
    return read_nonempty_lines(std::cin);
}

std::vector<std::string> read_paths_from_list_file(const std::string& list_path) {
    std::ifstream ifs(list_path);
    if (!ifs) {
        throw std::runtime_error("Cannot open list file: " + list_path);
    }
    return read_nonempty_lines(ifs);
}

std::vector<SketchPair> read_pairs_from_pairlist_file(const std::string& pairlist_path) {
    std::ifstream ifs(pairlist_path);
    if (!ifs) {
        throw std::runtime_error("Cannot open pairlist file: " + pairlist_path);
    }

    std::vector<SketchPair> pairs;
    std::string line;
    size_t lineno = 0;
    while (std::getline(ifs, line)) {
        ++lineno;
        if (line.empty()) {
            continue;
        }

        const auto tab = line.find('\t');
        if (tab == std::string::npos) {
            throw std::runtime_error("Invalid pairlist line " + std::to_string(lineno) + ": expected two tab-separated paths");
        }

        const std::string left = line.substr(0, tab);
        const std::string right = line.substr(tab + 1);
        if (left.empty() || right.empty() || right.find('\t') != std::string::npos) {
            throw std::runtime_error("Invalid pairlist line " + std::to_string(lineno) + ": expected exactly two tab-separated paths");
        }
        pairs.push_back({left, right});
    }
    return pairs;
}

std::vector<PairlistGroup> group_pairs_by_left(const std::vector<SketchPair>& pairs) {
    std::unordered_map<std::string, size_t> group_index_by_left;
    std::vector<PairlistGroup> groups;
    groups.reserve(pairs.size());

    for (size_t i = 0; i < pairs.size(); ++i) {
        const auto& left = pairs[i].left;
        const auto found = group_index_by_left.find(left);
        if (found == group_index_by_left.end()) {
            const size_t new_group_index = groups.size();
            group_index_by_left.emplace(left, new_group_index);
            groups.push_back({left, {i}});
        } else {
            groups[found->second].pair_indices.push_back(i);
        }
    }
    return groups;
}

const char* pos_mode_name(PosMode pos_mode) {
    switch (pos_mode) {
        case PosMode::Value:
            return "value";
        case PosMode::Mix:
            return "mix";
        case PosMode::Stripe:
            return "stripe";
    }
    return "value";
}

void run_dist_all_vs_all(const std::vector<std::string>& paths, const OddsketchOptions& options) {
    const auto sketches = load_all_sketches(paths);
    if (!sketches.empty()) {
        const auto& s0 = sketches.front();
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(all-to-all): nbits=" << nbits
                  << ", threads=" << normalize_thread_count(options.threads, sketches.size())
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::mutex output_mutex;
    run_in_parallel(sketches.size(), options.threads, [&](size_t i) {
        std::ostringstream chunk;
        for (size_t j = i + 1; j < sketches.size(); ++j) {
            const double d = jaccard_distance(sketches[i], sketches[j], options);
            chunk << paths[i] << '\t' << paths[j] << '\t' << d << "\n";
        }
        const std::string text = chunk.str();
        if (!text.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << text;
        }
    });
}

void run_dist_bipartite(const std::vector<std::string>& qpaths,
                        const std::vector<std::string>& dbpaths,
                        const OddsketchOptions& options) {
    const auto qsketches = load_all_sketches(qpaths);
    const auto dbsketches = load_all_sketches(dbpaths);
    if (!qsketches.empty()) {
        const auto& s0 = qsketches.front();
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(bipartite): nbits=" << nbits
                  << ", queries=" << qsketches.size()
                  << ", db=" << dbsketches.size()
                  << ", threads=" << normalize_thread_count(options.threads, qsketches.size())
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::mutex output_mutex;
    run_in_parallel(qsketches.size(), options.threads, [&](size_t qi) {
        std::ostringstream chunk;
        for (size_t di = 0; di < dbsketches.size(); ++di) {
            const double d = jaccard_distance(qsketches[qi], dbsketches[di], options);
            chunk << qpaths[qi] << '\t' << dbpaths[di] << '\t' << d << "\n";
        }
        const std::string text = chunk.str();
        if (!text.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << text;
        }
    });
}

void run_dist_pairlist(const std::vector<SketchPair>& pairs, const OddsketchOptions& options) {
    std::unordered_map<std::string, Sketch> cache;
    cache.reserve(pairs.size() * 2);

    for (const auto& pair : pairs) {
        if (cache.find(pair.left) == cache.end()) {
            cache.emplace(pair.left, load_sketch(pair.left));
        }
        if (cache.find(pair.right) == cache.end()) {
            cache.emplace(pair.right, load_sketch(pair.right));
        }
    }

    const auto groups = group_pairs_by_left(pairs);

    if (!pairs.empty()) {
        const auto& s0 = cache.at(pairs.front().left);
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(pairlist): nbits=" << nbits
                  << ", pairs=" << pairs.size()
                  << ", query_groups=" << groups.size()
                  << ", threads=" << normalize_thread_count(options.threads, groups.size())
                  << ", unique_sketches=" << cache.size()
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::vector<double> results(pairs.size(), 0.0);
    run_in_parallel(groups.size(), options.threads, [&](size_t group_index) {
        const auto& group = groups[group_index];
        const auto& left = cache.at(group.left);
        for (const size_t pair_index : group.pair_indices) {
            results[pair_index] = jaccard_distance(left, cache.at(pairs[pair_index].right), options);
        }
    });

    for (size_t i = 0; i < pairs.size(); ++i) {
        std::cout << pairs[i].left << '\t' << pairs[i].right << '\t' << results[i] << "\n";
    }
}

void set_kmer(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }
    options.kmer = static_cast<size_t>(std::stoul(value));
}

void set_sketch_size(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const size_t sketch_size = static_cast<size_t>(std::stoul(value));
    if (sketch_size == 0 || (sketch_size % 64) != 0) {
        throw std::runtime_error("sketch-size must be a positive multiple of 64");
    }
    options.sketch_size = sketch_size;
}

void set_j0(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const double j0 = std::stod(value);
    if (j0 <= 0.0 || j0 >= 1.0) {
        throw std::runtime_error("j0 must be in (0,1)");
    }
    options.j0 = j0;
}

void set_pos_mode(OddsketchOptions& options, const std::string& value) {
    if (value == "value") {
        options.pos_mode = PosMode::Value;
    } else if (value == "mix") {
        options.pos_mode = PosMode::Mix;
    } else if (value == "stripe") {
        options.pos_mode = PosMode::Stripe;
    } else {
        throw std::runtime_error("pos-mode must be one of: value, mix, stripe");
    }
}

void set_canonical(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        options.canonical = true;
        return;
    }

    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on") {
        options.canonical = true;
    } else if (normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off") {
        options.canonical = false;
    } else {
        throw std::runtime_error("canonical must be true/false/1/0");
    }
}

void set_threads(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const size_t threads = static_cast<size_t>(std::stoul(value));
    if (threads == 0) {
        throw std::runtime_error("threads must be positive");
    }
    options.threads = threads;
}

bool option_requires_value(const std::string& key) {
    return key == "qlist" ||
           key == "dblist" ||
           key == "pairlist" ||
           key == "kmer" ||
           key == "kmerlen" ||
           key == "sketch-size" ||
           key == "j0" ||
           key == "j-threshold" ||
           key == "pos-mode" ||
           key == "canonical" ||
           key == "threads";
}

void require_value(const std::string& key, const std::string& value) {
    if (value.empty()) {
        throw std::runtime_error("--" + key + " requires a value");
    }
}

CliArgs parse_options(int argc, char** argv) {
    CliArgs args;
    if (argc < 2) {
        return args;
    }

    args.mode = argv[1];
    if (args.mode == "--help" || args.mode == "-h") {
        return args;
    }

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--", 0) != 0) {
            continue;
        }

        const auto eq = arg.find('=');
        const std::string key = (eq == std::string::npos) ? arg.substr(2) : arg.substr(2, eq - 2);
        std::string value = (eq == std::string::npos) ? std::string() : arg.substr(eq + 1);

        if (eq == std::string::npos && option_requires_value(key)) {
            if (i + 1 < argc) {
                const std::string next = argv[i + 1];
                if (next.rfind("--", 0) != 0) {
                    value = next;
                    ++i;
                }
            }
            require_value(key, value);
        }

        if (key == "qlist") {
            require_value(key, value);
            args.qlist_path = value;
        } else if (key == "dblist") {
            require_value(key, value);
            args.dblist_path = value;
        } else if (key == "pairlist") {
            require_value(key, value);
            args.pairlist_path = value;
        } else if (key == "all-to-all" || key == "alltoall") {
            args.explicit_all_to_all = true;
        } else if (key == "bipartite") {
            args.explicit_bipartite = true;
        } else if (key == "kmer" || key == "kmerlen") {
            set_kmer(args.options, value);
        } else if (key == "sketch-size") {
            set_sketch_size(args.options, value);
        } else if (key == "j0" || key == "j-threshold") {
            set_j0(args.options, value);
        } else if (key == "pos-mode") {
            set_pos_mode(args.options, value);
        } else if (key == "canonical") {
            set_canonical(args.options, value);
        } else if (key == "threads") {
            set_threads(args.options, value);
        } else {
            throw std::runtime_error("unknown option: --" + key);
        }
    }

    return args;
}

}  // namespace

void print_usage() {
    std::cerr
        << "Usage:\n"
        << "  oddsketch sketch [options] < paths.list\n"
        << "  oddsketch dist --all-to-all [options] < sketches.list\n"
        << "  oddsketch dist --bipartite --qlist queries.list --dblist db.list [options]\n"
        << "  oddsketch dist --pairlist pairs.tsv [options]\n"
        << "\n"
        << "Options:\n"
        << "  --kmer=N, --kmerlen=N\n"
        << "  --sketch-size=M\n"
        << "  --canonical=0|1\n"
        << "  --pos-mode=value|mix|stripe\n"
        << "  --j0=F\n"
        << "  --threads=N\n"
        << "\n"
        << "Compatibility:\n"
        << "  oddsketch dist < sketches.list still runs all-to-all.\n"
        << "  --name=value and --name value are both accepted for value options.\n";
}

int oddsketch_cli_main(int argc, char** argv) {
    CliArgs args;
    try {
        args = parse_options(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Usage error: " << e.what() << "\n";
        return 1;
    }

    if (args.mode == "--help" || args.mode == "-h") {
        print_usage();
        return 0;
    }

    if (args.mode != "sketch" && args.mode != "dist") {
        print_usage();
        return 1;
    }
    if ((args.qlist_path.empty() && !args.dblist_path.empty()) ||
        (!args.qlist_path.empty() && args.dblist_path.empty())) {
        std::cerr << "Usage error: --qlist and --dblist must be specified together\n";
        return 1;
    }

    const bool pairlist_mode = !args.pairlist_path.empty();
    const bool bipartite_mode = args.explicit_bipartite || !args.qlist_path.empty();
    const bool all_to_all_mode = args.explicit_all_to_all;
    const int dist_mode_count = (pairlist_mode ? 1 : 0) + (bipartite_mode ? 1 : 0) + (all_to_all_mode ? 1 : 0);

    if (args.mode == "dist" && dist_mode_count > 1) {
        std::cerr << "Usage error: choose only one dist mode: --all-to-all, --bipartite, or --pairlist\n";
        return 1;
    }
    if (args.explicit_bipartite && (args.qlist_path.empty() || args.dblist_path.empty())) {
        std::cerr << "Usage error: --bipartite requires --qlist and --dblist\n";
        return 1;
    }

    try {
        if (args.mode == "sketch") {
            const auto paths = read_paths_from_stdin();
            std::cerr << "[oddsketch] sketch: nbits=" << args.options.sketch_size
                      << ", kbuckets=" << compute_num_buckets(args.options)
                      << ", threads=" << normalize_thread_count(args.options.threads, paths.size())
                      << ", j0=" << std::fixed << std::setprecision(6) << args.options.j0
                      << ", canonical=" << (args.options.canonical ? "true" : "false")
                      << ", pos-mode=" << pos_mode_name(args.options.pos_mode)
                      << "\n";

            // sketch は 1 行 1 FASTA パスを stdin から受け取り、同名 .sketch を出力する。
            std::vector<std::string> output_paths(paths.size());
            run_in_parallel(paths.size(), args.options.threads, [&](size_t i) {
                const SketchBuildResult sketch = make_odd_sketch_from_fasta(paths[i], args.options);
                write_sketch_binary(
                    sketch.words,
                    paths[i] + ".sketch",
                    static_cast<uint64_t>(args.options.sketch_size),
                    static_cast<uint64_t>(sketch.num_buckets),
                    args.options.j0
                );
                output_paths[i] = paths[i] + ".sketch";
            });

            for (const auto& output_path : output_paths) {
                std::cout << output_path << "\n";
            }
        } else if (bipartite_mode) {
            run_dist_bipartite(
                read_paths_from_list_file(args.qlist_path),
                read_paths_from_list_file(args.dblist_path),
                args.options
            );
        } else if (pairlist_mode) {
            run_dist_pairlist(read_pairs_from_pairlist_file(args.pairlist_path), args.options);
        } else {
            run_dist_all_vs_all(read_paths_from_stdin(), args.options);
        }
    } catch (const std::exception& e) {
        std::cerr << "[oddsketch] error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
