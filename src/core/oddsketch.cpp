
#include "oddsketch/oddsketch_cli.hpp"
#include "third_party/xxhash.hpp"
#include "third_party/libpopcnt.hpp"

#include <algorithm> 
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>
#include <array>
#include <queue>
#include <cmath>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <sys/resource.h>
#include <unistd.h>

// スケッチサイズ（実行時可変、64 の倍数）
static size_t G_SKETCH_SIZE = 8192;

// k-mer 長さ（実行時に上書き可能）
static size_t KMER = 64;
// canonical k-mer を使うか（実行時に上書き可能）
static bool G_CANONICAL = true;
static bool G_CANONICAL_SET = false;
static bool G_KMER_SET = false;

// odd sketchの閾値（実行時に上書き可能）
static double J0 = 0.75;

// odd sketchのハッシュ値の個数（One Permutation Hashingと組み合わせるため）
// 実装上、ハッシュ数はスケッチビット数と同スケールで扱う

// ビット位置決定のモード
enum class PosMode { Value = 0, Mix = 1, Stripe = 2 };
static PosMode G_POS_MODE = PosMode::Value; // 既定: 互換目的で従来挙動（hv のみ）


// kmerをハッシュ化(xxhash)
u_int64_t hash_kmer(const std::string_view &kmer){
    return XXH64(kmer.data(), kmer.length(), 0);
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

// canonical k-mer のハッシュ（forward/revcomp のうち辞書順で小さい方）
static inline uint64_t hash_kmer_canonical(const char* kmer, size_t k, std::string &rcbuf) {
    rcbuf.resize(k);
    for (size_t i = 0; i < k; ++i) {
        rcbuf[k - 1 - i] = comp_base(kmer[i]);
    }
    std::string_view fwd(kmer, k);
    std::string_view rev(rcbuf.data(), k);
    const std::string_view &canon = (rev < fwd) ? rev : fwd;
    return XXH64(canon.data(), canon.size(), 0);
}

// One-Permutation Hashing による MinHash 実装
// 2の冪への丸め（以下の densify ではパワーオブツーが望ましい）
static inline size_t floor_pow2(size_t x) {
    if (x == 0) return 1;
    size_t p = 1;
    while ((p << 1) <= x) p <<= 1;
    return p;
}

// ---- 2-universal hash (multiply-shift; k は 2 の冪を推奨) ----
struct UnivHash {
    uint64_t a, b;  // a は奇数
    size_t k;       // バケット数（2^m）
    unsigned shift; // 64 - m

    explicit UnivHash(size_t k_, uint64_t seed = 0x9e3779b97f4a7c15ULL) : k(k_) {
        a = 0x2545F4914F6CDD1DULL ^ seed;
        if ((a & 1ULL) == 0) a ^= 1ULL; // 奇数化
        b = 0x9E3779B185EBCA87ULL + (seed << 1);
        unsigned m = 0; while ((1ULL << m) < k) ++m; shift = 64 - m;
    }
    inline size_t operator()(uint64_t x) const {
        uint64_t y = a * x + b;      // mod 2^64
        return (size_t)(y >> shift); // 上位ビットを抽出 → [0, k)
    }
};

// （i, attempt）を 64bit にパック
static inline uint64_t pack_pair(uint32_t i, uint32_t attempt) {
    return (uint64_t(i) << 32) ^ uint64_t(attempt);
}

// ---- Optimal densification (Shrivastava 2017, Algorithm 1) ----
// buckets: OPH の各バケット最小値（空は UINT64_MAX）
// 戻り値: 稠密化後の配列（すべて非空）
static std::vector<uint64_t> densify_optimal(const std::vector<uint64_t>& buckets) {
    const size_t k = buckets.size();
    std::vector<uint64_t> out(k);
    if (k == 0) return out;

    bool any = false;
    for (auto v : buckets) { if (v != UINT64_MAX) { any = true; break; } }
    if (!any) { std::fill(out.begin(), out.end(), 0ULL); return out; }

    UnivHash huniv(k);
    for (uint32_t i = 0; i < (uint32_t)k; ++i) {
        if (buckets[i] != UINT64_MAX) {
            out[i] = buckets[i];
        } else {
            uint32_t attempt = 1;
            while (true) {
                size_t next = huniv(pack_pair(i, attempt));
                if (buckets[next] != UINT64_MAX) { out[i] = buckets[next]; break; }
                ++attempt;
            }
        }
    }
    return out;
}

// 直近の OPH バケット数を保持（書き出し用）
static size_t G_LAST_NUM_BUCKETS = 0;

std::vector<uint64_t> get_minhash_one_permutation(const std::string &seq) {
    // バケット数は Odd Sketch 論文の推奨式に基づき設定:
    //   L = n / (4 * (1 - J0))
    // ここで n はスケッチのビット数（G_SKETCH_SIZE）とし、J0 はしきい値（デフォルト 0.75）。
    double denom = 4.0 * (1.0 - J0);
    size_t desired = (denom > 0.0)
        ? static_cast<size_t>(std::max(1.0, std::round(static_cast<double>(G_SKETCH_SIZE) / denom)))
        : G_SKETCH_SIZE;
    // 実装簡素化と均等性のため、2 の冪に丸める
    size_t NUM_BUCKETS = floor_pow2(desired);
    if (NUM_BUCKETS == 0) NUM_BUCKETS = 1;
    G_LAST_NUM_BUCKETS = NUM_BUCKETS;

    // 各バケットに最小ハッシュを保持（空は UINT64_MAX）
    std::vector<uint64_t> bucket_min_hash(NUM_BUCKETS, UINT64_MAX);

    // k-mer を走査して各バケットの min を更新
    size_t num_kmers = (seq.size() >= KMER) ? (seq.size() - KMER + 1) : 0;
    // バケットIDは上位ビットから決める（下位ビットは後段の pos=hv%M と相関させないため）
    unsigned m = 0; while ((1ULL << m) < NUM_BUCKETS) ++m; unsigned shift = 64 - m;
    if (G_CANONICAL) {
        std::string rcbuf;
        for (size_t i = 0; i < num_kmers; ++i) {
            const char* kptr = &seq[i];
            uint64_t hv = hash_kmer_canonical(kptr, KMER, rcbuf);
            size_t bucket_id = static_cast<size_t>(hv >> shift); // 上位ビット → [0, NUM_BUCKETS)
            if (hv < bucket_min_hash[bucket_id]) bucket_min_hash[bucket_id] = hv;
        }
    } else {
        for (size_t i = 0; i < num_kmers; ++i) {
            std::string_view kmer(&seq[i], KMER);
            uint64_t hv = hash_kmer(kmer);
            size_t bucket_id = static_cast<size_t>(hv >> shift); // 上位ビット → [0, NUM_BUCKETS)
            if (hv < bucket_min_hash[bucket_id]) bucket_min_hash[bucket_id] = hv;
        }
    }

    // 最適稠密化（optimal densification）
    auto dense = densify_optimal(bucket_min_hash);

    // 参考: 空率の期待値 s*exp(-n/s)
    // size_t s = NUM_BUCKETS; size_t n = num_kmers;
    // double expected_empty = s * std::exp(-(double)n / (double)s);
    // std::cerr << "[OPH] buckets=" << s << ", kmers=" << n << ", expected_empty≈" << expected_empty << "\n";

    // 各バケット値をそのまま返す（後段で mod G_SKETCH_SIZE してビット反転）
    return dense;
}

// ビット操作（64ビットワード配列）
static inline void flip_bit(std::vector<uint64_t> &words, size_t idx) {
    size_t w = idx / 64; size_t b = idx % 64; words[w] ^= (uint64_t(1) << b);
}

// FASTA ファイルからシーケンスを読みつつ Odd Sketch を構築
std::vector<uint64_t> make_odd_sketch_from_fasta(const std::string &fname) {
    std::ifstream ifs(fname);
    if (!ifs) throw std::runtime_error("Cannot open " + fname);

    std::string line;
    std::string seq;
    // FASTAの場合、">" で始まる行をスキップして塩基列を一続きに読む
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        seq += line;
    }

    // One Permutation Hashingを使用してMinHashを取得
    std::vector<uint64_t> minhash_values = get_minhash_one_permutation(seq);

    if (G_SKETCH_SIZE == 0 || (G_SKETCH_SIZE % 64) != 0) {
        throw std::runtime_error("Invalid sketch size (must be a positive multiple of 64)");
    }
    std::vector<uint64_t> words(G_SKETCH_SIZE / 64, 0);

    // One Permutation Hashingで得られた値（バケットごとの最小値）をスケッチに入れる
    auto map_pos = [&](uint32_t idx, uint64_t hv, size_t nbits, size_t kbuckets)->size_t{
        if (G_POS_MODE == PosMode::Value) {
            return static_cast<size_t>(hv % nbits);
        } else if (G_POS_MODE == PosMode::Mix) {
            // ビン番号と値を混ぜて位置決定（位置情報を保持しつつ衝突を分散）
            uint64_t buf[2]; buf[0] = hv; buf[1] = static_cast<uint64_t>(idx);
            uint64_t h = XXH64(reinterpret_cast<const void*>(buf), sizeof(buf), 0x9E3779B97F4A7C15ULL);
            return static_cast<size_t>(h % nbits);
        } else { // Stripe
            size_t stride = (kbuckets > 0) ? (nbits / kbuckets) : 0;
            if (stride >= 2) {
                size_t base = static_cast<size_t>(idx) * stride;
                size_t offset = static_cast<size_t>(hv % stride);
                size_t pos = base + offset;
                if (pos >= nbits) pos = pos % nbits; // 念のため
                return pos;
            } else {
                // ストライドが小さすぎる場合は Mix にフォールバック
                uint64_t buf[2]; buf[0] = hv; buf[1] = static_cast<uint64_t>(idx);
                uint64_t h = XXH64(reinterpret_cast<const void*>(buf), sizeof(buf), 0x9E3779B97F4A7C15ULL);
                return static_cast<size_t>(h % nbits);
            }
        }
    };

    const size_t kbuckets = minhash_values.size();
    for (size_t i = 0; i < kbuckets; i++) {
        uint64_t hv = minhash_values[i];
        size_t pos = map_pos(static_cast<uint32_t>(i), hv, G_SKETCH_SIZE, kbuckets);
        flip_bit(words, pos); // 反転 (odd‐sketch の核心)
    }
    return words;
}

// スケッチヘッダ（可搬性簡易のためバイナリ直書き）
struct SketchHeader {
    char magic[4];      // 'O','D','S','K'
    uint32_t version;   // 1
    uint64_t nbits;     // スケッチビット数
    uint64_t kbuckets;  // OPHバケット数（稠密化後）
    double j0;          // しきい値
};

// ヘッダ付きで書き出し
void write_sketch_binary(const std::vector<uint64_t>& words,
                         const std::string& outfname,
                         uint64_t nbits,
                         uint64_t kbuckets,
                         double j0val) {
    std::ofstream ofs(outfname, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot open output file: " + outfname);
    SketchHeader h{};
    h.magic[0]='O'; h.magic[1]='D'; h.magic[2]='S'; h.magic[3]='K';
    h.version = 1;
    h.nbits = nbits;
    h.kbuckets = kbuckets;
    h.j0 = j0val;
    ofs.write(reinterpret_cast<const char*>(&h), sizeof(h));
    if (!ofs) throw std::runtime_error("Error writing header: " + outfname);
    ofs.write(reinterpret_cast<const char*>(words.data()), words.size()*sizeof(uint64_t));
    if (!ofs) throw std::runtime_error("Error writing body: " + outfname);
}

// スケッチをワード単位で保持する型
struct Sketch {
    std::vector<uint64_t> words;  // 動的長
    uint64_t bit_size = 0;         // 総ビット数（ヘッダ）
    uint64_t k_buckets = 0;        // OPHバケット数（ヘッダ。0 の場合は未知）
    double j0 = 0.75;              // 参照用（ヘッダ）
};

// バイナリファイル (*.sketch) を読み込んで Sketch オブジェクトを返す
Sketch load_sketch(const std::string &fname) {
    Sketch s;
    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open sketch file: " + fname);
    // ヘッダ判定
    SketchHeader h{};
    ifs.read(reinterpret_cast<char*>(&h), sizeof(h));
    bool has_header = false;
    if (ifs && h.magic[0]=='O' && h.magic[1]=='D' && h.magic[2]=='S' && h.magic[3]=='K' && h.version>=1) {
        has_header = true;
    } else {
        // ヘッダなしレガシー: 先頭に戻して全体を本文として読む
        ifs.clear();
        ifs.seekg(0, std::ios::beg);
    }
    if (has_header) {
        // 残りは本文
        ifs.seekg(0, std::ios::end);
        std::streampos end = ifs.tellg();
        std::streampos body = end - static_cast<std::streampos>(sizeof(SketchHeader));
        if (body < 0 || (body % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
            throw std::runtime_error("Invalid sketch body size: " + fname);
        }
        size_t num_words = static_cast<size_t>(body / static_cast<std::streampos>(sizeof(uint64_t)));
        s.words.resize(num_words);
        ifs.seekg(sizeof(SketchHeader), std::ios::beg);
        ifs.read(reinterpret_cast<char*>(s.words.data()), num_words*sizeof(uint64_t));
        if (!ifs) throw std::runtime_error("Error reading sketch body: " + fname);
        s.bit_size = h.nbits ? h.nbits : static_cast<uint64_t>(s.words.size()*64);
        s.k_buckets = h.kbuckets;
        s.j0 = (h.j0>0 && h.j0<1) ? h.j0 : 0.75;
    } else {
        // レガシー: 全体を本文として読む
        ifs.seekg(0, std::ios::end);
        std::streampos end = ifs.tellg();
        ifs.seekg(0, std::ios::beg);
        if (end <= 0 || (end % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
            throw std::runtime_error("Invalid sketch file size: " + fname);
        }
        size_t num_words = static_cast<size_t>(end / static_cast<std::streampos>(sizeof(uint64_t)));
        s.words.resize(num_words);
        ifs.read(reinterpret_cast<char*>(s.words.data()), num_words*sizeof(uint64_t));
        if (!ifs) throw std::runtime_error("Error reading sketch file: " + fname);
        s.bit_size = static_cast<uint64_t>(num_words*64);
        s.k_buckets = 0; // 不明
        s.j0 = 0.75;
    }
    return s;
}

// ファイル名リストからすべてのスケッチを読み込む
std::vector<Sketch> load_all_sketches(const std::vector<std::string> &paths) {
    std::vector<Sketch> sketches;
    sketches.reserve(paths.size());
    for (const auto &p : paths) {
        sketches.push_back(load_sketch(p));
    }
    return sketches;
}

/**
 *  odd-sketch 同士の Jaccard 類似度推定値を返す
 *  popcnt の定義上、2*popcnt > SKETCH_SIZE のとき
 *  log の引数が負になるケースは Jaccard = 0 として扱う
 */
double jaccard_distance(const Sketch &a,
                        const Sketch &b) {
    if (a.words.size() != b.words.size()) {
        throw std::runtime_error("Sketch size mismatch");
    }
    const size_t NUM_WORDS = a.words.size();
    const size_t SKETCH_BITS = NUM_WORDS * 64;

    // 1) XOR -> POPCNT
    uint64_t popcnt = 0;
    for (size_t w = 0; w < NUM_WORDS; ++w) {
        popcnt += __builtin_popcountll(a.words[w] ^ b.words[w]);
    }

    // 2) 推定式に代入
    // OddSketch 論文の式に基づき係数を n/(4*k) とする（k は OPH のハッシュ数）
    // ここでは k を L = n / (4*(1-J0)) から再推定し、2 の冪に丸める。
    const double n = static_cast<double>(SKETCH_BITS);
    // 実際にスケッチ生成時に使われた k を優先（ヘッダ）。両方に無ければ式から推定。
    uint64_t k_header = (a.k_buckets>0 && b.k_buckets>0) ? std::min(a.k_buckets, b.k_buckets)
                        : (a.k_buckets>0 ? a.k_buckets : (b.k_buckets>0 ? b.k_buckets : 0));
    double k;
    if (k_header > 0) {
        k = static_cast<double>(k_header);
    } else {
        const double denom = 4.0 * (1.0 - J0);
        size_t k_est = (denom > 0.0)
            ? floor_pow2(static_cast<size_t>(std::max(1.0, std::round(n / denom))))
            : SKETCH_BITS;
        k = static_cast<double>(k_est);
    }

    // 数値安定化: x = 1 - 2*popcnt/n の対数をとる前にクリップ
    const double x = 1.0 - (2.0 * static_cast<double>(popcnt)) / n; // 1 - 2p
    const double eps = 1e-12; // 小さすぎる負値/ゼロを避ける
    double term = std::log(std::max(eps, x));
    double jacc = 1.0 + (n / (4.0 * k)) * term;

    // 範囲クリップ（理論上 [0,1]）
    if (jacc < 0.0) jacc = 0.0;
    if (jacc > 1.0) jacc = 1.0;
    return jacc;
}

static std::vector<std::string> read_paths_from_stdin() {
    std::vector<std::string> paths;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        paths.push_back(line);
    }
    return paths;
}

static std::vector<std::string> read_paths_from_list_file(const std::string &list_path) {
    std::ifstream ifs(list_path);
    if (!ifs) {
        throw std::runtime_error("Cannot open list file: " + list_path);
    }

    std::vector<std::string> paths;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        paths.push_back(line);
    }
    return paths;
}

static void run_dist_all_vs_all(const std::vector<std::string> &paths) {
    auto sketches = load_all_sketches(paths);
    if (!sketches.empty()) {
        const auto &s0 = sketches.front();
        uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(all-vs-all): nbits=" << nbits
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << J0 << "\n";
    }
    for (size_t i = 0; i < sketches.size(); ++i) {
        for (size_t j = i + 1; j < sketches.size(); ++j) {
            double d = jaccard_distance(sketches[i], sketches[j]);
            std::cout << paths[i] << '\t' << paths[j] << '\t' << d << "\n";
        }
    }
}

static void run_dist_bipartite(const std::vector<std::string> &qpaths,
                               const std::vector<std::string> &dbpaths) {
    auto qsketches = load_all_sketches(qpaths);
    auto dbsketches = load_all_sketches(dbpaths);
    if (!qsketches.empty()) {
        const auto &s0 = qsketches.front();
        uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(bipartite): nbits=" << nbits
                  << ", queries=" << qsketches.size()
                  << ", db=" << dbsketches.size()
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << J0 << "\n";
    }
    for (size_t qi = 0; qi < qsketches.size(); ++qi) {
        for (size_t di = 0; di < dbsketches.size(); ++di) {
            double d = jaccard_distance(qsketches[qi], dbsketches[di]);
            std::cout << qpaths[qi] << '\t' << dbpaths[di] << '\t' << d << "\n";
        }
    }
}

// 簡易オプションパーサ（--name=value 形式）
static void parse_options(int argc, char** argv, std::string &mode, std::string &qlist_path, std::string &dblist_path) {
    if (argc < 2) {
        return;
    }
    mode = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a.rfind("--", 0) != 0) continue;
        auto eq = a.find('=');
        std::string key = (eq == std::string::npos) ? a.substr(2) : a.substr(2, eq - 2);
        std::string val = (eq == std::string::npos) ? std::string() : a.substr(eq + 1);
        if (key == "qlist") {
            qlist_path = val;
        } else if (key == "dblist") {
            dblist_path = val;
        } else if (key == "kmer" || key == "kmerlen") {
            if (!val.empty()) {
                try {
                    KMER = static_cast<size_t>(std::stoul(val));
                    G_KMER_SET = true;
                } catch (...) {}
            }
        } else if (key == "sketch-size") {
            if (!val.empty()) {
                size_t req = 0; try { req = static_cast<size_t>(std::stoul(val)); } catch (...) {}
                if (req % 64 == 0 && req > 0) {
                    G_SKETCH_SIZE = req;
                } else {
                    std::cerr << "[oddsketch] warning: sketch-size must be a positive multiple of 64; got " << req << " (ignored)\n";
                }
            }
        } else if (key == "j0" || key == "j-threshold") {
            if (!val.empty()) {
                try {
                    double v = std::stod(val);
                    if (v > 0.0 && v < 1.0) {
                        J0 = v;
                    } else {
                        std::cerr << "[oddsketch] warning: j0 must be in (0,1); got " << v << " (ignored)\n";
                    }
                } catch (...) {}
            }
        } else if (key == "pos-mode") {
            // value|mix|stripe
            if (!val.empty()) {
                if (val == "value") G_POS_MODE = PosMode::Value;
                else if (val == "mix") G_POS_MODE = PosMode::Mix;
                else if (val == "stripe") G_POS_MODE = PosMode::Stripe;
            }
        } else if (key == "canonical") {
            if (val.empty()) {
                G_CANONICAL = true;
                G_CANONICAL_SET = true;
            } else {
                std::string v = val;
                std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return std::tolower(c); });
                if (v == "1" || v == "true" || v == "yes" || v == "on") {
                    G_CANONICAL = true;
                    G_CANONICAL_SET = true;
                } else if (v == "0" || v == "false" || v == "no" || v == "off") {
                    G_CANONICAL = false;
                    G_CANONICAL_SET = true;
                } else {
                    std::cerr << "[oddsketch] warning: canonical must be true/false/1/0; got " << val << " (ignored)\n";
                }
            }
        }
    }
}

// pipeline_config.json から oddsketch 設定を読み取る（任意）
static void load_config_from_file_if_exists() {
    const char* candidates[] = {
        "test/pipeline_config.json",
        "src/test/pipeline_config.json",
        "pipeline_config.json"
    };
    for (auto c : candidates) {
        std::ifstream ifs(c);
        if (!ifs) continue;
        std::string json((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        // 極簡易パース（"oddsketch" セクション内を優先）
        auto find_num = [&](const std::string& key, size_t pos)->bool{
            auto kpos = json.find(key, pos);
            if (kpos == std::string::npos) return false;
            auto colon = json.find(':', kpos);
            if (colon == std::string::npos) return false;
            auto beg = json.find_first_of("0123456789", colon);
            if (beg == std::string::npos) return false;
            auto end = json.find_first_not_of("0123456789", beg);
            try {
                if (!G_KMER_SET) KMER = static_cast<size_t>(std::stoul(json.substr(beg, end - beg)));
            } catch (...) {}
            return true;
        };
        auto find_bool = [&](const std::string& key, size_t pos)->bool{
            auto kpos = json.find(key, pos);
            if (kpos == std::string::npos) return false;
            auto colon = json.find(':', kpos);
            if (colon == std::string::npos) return false;
            auto beg = json.find_first_not_of(" \t\r\n", colon + 1);
            if (beg == std::string::npos) return false;
            if (json.compare(beg, 4, "true") == 0) {
                if (!G_CANONICAL_SET) G_CANONICAL = true;
                return true;
            }
            if (json.compare(beg, 5, "false") == 0) {
                if (!G_CANONICAL_SET) G_CANONICAL = false;
                return true;
            }
            if (json[beg] == '1' || json[beg] == '0') {
                if (!G_CANONICAL_SET) G_CANONICAL = (json[beg] == '1');
                return true;
            }
            return false;
        };
        bool updated = false;
        auto pos = json.find("\"oddsketch\"");
        size_t base = (pos == std::string::npos) ? 0 : pos;
        if (find_num("\"kmerlen\"", base)) updated = true;
        if (find_bool("\"canonical\"", base)) updated = true;
        if (!updated) {
            // fallback to top-level search
            if (find_num("\"kmerlen\"", 0)) updated = true;
            if (find_bool("\"canonical\"", 0)) updated = true;
        }
        if (updated) return;
    }
}

int oddsketch_cli_main(int argc, char** argv) {
    std::string mode;
    std::string qlist_path;
    std::string dblist_path;
    parse_options(argc, argv, mode, qlist_path, dblist_path);
    if (mode != "sketch" && mode != "dist") {
        std::cerr << "Usage: oddsketch {sketch|dist} [--kmer=N] [--sketch-size=M] [--canonical=0|1] [--qlist=queries.txt --dblist=db.txt]\n";
        return 1;
    }
    if ((qlist_path.empty() && !dblist_path.empty()) || (!qlist_path.empty() && dblist_path.empty())) {
        std::cerr << "Usage error: --qlist and --dblist must be specified together\n";
        return 1;
    }
    // 設定ファイルからのロード（CLI未指定時の補助）
    load_config_from_file_if_exists();

    if (mode == "sketch") {
        auto paths = read_paths_from_stdin();
        // デバッグ: スケッチサイズとバケット数、J0 を表示
        std::cerr << "[oddsketch] sketch: nbits=" << G_SKETCH_SIZE
                  << ", kbuckets=" << G_LAST_NUM_BUCKETS
                  << ", j0=" << std::fixed << std::setprecision(6) << J0
                  << ", canonical=" << (G_CANONICAL ? "true" : "false")
                  << ", pos-mode=" << (G_POS_MODE==PosMode::Value?"value":(G_POS_MODE==PosMode::Mix?"mix":"stripe"))
                  << "\n";

        for (auto &f : paths) {
            auto S = make_odd_sketch_from_fasta(f);
            // 出力ファイル名は f+".sketch" 固定
            write_sketch_binary(S, f + ".sketch", static_cast<uint64_t>(G_SKETCH_SIZE), static_cast<uint64_t>(G_LAST_NUM_BUCKETS), J0);
            std::cout << f + ".sketch" << "\n";
        }
    }
    else if (mode == "dist") {
        if (!qlist_path.empty() && !dblist_path.empty()) {
            auto qpaths = read_paths_from_list_file(qlist_path);
            auto dbpaths = read_paths_from_list_file(dblist_path);
            run_dist_bipartite(qpaths, dbpaths);
        } else {
            auto paths = read_paths_from_stdin();
            run_dist_all_vs_all(paths);
        }
    }
    
    return 0;
}
