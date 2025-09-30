
#include "xxhash.hpp"
#include "xxhash_header_only.hpp"
#include "libpopcnt.hpp"

#include <algorithm> 
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

// odd sketchの閾値
constexpr double J0 = 0.75;

// odd sketchのハッシュ関数の数（One Permutation Hashingと組み合わせるため）
// 実装上、ハッシュ数はスケッチビット数と同スケールで扱う


// kmerをハッシュ化(xxhash)
u_int64_t hash_kmer(const std::string_view &kmer){
    return XXH64(kmer.data(), kmer.length(), 0);
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

std::vector<uint64_t> get_minhash_one_permutation(const std::string &seq) {
    // バケット数は G_SKETCH_SIZE 以下の 2 の冪（推奨）
    size_t desired = G_SKETCH_SIZE;
    size_t NUM_BUCKETS = floor_pow2(desired);
    if (NUM_BUCKETS == 0) NUM_BUCKETS = 1;

    // 各バケットに最小ハッシュを保持（空は UINT64_MAX）
    std::vector<uint64_t> bucket_min_hash(NUM_BUCKETS, UINT64_MAX);

    // k-mer を走査して各バケットの min を更新
    size_t num_kmers = (seq.size() >= KMER) ? (seq.size() - KMER + 1) : 0;
    // バケットIDは上位ビットから決める（下位ビットは後段の pos=hv%M と相関させないため）
    unsigned m = 0; while ((1ULL << m) < NUM_BUCKETS) ++m; unsigned shift = 64 - m;
    for (size_t i = 0; i < num_kmers; ++i) {
        std::string_view kmer(&seq[i], KMER);
        uint64_t hv = hash_kmer(kmer);
        size_t bucket_id = static_cast<size_t>(hv >> shift); // 上位ビット → [0, NUM_BUCKETS)
        if (hv < bucket_min_hash[bucket_id]) bucket_min_hash[bucket_id] = hv;
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

    // One Permutation Hashingで得られた値（バケットごとの最小値）をスケッチに入れる（直写像のベースライン）
    for (size_t i = 0; i < minhash_values.size(); i++) {
        uint64_t hv = minhash_values[i];
        size_t pos = static_cast<size_t>(hv % G_SKETCH_SIZE);  // 直写像: hv の下位ビットで位置決定
        flip_bit(words, pos);                                   // 反転 (odd‐sketch の核心)
    }
    return words;
}

// 生の uint64_t 配列をバイナリ書き出し
void write_sketch_binary(const std::vector<uint64_t>& words,
                         const std::string& outfname) {
    // バイナリモードでファイルを開いて write
    std::ofstream ofs(outfname, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open output file: " + outfname);
    }
    ofs.write(reinterpret_cast<const char*>(words.data()),
              words.size() * sizeof(uint64_t));
    if (!ofs) {
        throw std::runtime_error("Error writing to file: " + outfname);
    }
}

// スケッチをワード単位で保持する型
struct Sketch {
    std::vector<uint64_t> words;  // 動的長
};

// バイナリファイル (*.sketch) を読み込んで Sketch オブジェクトを返す
Sketch load_sketch(const std::string &fname) {
    Sketch s;
    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open sketch file: " + fname);
    ifs.seekg(0, std::ios::end);
    std::streampos end = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    if (end <= 0 || (end % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
        throw std::runtime_error("Invalid sketch file size: " + fname);
    }
    size_t num_words = static_cast<size_t>(end / static_cast<std::streampos>(sizeof(uint64_t)));
    s.words.resize(num_words);
    ifs.read(reinterpret_cast<char*>(s.words.data()), num_words * sizeof(uint64_t));
    if (!ifs) throw std::runtime_error("Error reading sketch file: " + fname);
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

    // 2*popcnt > sketch_size のときは 0 にクリップ
    if (2 * popcnt > SKETCH_BITS) {
        return 0.0;
    }

    // 2) 式に代入
    double d_pop    = static_cast<double>(popcnt);
    double d_size   = static_cast<double>(SKETCH_BITS);
    double d_hashes = static_cast<double>(SKETCH_BITS); // 近似: ハッシュ数 ≒ ビット数

    double ratio = (2.0 * d_pop) / d_size;            // = 2*popcnt/sketch_size
    double term  = std::log(1.0 - ratio);             // ≤ 0
    double jacc  = 1.0 + (d_size / (4.0 * d_hashes)) * term;
    return jacc;
}

// 簡易オプションパーサ（--name=value 形式）
static void parse_options(int argc, char** argv, std::string &mode) {
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
        if (key == "kmer" || key == "kmerlen") {
            if (!val.empty()) {
                try { KMER = static_cast<size_t>(std::stoul(val)); } catch (...) {}
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
        }
    }
}

// pipeline_config.json から kmer を読み取る（任意）
static void load_kmer_from_config_if_exists() {
    const char* candidates[] = {
        "test/pipeline_config.json",
        "src/test/pipeline_config.json",
        "pipeline_config.json"
    };
    for (auto c : candidates) {
        std::ifstream ifs(c);
        if (!ifs) continue;
        std::string json((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        // 極簡易パース（"oddsketch" セクション内の kmerlen または トップレベル kmerlen を探す）
        auto pos = json.find("\"oddsketch\"");
        auto find_num = [&](const std::string& key)->bool{
            auto kpos = json.find(key, pos == std::string::npos ? 0 : pos);
            if (kpos == std::string::npos) return false;
            auto colon = json.find(':', kpos);
            if (colon == std::string::npos) return false;
            auto beg = json.find_first_of("0123456789", colon);
            if (beg == std::string::npos) return false;
            auto end = json.find_first_not_of("0123456789", beg);
            try { KMER = static_cast<size_t>(std::stoul(json.substr(beg, end - beg))); } catch (...) {}
            return true;
        };
        if (find_num("\"kmerlen\"")) return;
        // fallback to top-level search
        pos = 0; if (find_num("\"kmerlen\"")) return;
    }
}

int main(int argc, char** argv) {
    std::string mode;
    parse_options(argc, argv, mode);
    if (mode != "sketch" && mode != "dist") {
        std::cerr << "Usage: oddsketch {sketch|dist} [--kmer=N] [--sketch-size=M]\n";
        return 1;
    }
    // 設定ファイルからのkmerロード（CLI未指定時の補助）
    if (KMER == 64) {
        load_kmer_from_config_if_exists();
    }

    // 1行ずつ受け取ってベクタに貯める
    std::vector<std::string> paths;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        paths.push_back(line);
    }

    if (mode == "sketch") {

        for (auto &f : paths) {
            auto S = make_odd_sketch_from_fasta(f);
            // 出力ファイル名は f+".sketch" 固定
            write_sketch_binary(S, f + ".sketch");
            std::cout << f + ".sketch" << "\n";
        }
    }
    else if (mode == "dist") {
        // すべての .sketch を読み込んでメモリに展開
        auto sketches = load_all_sketches(paths);
        // 二重ループで距離計算して stdout に TSV 出力
        for (size_t i = 0; i < sketches.size(); ++i) {
            for (size_t j = i + 1; j < sketches.size(); ++j) {
                double d = jaccard_distance(sketches[i], sketches[j]);
                std::cout << paths[i] << '\t' << paths[j] << '\t' << d << "\n";
            }
        }
    }
    
    return 0;
}
