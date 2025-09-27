
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
std::vector<uint64_t> get_minhash_one_permutation(const std::string &seq) {
    // 適切なバケット数を設定（k-merの数より小さく）
    const size_t HASH_NUM = G_SKETCH_SIZE;
    const size_t NUM_BUCKETS = std::min(static_cast<size_t>(HASH_NUM), seq.size() > KMER ? seq.size() - KMER + 1 : 1);
    
    // 各bucketに対して最小ハッシュ値を保持
    std::vector<uint64_t> bucket_min_hash(NUM_BUCKETS, UINT64_MAX);
    
    // スライディングウィンドウで k-mer を取り出しハッシュ化
    for (size_t i = 0; i + KMER <= seq.size(); i++) {
        std::string_view kmer(&seq[i], KMER);
        uint64_t hash_value = hash_kmer(kmer);
        
        // ハッシュ値をbucket番号に変換（上位ビットを使用）
        size_t bucket_id = (hash_value >> (64 - static_cast<int>(log2(NUM_BUCKETS)) - 1)) % NUM_BUCKETS;
        
        // 各bucketで最小値を更新
        if (hash_value < bucket_min_hash[bucket_id]) {
            bucket_min_hash[bucket_id] = hash_value;
        }
    }
    
    // 結果ベクターを準備（空bucketは除外してbottom-k相当に）
    std::vector<uint64_t> non_empty_values;
    for (size_t i = 0; i < NUM_BUCKETS; i++) {
        if (bucket_min_hash[i] != UINT64_MAX) {
            non_empty_values.push_back(bucket_min_hash[i]);
        }
    }
    
    // bottom-k相当に変換: ソートしてHASH_NUM個まで取る
    std::sort(non_empty_values.begin(), non_empty_values.end());
    std::vector<uint64_t> result(HASH_NUM);
    
    size_t copy_size = std::min(static_cast<size_t>(HASH_NUM), non_empty_values.size());
    for (size_t i = 0; i < copy_size; i++) {
        result[i] = non_empty_values[i];
    }
    
    // 不足分は最大のハッシュ値で埋める（OddSketchでは問題にならない）
    for (size_t i = copy_size; i < HASH_NUM; i++) {
        result[i] = copy_size > 0 ? non_empty_values[copy_size - 1] : 0;
    }
    
    return result;
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

    // One Permutation Hashingで得られた値をスケッチに入れる
    for (size_t i = 0; i < minhash_values.size(); i++) {
        uint64_t hash_val = minhash_values[i];
        size_t pos = static_cast<size_t>(hash_val % G_SKETCH_SIZE);  // ビット位置
        flip_bit(words, pos);                                        // 反転 (odd‐sketch の核心)
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
