
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
#include <bitset>
#include <sys/resource.h>
#include <unistd.h>

// スケッチサイズ
constexpr size_t SKETCH_SIZE = 8192;
static_assert(SKETCH_SIZE % 64 == 0);    // 64 の倍数であること

// スケッチをワード単位で保持する型
constexpr size_t NUM_WORDS = SKETCH_SIZE / 64;

// k-mer 長さ（ハッシュ）
constexpr size_t KMER = 64;

// odd sketchの閾値
constexpr double J0 = 0.75;

// odd sketchのハッシュ関数の数（One Permutation Hashingと組み合わせるため）
constexpr uint32_t HASH_NUM = SKETCH_SIZE; // n/{4(1-J0)}=n,since J0=0.75


// kmerをハッシュ化(xxhash)
u_int64_t hash_kmer(const std::string_view &kmer){
    return XXH64(kmer.data(), kmer.length(), 0);
}

// One-Permutation Hashing による MinHash 実装
std::vector<uint64_t> get_minhash_one_permutation(const std::string &seq) {
    // 適切なバケット数を設定（k-merの数より小さく）
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

// FASTA ファイルからシーケンスを読みつつ Odd Sketch ビットセットを構築
std::bitset<SKETCH_SIZE> make_odd_sketch_from_fasta(const std::string &fname) {
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

    // oddsketch本体, デフォルトで 0 に初期化
    std::bitset<SKETCH_SIZE> sketch;  
    
    // One Permutation Hashingで得られた値をスケッチに入れる
    for (size_t i = 0; i < minhash_values.size(); i++) {
        uint64_t hash_val = minhash_values[i];
        size_t pos = hash_val % SKETCH_SIZE;  // ビット位置
        sketch.flip(pos);                     // 反転 (odd‐sketch の核心)
    }
    return sketch;
}

// bitset を raw uint64_t 配列にしてバイナリ書き出す
void write_sketch_binary(const std::bitset<SKETCH_SIZE>& sketch,
                         const std::string& outfname) {
    // SKETCH_SIZE ビットを 64 ビットずつ区切ったワード数
    constexpr size_t NUM_WORDS = SKETCH_SIZE / 64;

    // ビットセットを uint64_t 配列に展開
    std::vector<uint64_t> words(NUM_WORDS, 0);
    for (size_t w = 0; w < NUM_WORDS; ++w) {
        uint64_t word = 0;
        for (size_t b = 0; b < 64; ++b) {
            if (sketch.test(w * 64 + b)) {
                word |= (uint64_t(1) << b);
            }
        }
        words[w] = word;
    }

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
    std::vector<uint64_t> words;  // size() == NUM_WORDS
};

// バイナリファイル (*.sketch) を読み込んで Sketch オブジェクトを返す
Sketch load_sketch(const std::string &fname) {
    Sketch s;
    s.words.resize(NUM_WORDS);
    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open sketch file: " + fname);

    // NUM_WORDS * 8 バイト分を一気読み
    ifs.read(reinterpret_cast<char*>(s.words.data()),
             NUM_WORDS * sizeof(uint64_t));
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
    if (a.words.size() != NUM_WORDS ||
        b.words.size() != NUM_WORDS) {
        throw std::runtime_error("Sketch size mismatch");
    }

    // 1) XOR -> POPCNT
    uint64_t popcnt = 0;
    for (size_t w = 0; w < NUM_WORDS; ++w) {
        popcnt += __builtin_popcountll(a.words[w] ^ b.words[w]);
    }

    // 2*popcnt > sketch_size のときは 0 にクリップ
    if (2 * popcnt > SKETCH_SIZE) {
        return 0.0;
    }

    // 2) 式に代入
    double d_pop    = static_cast<double>(popcnt);
    double d_size   = static_cast<double>(SKETCH_SIZE);
    double d_hashes = static_cast<double>(HASH_NUM);

    double ratio = (2.0 * d_pop) / d_size;            // = 2*popcnt/sketch_size
    double term  = std::log(1.0 - ratio);             // ≤ 0
    double jacc  = 1.0 + (d_size / (4.0 * d_hashes)) * term;
    // double jacc  = 1.0 + (d_size / (4.0 * d_hashes)) * term;

    return jacc;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: oddsketch {sketch|dist}\n";
        return 1;
    }
    std::string mode = argv[1];

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
    else {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 1;
    }
    
    return 0;
}


