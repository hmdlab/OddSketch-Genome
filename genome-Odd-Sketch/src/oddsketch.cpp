
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

// スケッチサイズ
constexpr size_t SKETCH_SIZE = 8192;
static_assert(SKETCH_SIZE % 64 == 0);    // 64 の倍数であること

// スケッチをワード単位で保持する型
constexpr size_t NUM_WORDS = SKETCH_SIZE / 64;

// k-mer 長さ（ハッシュ）
constexpr size_t KMER = 64;

// odd sketchの閾値
constexpr double J0 = 0.75;

// odd sketchのハッシュ関数の数（bottom-k sketchと組み合わせるため実際は１だが、計算上）
constexpr uint32_t HASH_NUM = SKETCH_SIZE; // n/{4(1-J0)}=n,since J0=0.75

/*
int main(){
    std::chrono::system_clock::time_point  start, end;
    start = std::chrono::system_clock::now();

    const uint32_t kmerlen = 64;
    const uint64_t sketch_size = 8192;
    const double J0 = 0.75;
    const uint32_t numHashFunc = sketch_size; // n/{4(1-J0)}=n,since J0=0.75
    

    size_t k,input_num;
    std::cin>>k>>input_num;
    std::vector<int> input_len;
    for(int i=0; i<input_num; i++){
        int t;
        std::cin>>t;
        input_len.push_back(t);

    }
    std::cout<<input_len[0]<<std::endl;

    std::vector<double> Jaccards;


    for (int i=0; i<input_num; i++){
        std::string input1, input2;
        std::cin>>input1>>input2;

        if (input1.length()<kmerlen){
            std::cout << "The length of the input string is smaller than kmerlen." << std::endl;
            return 1;
        }

        //std::cout<<input1.length()<<" "<< input2.length()<<std::endl;
        std::array<int, sketch_size> oddsketch1{};
        std::array<int, sketch_size> oddsketch2{};

//ここから計算
        std::priority_queue<uint64_t> sortedkmer;
        for (size_t i = 0; i <= input1.length() - kmerlen; i++) {
                std::string kmer = input1.substr(i, kmerlen);
                //std::cout << kmer << std::endl;

                uint64_t hash_value = XXH64(kmer.c_str(), kmerlen, 0);
                //重複するkmerは除外しない（kmerの重複度が同じほど良いため）
                if (sortedkmer.size() < sketch_size){
                    sortedkmer.push(hash_value);        
                }
                else if (sortedkmer.top() > hash_value){
                    sortedkmer.pop();
                    sortedkmer.push(hash_value);
                }
                
        }
        // bottom-numHashFunc個のminを取り出して、スケッチに入れる
        for (uint32_t i = 0; i <= numHashFunc; i++){
            uint64_t min_i = sortedkmer.top();
            sortedkmer.pop();
            
            //std::cout << mkmer << min_i << std::endl;
            int idx_odd1 = (min_i & (sketch_size-1));
            
            oddsketch1[idx_odd1] ^= 1;
            //std::cout<<idx_odd1<<oddsketch1[idx_odd1]<< std::endl;
        }

//ここから計算
        std::priority_queue<uint64_t> sortedkmer2;

        for (size_t i = 0; i <= input2.length() - kmerlen; i++) {
                std::string kmer = input2.substr(i, kmerlen);
                //std::cout << kmer << std::endl;

                uint64_t hash_value = XXH64(kmer.c_str(), kmerlen, 0);
                //std::cout <<"1 "<< kmer <<hash_value<< std::endl;
                if (sortedkmer2.size() < sketch_size){
                    sortedkmer2.push(hash_value);        
                }
                else if (sortedkmer2.top() > hash_value){
                    sortedkmer2.pop();
                    sortedkmer2.push(hash_value);
                }
                
        }
        // bottom-numHashFunc個のminを取り出して、スケッチに入れる
        for (uint32_t i = 0; i <= numHashFunc; i++){
            uint64_t min_i = sortedkmer2.top();
            sortedkmer2.pop();
            
            //std::cout << mkmer << min_i << std::endl;
            int idx_odd2 = (min_i & (sketch_size-1));
            
            oddsketch2[idx_odd2] ^= 1;
            //std::cout<<idx_odd1<<oddsketch1[idx_odd1]<< std::endl;
        }
        //std::cout<<"sketch2 finish"<<std::endl;

        // oddsketch同士でXOR
        std::array<int, sketch_size> oddsketch{};
        uint64_t popcnt = 0;
        for (uint64_t i=0; i<sketch_size; i++){
            oddsketch[i] = oddsketch1[i] ^ oddsketch2[i];
            popcnt += oddsketch[i];
            //std::cout<<oddsketch[i]<<std::endl;
        }

        double_t d_popcnt = (double)popcnt;
        double_t d_sketch_size = (double)sketch_size;
        double_t d_numHashFunc = (double)numHashFunc;

        double_t Jaccard = 1 + d_sketch_size/(2*d_numHashFunc)*std::log(1 - (2 * d_popcnt / d_sketch_size));
        //std::cout<<"popcnt"<<" "<<popcnt<<std::endl;
        //std::cout<<"login"<<" "<<(1 - (2 * d_popcnt / d_sketch_size))<<std::endl;
        //std::cout<<"log"<<" "<<std::log(1 - (2 * d_popcnt / d_sketch_size))<<std::endl;
        if (2 * popcnt > sketch_size){
            Jaccard = 0;
        } 

        //std::cout<<Jaccard<<std::endl;
        Jaccards.push_back(Jaccard);


        end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        std::cout<<elapsed<<std::endl;
    
    }


       double_t Jaccard = 1 + sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size);
    std::cout<<"popcnt"<<" "<<popcnt(&oddsketch, sketch_size)<<std::endl;
    std::cout<<"log"<<" "<<log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<"popcnt"<<" "<< sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<Jaccard<<std::endl;


    for(int i=0; i<Jaccards.size(); i++){
        std::cout<<Jaccards[i]<<std::endl;
    }
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout<<elapsed<<std::endl;
}

*/
// kmerをハッシュ化(xxhash)
u_int64_t hash_kmer(const std::string_view &kmer){
    return XXH64(kmer.data(), kmer.length(), 0);
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

    // bottom-k sketch のための優先度付きキュー
    std::priority_queue<uint64_t> sortedkmer;

    // スライディングウィンドウで k-mer を取り出しハッシュ化、優先度付きキューに格納
    for (size_t i = 0; i + KMER <= seq.size(); i++) {
        std::string_view kmer(&seq[i], KMER);
        uint64_t hash_value = hash_kmer(kmer);// ハッシュ値
        // 重複するkmerは除外しない（kmerの重複度が同じほど良いため）
        if (sortedkmer.size() < HASH_NUM) {
            sortedkmer.push(hash_value);
        } else if (sortedkmer.top() > hash_value) {
            sortedkmer.pop();
            sortedkmer.push(hash_value);
        }
    }

    // oddsketch本体, デフォルトで 0 に初期化
    std::bitset<SKETCH_SIZE> sketch;  
    

    // HASH_NUM 個の最小値を取り出して、スケッチに入れる
    for (size_t i = 0; i < HASH_NUM; i++) {
        uint64_t min_i = sortedkmer.top();
        sortedkmer.pop();
        size_t pos = min_i % SKETCH_SIZE;                // ビット位置
        sketch.flip(pos);                            // 反転 (odd‐sketch の核心)
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
    double jacc  = 1.0 + (d_size / (2.0 * d_hashes)) * term;

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


