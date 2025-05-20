
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

/*
/Volumes/SSD/seq/reference_genomes/GCF_000009605.1_ASM960v1_genomic.fna
/Volumes/SSD/seq/reference_genomes/GCF_000008685.2_ASM868v2_genomic.fna
...
*/


int main(){
    std::chrono::system_clock::time_point  start, end;
    start = std::chrono::system_clock::now();

    const uint32_t kmerlen = 32;
    const uint64_t sketch_size = 1024;
    const double J0 = 0.75;
    const uint32_t numHashFunc = sketch_size; // n/{4(1-J0)}=n,since J0=0.75
    
    std::string line;
    std::vector<std::string> input_genomes;
    while (getline(std::cin, line)) {
        // 読み込んだ行を処理
        //std::cout << "読み込んだ行: " << line << std::endl;
        input_genomes.push_back(line);
    }

    size_t num = input_genomes.size();
    const size_t input_num = num;

    std::vector<std::array<int, sketch_size>> oddsketch(input_num);

    for (size_t first=0; first<input_num; first++){
            //std::cout<< input_genomes[first]<<" " << input_genomes[second]<< std::endl;
            std::ifstream input1_fna(input_genomes[first]);
            if (!input1_fna) {
                std::cout << "ファイルを開けませんでした:"<< std::endl;
                return 1;  // エラー終了
            }

            std::string input1;
            while (std::getline(input1_fna, line)) {
                // FASTA形式のヘッダ行（'>'で始まる）
                if (line[0] == '>') {
                    //std::cout << "ヘッダ: " << line << std::endl;  // ヘッダ行を表示
                } else {
                    input1 += line;
                }
                
            }
            //std::cout << "配列: " << input1 << std::endl;  
//ここから計算
            if (input1.length()<kmerlen){
                std::cout << "The length of the input string is smaller than kmerlen." << std::endl;
                return 1;
            }
            std::cout<<input1.length()<<std::endl;
            
            std::priority_queue<uint64_t> sortedkmer;
            for (size_t i = 0; i <= input1.length() - kmerlen; i++) {
                    std::string kmer = input1.substr(i, kmerlen);
                    //std::cout << kmer << std::endl;

                    uint64_t hash_value = XXH64(kmer.c_str(), kmerlen, 0);
                    //std::cout <<"1 "<< kmer <<hash_value<< std::endl;
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
                
                oddsketch[first][idx_odd1] ^= 1;
                //std::cout<<idx_odd1<<oddsketch1[idx_odd1]<< std::endl;
            }

            //!!!!!!
            std::cout<<first+1<<"sketch finish"<<std::endl;
            end = std::chrono::system_clock::now();
            double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
            std::cout<<"time:"<<time<<std::endl;
            
    }

            /*
            for (uint64_t i=0; i<sketch_size; i++){
                std::cout<<oddsketch1[i]<<std::endl;
            }
            */

                // numHashFunc個のminを取り出して、セットに入れる

            // oddsketch同士でXOR
    std::vector<std::vector<double_t>> Jaccard_list(input_num, std::vector<double_t>(input_num));


    for (size_t first=0; first<input_num; first++){
        for (size_t second=0; second<input_num; second++){
            std::array<int, sketch_size> xor_sketch{};
            uint64_t popcnt = 0;
            for (uint64_t i=0; i<sketch_size; i++){
                xor_sketch[i] = oddsketch[first][i] ^ oddsketch[second][i];
            }
            //std::cout<<oddsketch[i]<<std::endl;
            for (uint64_t i=0; i<sketch_size; i++){
                popcnt += xor_sketch[i];
            }
            double_t d_popcnt = (double)popcnt;
            double_t d_sketch_size = (double)sketch_size;
            double_t d_numHashFunc = (double)numHashFunc;

            double_t Jaccard = 1 + d_sketch_size/(4*d_numHashFunc)*std::log(1 - (2 * d_popcnt / d_sketch_size));
            //std::cout<<"popcnt"<<" "<<popcnt<<std::endl;
            //std::cout<<"login"<<" "<<(1 - (2 * d_popcnt / d_sketch_size))<<std::endl;
            //std::cout<<"log"<<" "<<std::log(1 - (2 * d_popcnt / d_sketch_size))<<std::endl;
            if (2 * popcnt >= sketch_size || Jaccard<J0){
                Jaccard = 0;
            } 
            Jaccard_list[first][second] = Jaccard;
            Jaccard_list[second][first] = Jaccard;
            //std::cout<<"genome"<<first<<" genome"<<second<<" Jaccard係数:"<<Jaccard<<" popcount:"<<popcnt<<std::endl;
            end = std::chrono::system_clock::now();
            double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
            //std::cout<<"time:"<<time<<std::endl;
        }

    }
    for (size_t first=0; first<input_num; first++){
        for (size_t second=0; second<input_num; second++){
            std::cout<<Jaccard_list[first][second];
        }
        std::cout<<" "<<std::endl;
    }

    /*    double_t Jaccard = 1 + sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size);
    std::cout<<"popcnt"<<" "<<popcnt(&oddsketch, sketch_size)<<std::endl;
    std::cout<<"log"<<" "<<log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<"popcnt"<<" "<< sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<Jaccard<<std::endl;*/
    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    //std::cout<<time<<std::endl;
}