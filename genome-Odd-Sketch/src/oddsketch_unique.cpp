
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
#include <unordered_set>
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
kmerlen input_num
inputlen1
...
inputlen(input_num)
input1-1
input1-2
*/

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

//ここからスケッチ1計算
        std::priority_queue<uint64_t> sortedkmer;
        std::unordered_set<uint64_t> unique_hashes;
        for (size_t i = 0; i <= input1.length() - kmerlen; i++) {
                std::string kmer = input1.substr(i, kmerlen); // 
                //std::cout << kmer << std::endl;

                uint64_t hash_value = XXH64(kmer.c_str(), kmer.length(), 0); // kmerlenの長さの部分文字列を取り出し、ハッシュ化する
                //std::cout <<"1 "<< kmer <<hash_value<< std::endl;

                if (unique_hashes.find(hash_value) != unique_hashes.end()) {
                    continue;
                }

                if (sortedkmer.size() < sketch_size){
                    sortedkmer.push(hash_value); 
                      
                }
                else if (sortedkmer.top() > hash_value){
                    sortedkmer.pop();
                    sortedkmer.push(hash_value);
                    
                }
                unique_hashes.insert(hash_value);    //もう見た要素は考えないので
        }
        // bottom-numHashFunc個のminを取り出して、スケッチに入れる
        for (uint32_t i = 0; i <= sortedkmer.size(); i++){
            uint64_t min_i = sortedkmer.top();
            sortedkmer.pop();
            
            //std::cout << mkmer << min_i << std::endl;
            int idx_odd1 = (min_i & (sketch_size-1));
            
            oddsketch1[idx_odd1] ^= 1;
            //std::cout<<idx_odd1<<oddsketch1[idx_odd1]<< std::endl;
        }

        //!!!!!!
        //std::cout<<"sketch1 finish"<<std::endl;



//ここからスケッチ２計算
        std::priority_queue<uint64_t> sortedkmer2;
        std::unordered_set<uint64_t> unique_hashes2;
        for (size_t i = 0; i <= input2.length() - kmerlen; i++) {
                std::string kmer = input2.substr(i, kmerlen);
                //std::cout << kmer << std::endl;

                uint64_t hash_value = XXH64(kmer.c_str(), kmer.length(), 0);
                //std::cout <<"1 "<< kmer <<hash_value<< std::endl;

                if (unique_hashes2.find(hash_value) != unique_hashes2.end()) {
                    continue;
                }

                if (sortedkmer2.size() < sketch_size){
                    sortedkmer2.push(hash_value);        
                }
                else if (sortedkmer2.top() > hash_value){
                    sortedkmer2.pop();
                    sortedkmer2.push(hash_value);
                }
                unique_hashes2.insert(hash_value);
        }
        // bottom-numHashFunc個のminを取り出して、スケッチに入れる
        for (uint32_t i = 0; i <= sortedkmer2.size(); i++){
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

        double_t Jaccard = 1 + d_sketch_size/(4*d_numHashFunc)*std::log(1 - (2 * d_popcnt / d_sketch_size));
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


    /*    double_t Jaccard = 1 + sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size);
    std::cout<<"popcnt"<<" "<<popcnt(&oddsketch, sketch_size)<<std::endl;
    std::cout<<"log"<<" "<<log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<"popcnt"<<" "<< sketch_size/(4*numHashFunc)*log(1 - (2 * popcnt(&oddsketch, sketch_size))/sketch_size)<<std::endl;
    std::cout<<Jaccard<<std::endl;*/


    for(int i=0; i<Jaccards.size(); i++){
        std::cout<<Jaccards[i]<<std::endl;
    }
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout<<elapsed<<std::endl;
}