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
int main() {
    std::string line;
    std::vector<std::string> input_genomes;
    while (getline(std::cin, line)) {
        // 読み込んだ行を処理
        std::cout << "読み込んだ行: " << line << std::endl;
        input_genomes.push_back(line);
    }

    for (size_t first=0; first<input_genomes.size(); first++){
        for (size_t second=first+1; second<input_genomes.size(); second++){
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
                    std::cout << "ヘッダ: " << line << std::endl;  // ヘッダ行を表示
                } else {
                    input1 += line;
                }
                
            }
            //std::cout << "配列: " << input1 << std::endl;  
        }
    }
}
