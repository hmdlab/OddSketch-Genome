#include<set>
#include <unordered_set>
#include<iostream>
#include<vector>
#include <map>

using namespace std;



int main() {
	//ローカル変数として、stを生成
	int k,input_num;
    cin>>k>>input_num;
    vector<int> len_genome;
    for(int i=0; i<input_num; i++){
        int l;
        cin>>l;
        len_genome.push_back(l);
    }

    for(int i=0; i<input_num; i++){
        string input1,input2;
        cin>>input1>>input2;
        unordered_set<string>bunbo,bunshi,kmer1;


        for (int j=0; j<= len_genome[i] - k; j++) {
            string kmer = input1.substr(j, k);
            bunbo.insert(kmer);
            kmer1.insert(kmer);

        }

        for (int j=0; j<= len_genome[i] - k; j++) {
            string kmer = input2.substr(j, k);
            //bunbo.insert(kmer);
            if (kmer1.count(kmer)){
                bunshi.insert(kmer);
            }
            bunbo.insert(kmer);
        }

        double bunshi_size,bunbo_size;
        bunshi_size = double(bunshi.size());
        bunbo_size = double(bunbo.size());

        double Jaccard = bunshi_size/bunbo_size;

        //cout<<bunshi_size<<" "<<bunbo_size<<endl;
        cout<<Jaccard<<endl;


    }
}
