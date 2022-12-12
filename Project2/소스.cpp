// reference
// https://www.youtube.com/watch?v=wfi_KimrNQM

#include <iostream>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>
#include <unordered_map>

using namespace std;

int DBlength, merSize;

map<int, pair<string, vector<int>>> kmerMap;
//unordered_map<int, vector<int>> kmerTree;

bool isvisit[3];
int arr[8];
int N;
int kmer_len = 3;
vector<int> K_mer;
vector<pair<int, int>> Query_result;//first: 쿼리 시퀀스에서의 위치 index, second: Key값
vector<pair<int, int>> Seed_idx;
string QuerySeq = "ATCG"; //Query Sequence
enum DNA {
    A, C, G, T
};
const char* enum_DNA[] = {
   "A","C","G","T"
};
//쿼리시퀀스 처리를 위한 메소드
char getEnumDna(int);
bool a(vector<int>, int,int);
void Get_Permutation(int);
int main()
{
    for (int i = 0; i < 4; i++)
        arr[i] = i;
    string DataBase = ""; //DataBase Sequence
    string kmer = "";
    int k, i, j, key, p;
    vector<int> v;
    map<int, pair<string, vector<int>>>::iterator it;
    int vecSize;

    k = 3;
    DataBase = "GGACGGATTCCATGGATA";
    DBlength = DataBase.length();
    merSize = DBlength - k + 1;
    //MER kmers[merSize];

    for (i = 0; i < merSize; i++) {
        key = 0;
        kmer = DataBase.substr(i, k);
        p = i + 1;//position of k-mer
        for (j = 0; j < k; j++) {// calculate key value of k-mer
            switch (DataBase[i + j]) {
            case 'A':
                break;
            case 'C':
                key += (int)pow(4, k - j - 1);
                break;
            case 'G':
                key += 2 * (int)pow(4, k - j - 1);
                break;
            case 'T':
                key += 3 * (int)pow(4, k - j - 1);
                break;
            }
        }

        if (kmerMap.find(key) != kmerMap.end()) { //kmer exists in kmerMap
            kmerMap[key].second.push_back(p); //insert another position
        }
        else {
            v.push_back(p);
            //Automatically sorted in ascending order
            kmerMap.insert(make_pair(key, make_pair(kmer, v))); //add new kmer
        }
        v.clear();
        v.shrink_to_fit();
    }

    //print sorted k-mers
    for (it = kmerMap.begin(); it != kmerMap.end(); it++) {
        cout << "kmer : " << it->second.first << " || position : " << it->second.second[0];
        vecSize = it->second.second.size();
        for (i = 1; i < vecSize; i++) {
            cout << ", " << it->second.second[i];
        }
        //포지션 출력
        cout << " || key : " << it->first << '\n';
    }
    //여기서부터 김균호 학우가 추가하였습니다.
    Get_Permutation(0);

    for (int i = 0; i < Query_result.size(); i++) {
        map<int, pair<string, vector<int>>>::iterator it;
        it = kmerMap.find(Query_result[i].second);
        if(it != kmerMap.end()) {
            //char it->second.first[0];
            int size = it->second.second.size();
            cout << "일치하는 키: " << Query_result[i].second << " 위치: ("<<Query_result[i].first<<","<<it->second.second[0];
            if (size == 1) {
                cout << ")" << endl;
            }
            else {
                for (int i = 1; i < size; i++) {
                    cout << it->second.second[i] << ",";
                }
            cout << endl;
            }

        }
    }
    //for (int i = 0; i < Query_result.size(); i++) {
    //    cout << Query_result[i].first<<" "<<Query_result[i].second<<endl;
    //}

    return 0;
}

//enum과 문자열을 이용해 문자를 바로 반환하게끔함.
char getEnumDna(int enum_num) {
    char str = *enum_DNA[enum_num];
    return str;
}
//HSSP score를 벗어나는지 판단
bool a(vector<int> vec, int num, int idx) {
    int score = 0;
    for (int i = idx; i < ((idx == 0) ? num : (num + idx)); i++) {
        if (idx == 0) {
            if (QuerySeq[i] != getEnumDna(vec[i]))
                score -= 1;
            else
                score += 1;
        }
        else {
            if (QuerySeq[i] != getEnumDna(vec[i - idx]))
                score -= 1;
            else
                score += 1;
        }

    }
    if (score >= 1) {
        return true;
    }
    else
        return false;
}
//kmer 길이만큼의 가능한한 문자열을 모두 뱉어냄(a 함수에 의해 걸러진다)
void Get_Permutation(int L) {

    if (L == kmer_len) {
        for (int i = 0; i < kmer_len - 1; i++) {
            //for (int i = 0; i < result.size(); i++) {
            //    int key = 0;
            //    for (int j = 0; j < kmer_len; j++) {

            //        cout << getEnumDna(result[i][j]) << " ";
            //        key += result[i][j] * pow(4, kmer_len - 1 - j);
            //    }
            //    Query_Keys.push_back(key);
            //    cout << endl;

            //}
            if (a(K_mer, kmer_len,i)) {
                int key = 0;
                for (int i = 0; i < kmer_len; i++) {
                    key += K_mer[i] * pow(4, kmer_len - 1 - i);
                }
                //result.push_back(K_mer);
                Query_result.push_back(make_pair(i+1,key));
            }
        }
        

        return;
    }
    for (int i = 0; i < 4; ++i) {
        K_mer.push_back(arr[i]);
        Get_Permutation(L + 1);
        K_mer.pop_back();
    }

    return;
}