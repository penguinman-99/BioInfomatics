#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <tuple>
#include <unordered_map>
using namespace std;
bool isvisit[3];
int arr[8];
int N;
int kmer = 3;
string t = "ATCG";
vector<int> K_mer;//k-mer에서 
vector<vector<int>> result;
enum DNA {
	A, C, G, T
};
 const char *enum_DNA[] = {
	"A","C","G","T"
};
char getEnumDna(int enum_num) {
	 char str = *enum_DNA[enum_num];
	 return str;
}
//HSSP score에 벗어나지 않는 kmer 참 거짓 판단
 bool a(vector<int> vec,int num,int idx) {
	 int score = 0;
	 for (int i = idx; i < ((idx==0)?num:(num+idx)); i++) {
		 if (idx == 0) {
			 if (t[i] != getEnumDna(vec[i]))
				 score -= 1;
			 else
				 score += 1;
		 }
		 else {
			 if (t[i] != getEnumDna(vec[i-idx]))
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
void Get_Permutation(int L) {

	if (L==kmer) {
		for (int i = 0; i < kmer-1; i++) {
			if (a(K_mer, kmer, i)) {
				result.push_back(K_mer);
			}
		}


		return;
	}
		for (int i = 0; i <4; ++i) {
			K_mer.push_back(arr[i]);
			Get_Permutation(L + 1);
			K_mer.pop_back();
		}
	
	return;
}
int main() {
	for (int i = 0; i < 4; i++)
		arr[i] = i;
	Get_Permutation(0);

	for (int i = 0; i < result.size(); i++) {
		int key = 0;
		for (int j = 0; j < kmer; j++) {

			cout << getEnumDna(result[i][j]) << " ";
			key += result[i][j] * pow(4, kmer-1-j);
		}

		cout << key;
		cout << endl;

	}
	cout << K_mer.size() << endl;

}