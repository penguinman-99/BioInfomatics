#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
using namespace std;
//2018112008 김균호
int LCS[100][100] = { 0, };
int tmp1[100][100] = { -987654321, };
vector<int> info;
string v, w;
string result1 = "";
string result2 = "";
//TTACCTACG
//ACTACTAT
void PrintLCS(int i, int j) {
	if (LCS[i][j] == 0)
	{
		cout << "END" << endl;
		return;
	}

	if ((i > 0 && j > 0) && (LCS[i - 1][j - 1] + 1 == LCS[i][j]) && (v[i - 1] == w[j - 1])) {
		cout << LCS[i - 1][j - 1] << " ";
		PrintLCS(i - 1, j - 1);
		result1 += v[i - 1];
		result2 += w[j - 1];
	}
	if ((i > 0) && LCS[i - 1][j] - 1 == LCS[i][j]) {
		cout << LCS[i - 1][j] << " ";
		PrintLCS(i - 1, j);
		result1 += v[i - 1];
		result2 += "-";
	}
	if ((j > 0) && LCS[i][j - 1] - 1 == LCS[i][j]) {
		cout << LCS[i][j - 1] << " ";
		PrintLCS(i, j - 1);
		result1 += "-";
		result2 += w[j - 1];
	}
	if ((i > 0 && j > 0) && LCS[i - 1][j - 1] - 1 == LCS[i][j] && (v[i - 1] != w[j - 1])) {
		cout << LCS[i - 1][j - 1] << " ";
		PrintLCS(i - 1, j - 1);
		result1 += v[i - 1];
		result2 += w[j - 1];
	}


}
int main() {
	string tmp;
	cout << "Press string v, w : " << endl;
	//getline(cin, v);
	//getline(cin, w);
	v = " ATCG";

	for (int i = 1; i <= v.length() + 1; i++) {
		for (int j = 1; j <= w.length() + 1; j++) {
			//문자가 같은 경우 기존 길이에서 +1
			int p1, p2 = -987654321;
			if (v[i - 1] != w[j - 1]) {
				p1 = LCS[i - 1][j - 1] - 1;
			}
			else if (v[i - 1] == w[j - 1]) {
				p2 = LCS[i - 1][j - 1] + 1;
				//LCS[i][j] = LCS[i - 1][j - 1] + 1;
			}
			LCS[i][j] = max(0, max(p2, max(p1, max(LCS[i - 1][j] - 1, LCS[i][j - 1] - 1))));

		}
	}
	for (int i = 0; i < v.length(); i++) {
		for (int j = 0; j < w.length(); j++) {
			cout << LCS[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	int max_length = -1;
	for (int i = 0; i < v.length(); i++) {
		for (int j = 0; j < w.length(); j++) {
			if (max_length < LCS[i][j])
				max_length = LCS[i][j];
		}
	}
	for (int i = 0; i < v.length(); i++) {
		for (int j = 0; j < w.length(); j++) {
			if (max_length == LCS[i][j])
				PrintLCS(i, j);
		}
	}
	cout << endl;

	cout << result1 << " " << result2 << " ";
	cout << endl;
	cout << result1.length() << " " << result2.length() << endl;

	return 0;


}