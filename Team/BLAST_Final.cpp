#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <iterator>
#include <cmath>
#include <map>
#include <unordered_map>
using namespace std;
//2018112008 ���ȣ

int DBlength, merSize;

map<int, pair<string, vector<int>>> kmerMap;
//unordered_map<int, vector<int>> kmerTree;

bool isvisit[3];
int arr[8];
int N;
int kmer_len = 3;
vector<int> K_mer;
vector<pair<int, int>> Query_result;//first: ���� ������������ ��ġ index, second: Key��
string QuerySeq = "GACAGC"; //Query Sequence

//���� �� ���ϱ�
int LCS[100][100] = { 0, };
int tmp1[100][100] = { -987654321, };
vector<int> info;
vector<pair<int, int>> location;
vector<vector<pair<int, int>>> location_result;
vector<pair<int, int>> seed_index;
string v_i, w_i;

enum DNA {
	A, C, G, T
};
const char* enum_DNA[] = {
   "A","C","G","T"
};
//���������� ó���� ���� �޼ҵ�
char getEnumDna(int);
bool a(vector<int>, int, int);
void Get_Permutation(int);

void Return_Sequence(string&, string&, vector<pair<int, int>>);
void PrintLCS(int, int);



int main() {
	for (int i = 0; i < 4; i++)
		arr[i] = i;
	string DataBase = ""; //DataBase Sequence
	string kmer = "";
	int k, i, j, key, p;
	vector<int> v;
	map<int, pair<string, vector<int>>>::iterator it;
	int vecSize;

	k = 3;
	DataBase = "ACGGATTCCATATT";
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
		if (kmerMap.find(key) != kmerMap.end()) { //kmer existsin kmerMap
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
		//������ ���
		cout << " || key : " << it->first << '\n';
	}
	//���⼭���� ���ȣ �п찡 �߰��Ͽ����ϴ�.
	
	//Query Sequence�� HSSP score�� 1�� ���յ��� ã�� ���� �Լ�
	Get_Permutation(0);

	for (int i = 0; i < Query_result.size(); i++) {
		map<int, pair<string, vector<int>>>::iterator it;
		it = kmerMap.find(Query_result[i].second);
		if (it != kmerMap.end()) {
			//char it->second.first[0];
			int size = it->second.second.size();
			cout << "Database Sequence�� ��ġ�ϴ� Ű: " << Query_result[i].second 
				<< " ��ġ: (" << Query_result[i].first << "," << it->second.second[0];
			//���� seed index ����
			//seed_index.push_back(make_pair(Query_result[i].first, it->second.second[0]));
			if (size == 1) {
				seed_index.push_back(make_pair(Query_result[i].first, it->second.second[0]));
				cout << ")" << endl;
			}
			else {
				cout << ",";
				for (int i = 1; i < size; i++) {
					seed_index.push_back(make_pair(Query_result[i].first, it->second.second[i]));
					cout << it->second.second[i] << ",";
				}
				cout << ")";
				cout << endl;
			}

		}
	}




//���̽� ���͸� �˰����� �̿��� ���̺� �����
	string a = " ";
	string tmp;
	v_i = "GACAGC";
	w_i = "ACGGATTCCATATT";
	v_i = a + v_i; //ATCG���� ����ִ� ���� �ϳ� �߰� ��������.
	w_i = a + w_i;
	//w = w;
	//ù �� ���� 0���� �ʱ�ȭ. 
	for (int i = 0; i < v_i.length(); i++) {
		LCS[i][0] = 0;
		if (v_i[i] == w_i[0]) {
			LCS[i][0] = 1;
		}
	}
	for (int i = 0; i < w_i.length(); i++) {
		LCS[0][i] = 0;
	}
	for (int i = 0; i < v_i.length(); i++) {
		for (int j = 0; j < w_i.length(); j++) {
			//���ڰ� ���� ��� ���� ���̿��� +1
			int p1, p2 = -987654321;
			if (i == 0 || j == 0) {
				continue;
			}
			else if (v_i[i] != w_i[j]) {
				p1 = LCS[i-1][j-1]- 1;
			}
			else if (v_i[i] == w_i[j]) {
				p2 = LCS[i - 1][j - 1] + 1;
			}
			//���̺��� ��ȭ��. 0, ������+���ڽ��� ��ġ�Ѱ��,+1, ������+�� �ڽ��� ����ġ�� ���,-1, indel -1
			LCS[i][j] = max(0, max(p2, max(p1, max(LCS[i - 1][j] - 1, LCS[i][j - 1] - 1))));

		}
	}
	//���̺� �����ֱ�
	for (int i = 0; i < v_i.length(); i++) {
		for (int j = 0; j < w_i.length(); j++) {
			cout << LCS[i][j] << " ";
		}
		cout << endl;
	}

	//���� ��� ���
	sort(seed_index.begin(), seed_index.end());
	for (int k = 0; k < seed_index.size(); k++) {
		PrintLCS(seed_index[k].first,seed_index[k].second);
		for (int i = 0; i < location_result.size(); i++) {
			if (LCS[seed_index[k].first][seed_index[k].second] == 0) {

			}
			else {
				location_result[i].push_back(make_pair(seed_index[k].first, seed_index[k].second));
				rotate(location_result[i].rbegin(), location_result[i].rbegin() + 1, location_result[i].rend());
			}

			//��� ���
			for (int j = 0; j < location_result[i].size(); j++) {
				cout << "(" << location_result[i][j].first << "," << location_result[i][j].second << ")" << " ";
			}
			string tmp1 = "";
			string tmp2 = "";
			//��ο� ���� ���� ��� ���
			Return_Sequence(tmp1, tmp2, location_result[i]);
			cout << tmp1 << " " << tmp2;
			cout << endl;
		}
		location_result.clear();
		location.clear();
	}

	return 0;


}

//enum�� ���ڿ��� �̿��� ���ڸ� �ٷ� ��ȯ�ϰԲ���.
char getEnumDna(int enum_num) {
	char str = *enum_DNA[enum_num];
	return str;
}
//HSSP score�� ������� �Ǵ�
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
//kmer ���̸�ŭ�� �������� ���ڿ��� ��� ��(a �Լ��� ���� �ɷ�����)
void Get_Permutation(int L) {
	if (L == kmer_len) {
		for (int i = 0; i < kmer_len - 1; i++) {
			if (a(K_mer, kmer_len, i)) {
				int key = 0;
				for (int i = 0; i < kmer_len; i++) {
					key += K_mer[i] * pow(4, kmer_len - 1 - i);
				}
				Query_result.push_back(make_pair(i + 1, key));
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

//��ǥ ������ �� �ľ��Ͽ�, ���� DNA ���⼭���� ����մϴ�.
void Return_Sequence(string& a, string& b, vector<pair<int, int>> location) {
	string tmp_i = "";
	tmp_i += v_i[location[0].first];
	string tmp_j = "";
	tmp_j += w_i[location[0].second];
	for (int i = 0; i < location.size() - 1; i++) {
		//��ġ, ����ġ 
		if (location[i + 1].first - location[i].first == 1 && location[i + 1].second - location[i].second == 1) {
			tmp_i += v_i[location[i + 1].first];
			tmp_j += w_i[location[i + 1].second];
		}
		//���η� �� ���
		if (location[i + 1].first - location[i].first == 0 && location[i + 1].second - location[i].second == 1) {

			tmp_i += "-";
			tmp_j += w_i[location[i + 1].second];
		}
		//���η� �����
		if (location[i + 1].first - location[i].first == 1 && location[i + 1].second - location[i].second == 0) {
			tmp_i += v_i[location[i + 1].first];
			tmp_j += "-";
		}
	}
	a = tmp_i;
	b = tmp_j;
}
void PrintLCS(int i, int j) {
	//�� �̻� �����ٷ� �� �� ����, �����ε� ������ �� ���� ���
	if (i == v_i.length() - 1 && (LCS[i][j + 1] + 1 != LCS[i][j])) {
		location_result.push_back(location);
	}
	//�� �̻� ���� �ٷ� �� �� ����, �����ε� ������ �� ���� ���
	else if (j == w_i.length() - 1 && (LCS[i + 1][j] + 1 != LCS[i][j])) {
		location_result.push_back(location);
	}
	//1. ���� 1�̰�, ����, ���� �� ĭ �� �� ���� 2�� �ƴ� ���
	else if ((LCS[i][j] == 1 && (LCS[i + 1][j + 1] - 1 != LCS[i][j])))
	{
		location_result.push_back(location);
	}

	else {
		//��ġ�߰ų�, ����ġ�߰ų�
		if ((i < v_i.length() - 1 && j < w_i.length() - 1) && (i > 0 && j > 0) && 
			(LCS[i][j] - 1 == LCS[i + 1][j + 1] || (LCS[i][j] + 1 == LCS[i + 1][j + 1]))) {
			location.push_back(make_pair(i + 1, j + 1));
			PrintLCS(i + 1, j + 1);
			location.pop_back();
		}
		//���η� ������
		if (LCS[i + 1][j] != 0 && (i < v_i.length() - 1) && LCS[i + 1][j] + 1 == LCS[i][j]) {
			location.push_back(make_pair(i + 1, j));
			PrintLCS(i + 1, j);
			location.pop_back();
		}
		//���η� ������
		if (LCS[i][j + 1] != 0 && (j < w_i.length() - 1) && LCS[i][j + 1] + 1 == LCS[i][j]) {
			location.push_back(make_pair(i, j + 1));
			PrintLCS(i, j + 1);
			location.pop_back();
		}
	}
}