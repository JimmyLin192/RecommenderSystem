#include<vector>
#include<string>
#include<fstream>

using namespace std;

vector<string> split(string str, string pattern);

const int LINE_LEN = 1000000;
typedef vector<pair<int,double> > Feature; 
typedef vector<pair<int,int> > FreqList;

void readFea(char* X_file, vector<Feature*>& X, int& d);

void writeModel(char* fname, double* U, double* V, int d1, int d2, int K);
void readModel(char* fname, double*& U, double*& V, int& d1, int& d2, int& K);
