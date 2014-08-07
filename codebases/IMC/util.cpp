#include "util.h"
#include <stdlib.h>

vector<string> split(string str, string pattern){

	vector<string> str_split;
	size_t i=0;
	size_t index=0;
	while( index != string::npos ){

		index = str.find(pattern,i);
		str_split.push_back(str.substr(i,index-i));

		i = index+1;
	}

	if( str_split.back()=="" )
		str_split.pop_back();

	return str_split;
}

void readFea(char* X_file, vector<Feature*>& X, int& d){
	
	X.clear();
	d = -1;
	
	ifstream fin(X_file);
	char* line_cstr = new char[LINE_LEN];
	vector<string> tokens;
	vector<string> k_v;
	while( !fin.eof() ){
		fin.getline(line_cstr, LINE_LEN);
		
		if( fin.eof() )
			break;
		
		string line(line_cstr);
		tokens = split(line, " ");
		
		Feature* fea = new Feature();
		for(int i=0;i<tokens.size();i++){
			k_v = split( tokens[i], ":" );
			int ind = atoi( k_v[0].c_str() );
			double val = atof( k_v[1].c_str() );
			
			fea->push_back(make_pair(ind,val));
			if( ind > d )
				d = ind;
		}
		
		X.push_back(fea);
	}
	fin.close();
	d = d+1;
	
	delete[] line_cstr;
}

void writeModel(char* fname, double* U, double* V, int d1, int d2, int K){
	
	ofstream fout(fname);
	fout << d1 << " " << K << endl;
	for(int i=0;i<d1*K;i++){
		fout << U[i] << " ";
	}
	fout << endl;
	
	fout << d2 << " " << K << endl;
	for(int i=0;i<d2*K;i++){
		fout << V[i] << " ";
	}
	
	fout.close();
}

void readModel(char* fname, double*& U, double*& V, int& d1, int& d2, int& K){
	
	ifstream fin(fname);
	fin >> d1 >> K;
	U = new double[d1*K];
	for(int i=0;i<d1*K;i++){
		fin >> U[i];
	}

	fin >> d2 >> K;
	V = new double[d2*K];
	for(int i=0;i<d2*K;i++){
		fin >> V[i];
	}
	fin.close();
}
