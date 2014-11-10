#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <cmath>
#include <algorithm>
//#include <omp.h>

#include "util.h"
#include "tron.h"

using namespace std;

void matrix_vec_mul(Feature* x, double* U, int d1, int K, double* xTU){
	
	for(int k=0;k<K;k++)
		xTU[k] = 0.0;
	
	for(int j=0;j<x->size();j++){
		int ind = x->at(j).first;
		double val = x->at(j).second;
		if( ind >= d1 )
			continue;
		for(int k=0;k<K;k++){
			xTU[k] += val*U[ind*K + k];
		}
	}
}

double prod(double* xTU, double* h, int K){
	
	double sum =0.0;
	for(int i=0;i<K;i++)
		sum += xTU[i]*h[i] ;

	return sum;
}


void readModel(char* fname, double*& U, int& U_size, double*& V, int& V_size){
	
	ifstream fin(fname);
	fin >> U_size;
	U = new double[U_size];
	for(int i=0;i<U_size;i++){
		fin >> U[i];
	}
	
	fin >> V_size;
	V = new double[V_size];
	for(int i=0;i<V_size;i++){
		fin >> V[i];
	}
	fin.close();
}

bool scoreCompare(pair<int,double> a, pair<int,double> b) { return (a.second > b.second); }

int main(int argc, char** argv){
	
	if( argc < 6 ){
		cerr << "./predict [X] [Y] [model_file] [top#] [pred_file]" << endl;
		exit(0);
	}
	
	char* X_file = argv[1];
	char* Y_file = argv[2];
	char* model_file = argv[3];
	int TOP_K = atoi(argv[4]);
	char* pred_file = argv[5];
	
	vector<Feature*> X;
	vector<Feature*> Y;
	int tmp;
	readFea(X_file, X, tmp);
	readFea(Y_file, Y, tmp);
	
	int n1 = X.size();
	int n2 = Y.size();
	
	int d1, d2, K;
	double* U;
	double* V;
	readModel(model_file, U, V, d1, d2, K);
	
	cerr << "n1=" << n1 << ", n2=" << n2 << ", d1=" << d1 << ", d2=" << d2 << ", K=" << K << endl;
	//construct latent features
	Feature* x;
	vector<double*> UTx_list;
	double* UTx;
	for(int i=0;i<n1;i++){
		
		x = X[i];
		UTx = new double[K];
		matrix_vec_mul(x,U,d1,K,UTx); //#
		UTx_list.push_back(UTx);
	}
	Feature* y;
	vector<double*> VTy_list;
	double* VTy;
	for(int i=0;i<n2;i++){
		
		y = Y[i];
		VTy = new double[K];
		matrix_vec_mul(y,V,d2,K,VTy); //#
		VTy_list.push_back(VTy);
	}
	//predictions
	vector<vector<pair<int,double> >* > predict_list;
	predict_list.resize(n1);
	#pragma omp parallel for
	for(int i=0;i<n1;i++){
		
		double* UTx = UTx_list[i];
		double sim_score, pred;
		vector<pair<int,double> > item_score_list;
		item_score_list.resize(n2);
		item_score_list.clear();
		for(int j=0;j<n2;j++){
			
			double* VTy = VTy_list[j];
			pred = prod(UTx, VTy, K); //#
			item_score_list.push_back( make_pair(j,pred) );
		}
		sort( item_score_list.begin(), item_score_list.end(), scoreCompare );
		vector<pair<int,double> >* tmp = new vector<pair<int,double> >();
		tmp->resize(TOP_K);
		for(int j=0;j<TOP_K;j++){
			pair<int,double> p = item_score_list[j];
			(*tmp)[j] = p;
		}
		predict_list[i] = tmp;
	}
	
    // output pred file
	ofstream fout_pred(pred_file);
	for(int i=0;i<n1;i++){
		vector<pair<int,double> >* tmp = predict_list[i];
		for(int j=0;j<TOP_K;j++){
			pair<int,double> p = tmp->at(j);
			fout_pred << p.first << ":" << p.second << " ";
		}
		fout_pred << endl;
	}
	fout_pred.close();

    // output prec_recall file

    /*
    ofstream fout_prec_recall (prec_recall_file);
    fout_prec_recall << "precision recall" << endl;
    // compute precision and recall

    fout_prec_recall.close();
    */
	
	return 0;
}
