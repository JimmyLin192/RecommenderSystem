#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <cmath>
#include <algorithm>

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
		cerr << "./validate [X] [Y] [model_file] [ground_truth_file] [prec_recall_file]" << endl;
		exit(0);
	}
	
	char* X_file = argv[1];
	char* Y_file = argv[2];
	char* model_file = argv[3];
	char* ground_truth_file = argv[4];
	char* prec_recall_file = argv[5];
	
	vector<Feature*> X;
	vector<Feature*> Y;
	int tmp;
	readFea(X_file, X, tmp);
	readFea(Y_file, Y, tmp);

    // read ground truth
    vector<FreqList> G;
    readMat(ground_truth_file, G, tmp);
	
	int n1 = X.size();
	int n2 = Y.size();
    /*
    if (G.size() != n1) {
        cerr << "inconsistent G.size() and X.size()" << endl;
        exit(-1);    
    }
    */
	
	int d1, d2, K;
	double* U;
	double* V;
	readModel(model_file, U, V, d1, d2, K);
	
	cerr << "n1=" << n1 << ", n2=" << n2 << ", d1=" << d1 << ", d2=" << d2 << ", K=" << K << endl;
	//construct latent topic features
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
    cerr << "for prediction list" << endl;
    // n1 = 10;
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
		tmp->resize(n2);
		for(int j=0;j<n2;j++){
			pair<int,double> p = item_score_list[j];
			(*tmp)[j] = p;
		}
		predict_list[i] = tmp;
	}

    // get max energy and min energy
    double max_value = -1e300, min_value = 1e300;
    for (int i = 0; i < n1; i ++) {
        double max_energy = (*predict_list[i])[0].second;
        double min_energy = (*predict_list[i])[n2-1].second;
        if (max_energy > max_value) 
            max_value = max_energy;
        if (min_energy < min_value) 
            min_value = min_energy;
    }
    cerr << "max_value: " << max_value << endl;
    cerr << "min_value: " << min_value << endl;

    // construct hashmap for each user
    int T = 0;
    vector< set<int> > true_indices (n1, set<int>());
    for(int i=0;i<n1;i++){
        int T_user = G[i].size();
        T += T_user;
        for (int j=0;j<T_user;j++) {
            FreqList tmp = G[i];
            true_indices[i].insert(tmp[j].first);
        }
    }

    // aux: check the position of known entries in predict_list
   // ofstream aux_out ("aux");

    // aux_out.close();

    cerr << "True Label in total = " << T << endl;
    ofstream fout_prec_recall (prec_recall_file);
    // fout_prec_recall << "Recall Precision " << endl;
     fout_prec_recall << "Value Recall Precision " << endl;
    fout_prec_recall.flush();
    int TP = 0;
    int P = 0;
    vector<int> current_index (n1, 0);
    for (double value = max_value; value + 1 >= min_value; ) {
        for (int i=0; i<n1; i++) {
            while (current_index[i] < n2) {
                int j = current_index[i];
                int job_index = (*predict_list[i])[j].first;
                double energy = (*predict_list[i])[j].second;
                if (energy >= value) { // consider this entry
                    set<int>::iterator it = true_indices[i].find(job_index);
                    if (it != true_indices[i].end()) 
                        TP += 1;
                    P += 1;
                    ++ current_index[i];
                } else break;
            }
        }
        int FP = P - TP;  
        int TN = T - TP;  
        double precision = 1.0 * TP / P;
        double recall = 1.0 * TP / T;
        // cerr << recall << " " << precision  << endl;
     //   fout_prec_recall <<  recall  << " " << precision << endl;
        fout_prec_recall << value << " " << recall  << " " << precision << endl;
        fout_prec_recall.flush();
        if (value > -10 && value < 10)
            value -= 0.01;
        else
            value -= 1;
    }
    fout_prec_recall.close();
    return 0;
}
