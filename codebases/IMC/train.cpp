#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <cmath>
#include <omp.h>

#include "util.h"
#include "tron.h"

using namespace std;


void matrix_create(double**& A, int R, int C){
	
	A = new double*[R];
	for(int k=0;k<R;k++)
		A[k] = new double[C];
}

void matrix_delete(double** A, int R, int C){
	
	for(int k=0;k<R;k++)
		delete[] A[k];

	delete[] A;
}

void matrix_vec_mul(Feature* x, double* U, int K, double* xTU){
	
	for(int k=0;k<K;k++)
		xTU[k] = 0.0;
	
	for(int j=0;j<x->size();j++){
		int ind = x->at(j).first;
		double val = x->at(j).second;
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

void rank_one_add(double* mat, double scalar, Feature* x, double* h,  int K){
	
	for(int i=0;i<x->size();i++){
		int ind = x->at(i).first;
		double val = x->at(i).second;
		for(int k=0;k<K;k++){
			mat[ind*K + k] += scalar * (val * h[k]);
		}
	}
}

void transpose(vector<FreqList>& A, int n1, int n2, vector<FreqList>& A_tp){
	
	A_tp.clear();
	A_tp.resize(n2);
	
	for(int i=0;i<n1;i++){
		FreqList::iterator it;
		for(it=A[i].begin();it!=A[i].end();it++){
			A_tp[it->first].push_back(make_pair(i,it->second));
		}
	}
}

void transpose(double** A, int R, int C, double** A_tr){
	
	for(int i=0;i<R;i++)
		for(int j=0;j<C;j++)
			A_tr[j][i] = A[i][j];
}

// Z = XY; X:n1*d1  Y:d1*K 
void matrix_mul( double** Z,    vector<Feature*>* X, double* Y, int K){
	
	Feature* fea;
	Feature::iterator it;
	
	#pragma omp parallel for
	for(int i=0;i<X->size();i++){
		fea = X->at(i);
		for(int k=0;k<K;k++)
			Z[i][k] = 0.0;
		for(it=fea->begin();it!=fea->end();it++){
			for(int k=0;k<K;k++)
				Z[i][k] += it->second* Y[it->first*K + k];
		}
	}
}


void matrix_mul( double** Z,  double** X, double** Y, int n1, int n2, int K){
	
	#pragma omp parallel for 
	for(int i=0;i<n1;i++){
		for(int j=0;j<n2;j++){
			Z[i][j] = 0.0;
			for(int k=0;k<K;k++){
				Z[i][j] += X[i][k] * Y[k][j];
			}
		}
	}
}

double trace(double** A, double** B, double** C, int n1, int K, int n2){
	
	double* tmp = new double[K*K];
	for(int i=0;i<K;i++)
		for(int j=0;j<K;j++)
			tmp[i*K+j] = 0.0;

	#pragma omp parallel for 
	for(int i=0;i<n1;i++){
		for(int j=0;j<K;j++)
			for(int k=0;k<K;k++){
				tmp[j*K+k] += A[i][j] * C[i][k];
			}
	}
	
	double sum = 0.0;
	for(int i=0;i<K;i++){
		for(int j=0;j<K;j++)
			sum += tmp[i*K+j]*B[i][j];
	}

	delete[] tmp;

	return sum;
}

/**   Obj= (2a-1)/2 \sum_{(i,j)\in \Omega} (x_i'Wy_j-1)^2  + (1-a)/2 \|XWY'-A\|^2  + \lambda\|W\|_*
  *   s.t. W = UV' (low-rank)
  */

class ALS_fun : public function {
	
	/**   Min_{U} (2a-1)/2 \sum_{(i,j)\in \Omega} (x_i'U h_j-1)^2  + (1-a)/2 \|XUH'-A\|^2
	  *   where h_j' = y_j'V
	  *   	    (H = YV)
	  *   
	  *   Min_{V} (2a-1)/2 \sum_{(j,i)\in \Omega} (y_j'V h_i-1)^2  + (1-a)/2 \|YVH'-A'\|^2
	  *   where h_i' = x_i'U
	  *	    (H = XU)
	  */
          
	public:
        double alpha;
	ALS_fun(vector<Feature*>* _X, int _d1, double** _H, int _n2, int _K, vector<FreqList >* _clicks, double _lambda){
		
		X = _X;
		n1 = X->size(); //number of queries (or items)
		d1 = _d1; //number of query features (or item features)
		
		H = _H;
		n2 = _n2;
		K = _K; //number of latent dimension
		
		A = _clicks;
		
		//compute alpha = 1 - |click|/|item|  (average)
		int num_clicks = 0;
		for(int i=0;i<A->size();i++){
			FreqList::iterator it;
			for(it=A->at(i).begin(); it!=A->at(i).end(); it++){
				num_clicks += it->second;
			}
		}
		alpha =  1.0; //1.0 - (double)num_clicks / n1 / n2 ;
		cerr << "1-alpha=" << 1.0-alpha << endl;

		lambda = _lambda;
	}
	
	
	double fun(double* U){
		
		// (1-alpha) \sum_{A_ij=1} (f(xi,yj)-1)^2
		Feature* x;
		double* h;
		int j, freq;
		double sum_1=0.0;
		double pred;
		double* tmp = new double[K];
		//#pragma omp parallel for
		for(int i=0;i<A->size();i++){
			x = X->at(i);
			matrix_vec_mul(x, U, K, tmp);
			for(int r=0;r<A->at(i).size();r++){
				j = A->at(i)[r].first;
				freq = A->at(i)[r].second;
				h = H[j];
				
				pred = prod(tmp,h, K);
				
				sum_1 += freq*(pred-1.0)*(pred-1.0);
			}
		}
		sum_1 = (1.0 - alpha)*sum_1;

		// \frac{1}{2}(1-a) \|XWH'-A\|_F^2
		double sum_2 = 0.0;
		
		// term1 = tr(XWH'HW'X')
		double** HT, **HTH;
		matrix_create(HT, K, n2);
		transpose(H, n2, K, HT);
		
		matrix_create(HTH, K, K);
		matrix_mul(HTH, HT, H, K, K, n2);
		matrix_delete(HT, K, n2);
		
		double** XU;
		matrix_create(XU, n1, K);
		matrix_mul( XU, X, U, K);
		
		sum_2 += trace(XU, HTH, XU, n1, K, n1);
		matrix_delete(HTH, K, K);
		matrix_delete(XU, n1, K);
		
		
		// term2 = -2*XWH'A + |A|_F^2
		//#pragma omp parallel for
        /*
		for(int i=0;i<n1;i++){
			x = X->at(i);
			matrix_vec_mul(x,U,K, tmp);
			for(int r=0;r<A->at(i).size();r++){
				j = A->at(i)[r].first;
				h = H[j];
				
				pred = prod(tmp,h, K);
				sum_2 += -2*pred*1.0 + 1.0*1.0;
			}
		}
		delete[] tmp;
        */
		
		sum_2 = (1.0-alpha)*(sum_2)/2.0;
		
		//regualrization
		int d = get_nr_variable();
		double sum_3 = 0.0;
		for(int i=0;i<d;i++)
			sum_3 += U[i]*U[i];
		sum_3 *= lambda/2.0;
		
		return sum_1 + sum_2 + sum_3;
	}

	void grad(double* U, double* g){
		
		int d = get_nr_variable();
		
		double* g1 = new double[d];
		double* g2 = new double[d];
		for(int i=0;i<d;i++){
			g1[i] = 0.0;
			g2[i] = 0.0;
		}
		// (2a-1) \sum_{A_ij=1} (f(xi,yj)-1) (xi'yj) --> g[0]
		// (2a-1) \sum_{A_ij=1} (f(xi,yj)-1) 1[j==j'] --> g[j']
		Feature* x;
		double* h;
		int j, freq;
		double pred, sim_score;
		double* tmp = new double[K];
		//#pragma omp parallel for
		for(int i=0;i<A->size();i++){
			x = X->at(i);
			matrix_vec_mul( x, U, K, tmp );
			for(int r=0;r<A->at(i).size();r++){
				j = A->at(i)[r].first;
				freq = A->at(i)[r].second;
				h = H[j];
				
				pred = prod(tmp, h, K); //item bias
				
				//g1 += (pred-1.0)*x_i h_j';
				rank_one_add(g1, freq*(pred-1.0), x, h, K);
			}
		}
		delete[] tmp;
		for(int i=0;i<d;i++)
			g1[i] *= (2.0*alpha-1.0);
		
		
		// (1-a) X'(XWH'-A)H
		// term1 = X'XWH'H
		double** HT, **HTH;
		matrix_create(HT, K, n2);
		transpose(H, n2, K, HT);
		
		matrix_create(HTH, K, K);
		matrix_mul(HTH, HT, H, K, K, n2);
		matrix_delete(HT, K, n2);
		
		double** XU;
		matrix_create(XU, n1, K);
		matrix_mul( XU, X, U, K);
		
		double** XUHTH;
		matrix_create(XUHTH, n1, K);
		matrix_mul( XUHTH, XU, HTH, n1, K, K);
		matrix_delete(HTH, K, K);
		matrix_delete(XU, n1, K);
		//#pragma omp parallel for 
		for(int i=0;i<n1;i++){
			rank_one_add(g2, 1.0, X->at(i), XUHTH[i], K);
		}
		matrix_delete(XUHTH, n1, K);
		//term 2 = -X'AH
		//#pragma omp parallel for 
		for(int i=0;i<n1;i++){
			x = X->at(i);
			for(int r=0;r<A->at(i).size();r++){
				j = A->at(i)[r].first;
				h = H[j];
				
				//g2 += -1.0 * x_i h_j';
				rank_one_add(g2, -1.0, x, h, K);
			}
		}
		
		for(int i=0;i<d;i++)
			g2[i] *= (1.0-alpha);
		
		for(int i=0;i<d;i++)
			g[i] = g1[i] + g2[i] + lambda*U[i];
		
		
		delete[] g1;
		delete[] g2;
	}
	
	void Hv(double* S, double* Hs){
		
		int d = get_nr_variable();
		
		double* Hs1 = new double[d];
		double* Hs2 = new double[d];
		for(int i=0;i<d;i++){
			Hs1[i] = 0.0;
			Hs2[i] = 0.0;
		}
		
		// (2a-1) \sum_{A_ij=1} (f(xi,yj)-1) (xi'yj) --> Hs[0]
		Feature* x;
		double* h;
		int j, freq;
		double pred_change, sim_score;
		double* tmp = new double[K];
		//#pragma omp parallel for
		for(int i=0;i<A->size();i++){
			x = X->at(i);
			matrix_vec_mul(x, S, K, tmp);
			for(int r=0;r<A->at(i).size();r++){
				j = A->at(i)[r].first;
				freq = A->at(i)[r].second;
				h = H[j];
				
				pred_change = prod( tmp, h,  K );
				
				//Hs1 += pred * x_i h_j';
				rank_one_add( Hs1, freq*pred_change, x, h,  K);
			}
		}
		delete[] tmp;
		for(int i=0;i<d;i++)
			Hs1[i] *= (2.0*alpha-1.0);
		
		// (1-a) X'(XSH')H
		double** HT, **HTH;
		matrix_create(HT, K, n2);
		transpose(H, n2, K, HT);
		
		matrix_create(HTH, K, K);
		matrix_mul(HTH, HT, H, K, K, n2);
		matrix_delete(HT, K, n2);
		
		double** XS;
		matrix_create(XS, n1, K);
		matrix_mul( XS, X, S, K);
		
		double** XSHTH;
		matrix_create(XSHTH, n1, K);
		matrix_mul( XSHTH, XS, HTH, n1, K, K);
		matrix_delete(HTH, K, K);
		matrix_delete(XS, n1, K);
		//#pragma omp parallel for 
		for(int i=0;i<n1;i++){
			rank_one_add( Hs2, 1.0, X->at(i), XSHTH[i], K);
		}
		matrix_delete(XSHTH, n1, K);

		for(int i=0;i<d;i++)
			Hs2[i] *= (1.0-alpha);
		
		for(int i=0;i<d;i++)
			Hs[i] = Hs1[i] + Hs2[i] + lambda*S[i];
		
		delete[] Hs1;
		delete[] Hs2;
	}

	int get_nr_variable(void){

		return d1 * K ;
	}
	
	private:
	
	vector<Feature*>* X;
	double** H;
	vector<FreqList>* A;
	double lambda;
	int n1;
	int d1;
	int n2;
	int K;
};

void train(vector<Feature*>& X, int d1, vector<Feature*>& Y, int d2, vector<FreqList>& A, int K, double* U, double* V, double lambda, int max_iter=30){
	
	int n1 = X.size();
	int n2 = Y.size();	
    int dd1 = d1 * K;
    int dd2 = d2 * K;
	
	//transposed click matrix in sparse format
	vector<FreqList> A_tp;
	transpose(A, n1, n2, A_tp);
	
	//Hq = X*U
	double** Hq = new double*[n1]; 
	for(int i=0;i<n1;i++)
		Hq[i] = new double[K];
	matrix_mul(Hq, &X, U, K);
	
	//Hi = Y*V
	double** Hi = new double*[n2];
	for(int i=0;i<n2;i++)
		Hi[i] = new double[K];
	matrix_mul(Hi, &Y, V, K);
	
	ALS_fun* als_fun_U = new ALS_fun( &X, d1, Hi, n2, K, &A, lambda);
	ALS_fun* als_fun_V = new ALS_fun( &Y, d2, Hq, n1, K, &A_tp, lambda);
	
	TRON* solver_U = new TRON(als_fun_U, 1e-1);
	TRON* solver_V = new TRON(als_fun_V, 1e-1);
	
	//Alternating Least-Square Iterations
	int iter = 0;
    cerr << "iter: " << iter << endl;
	while( iter < max_iter ){
		
		// cerr << "solve U..." << endl;
		solver_U->tron_quad(U);
		//update Hq=XU
		matrix_mul(Hq, &X, U, K);
		
		// cerr << "solve V..." << endl;
		solver_V->tron_quad(V);
		//update Hi=YV
		matrix_mul(Hi, &Y, V, K);
		//cerr << "solved U, obj=" << als_fun_U->fun(U) + ||V||^2<< endl;
		// cout << "iter=" << iter << ", obj=" << als_fun_V->fun(V) + ||U||^2<< endl;
        cerr << "UV solved" << endl;
        //-----------------------------------------------------------
        vector<double> sum1 (n1, 0.0);
        vector<double> sum2 (n1, 0.0);
        // #pragma omp parallel for
        for(int i=0;i<n1;i++){
            double tmp_sum1 = 0.0;
            double tmp_sum2 = 0.0;

            double* UTx = Hq[i];
            double sim_score, pred;
            vector<double> item_score_list;
            item_score_list.resize(n2);
            item_score_list.clear();
            for(int j=0;j<n2;j++){
                double* VTy = Hi[j];
                pred = prod(UTx, VTy, K); //#
                item_score_list.push_back( pred );
                cerr << pred << endl;
            }

            set<int> true_index;
            for (int j=0;j<A[i].size();j++) {
                int index = A[i][j].first;
                tmp_sum1 += (item_score_list[index]-1)*(item_score_list[index]-1);
                true_index.insert(index);
            }
            for (int j=0;j<n2;j++) {
                int index = j;
                set<int>::iterator it = true_index.find(index);
                if (it != true_index.end())
                    tmp_sum2 += item_score_list[index]*item_score_list[index];
            }
            sum1[i] = tmp_sum1;
            sum2[i] = tmp_sum2;
        }
        double tsum1 = 0.0, tsum2 = 0.0;
        for (int i=0;i<n1;i++) tsum1 += sum1[i];
        for (int i=0;i<n1;i++) tsum2 += sum2[i];
        tsum1 *= 1-als_fun_V->alpha;
        tsum2 *= als_fun_V->alpha;
		double reg_U = 0.0, reg_V = 0.0;
        for(int i=0;i<dd1;i++) reg_U += U[i]*U[i];
		reg_U *= lambda/2.0;
        for(int i=0;i<dd2;i++) reg_V += V[i]*V[i];
		reg_V *= lambda/2.0;
        cerr << "sum1=" << tsum1 << ", sum2=" << tsum2 << ", "
            << "reg_U=" << reg_U << ", reg_V=" << reg_V << endl;
        //-----------------------------------------------------------
		iter++;
	}
}


int main(int argc, char** argv){
	
	if( argc < 6 ){
		cerr << "train [A] [X] [Y] [K] [lambda] (model)" << endl;
		cerr <<  endl;
		cerr << "min_{U,V} (1-alpha)* sum_{A_ij=1} (x_i'UV'y_j - 1)^2" << endl
		     << "              alpha* sum_{A_ij=0} (x_i'UV'y_j - 0)^2 + lambda/2*(|U|_F^2 + |V|_F^2)" << endl
		     << "A: n1*n2" << endl
		     << "X: n1*d1" << endl
		     << "Y: n2*d2" << endl
		     << "U: d1*K" << endl
		     << "V: d2*K" << endl;
		exit(0);
	}
	
	srand(1000);
	char* A_file = argv[1];
	char* X_file = argv[2];
	char* Y_file = argv[3];
	int K = atoi(argv[4]);
	double lambda = atof(argv[5]);
	
	char* model_file;
	if( argc > 6 )
		model_file = argv[6];
	else
		model_file = "model";
	

	//Read Features of X and Y
	vector<Feature*> X;
	vector<Feature*> Y;
	int n1, n2, d1, d2;
	readFea(X_file, X, d1);
	readFea(Y_file, Y, d2);
	n1 = X.size();
	n2 = Y.size();
	
	vector< FreqList >  A;
	readMat( A_file, A, n1 );
	
	cerr << "n1=" << n1 << ", n2=" << n2 << ", d1=" << d1 << ", d2=" << d2 << ", K=" << K << ", lambda=" << lambda << endl;
	int U_size = d1*K;
	int V_size = d2*K;
	double* U = new double[U_size];
	double* V = new double[V_size];
	for(int i=0;i<U_size;i++)
		U[i] = (double)rand()/RAND_MAX;
	for(int j=0;j<V_size;j++)
		V[j] = (double)rand()/RAND_MAX;
	
	train( X, d1, Y, d2, A, K,    U, V , lambda);
	
    cerr << "out train.." << endl;
	writeModel( model_file, U, V, d1, d2, K );
	
	return 0;
}
