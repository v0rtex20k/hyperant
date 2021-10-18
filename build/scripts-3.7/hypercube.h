#ifndef ANT_H
#define ANT_H

#include <new>
#include <omp.h>
#include <string>
#include <random>
#include <cstdio>
#include <math.h>
#include <vector>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <nmmintrin.h>

#define ARMA_WARN_LEVEL 1
#include <armadillo>

namespace BMC {
	void no_memory () {
	  printf("Sorry, we can't operate at this dimension - ");
	  printf("too many wormholes and tears in space-time-continuum.");
	  std::exit(1);
	}
    class HyperCube {
    public:
    	int v;
        int dim;
        double* eta;
        HyperCube(int src, int dst, int n_threads=50, int n_dimensions=3){
        	dim = n_dimensions;
        	v = 1 << std::max(n_dimensions,1);

        	omp_set_num_threads(std::max(n_threads,1));
	        std::set_new_handler(no_memory);

        	q   = new double* [v-1];
        	l   = new double* [v-1];
        	u   = new double* [v-1];
        	eta = new double[v-1]();

        	if (n_dimensions < 1 || n_threads < 1){
    			printf("Received invalid parameters.\n");return;
        	}
        	if (src == dst){printf("You're already here! (0s)\n");}
			else if (dst >= v){
				printf("Oops, %d doesn't exist in this dimension (%d is max)", dst, v-1);
				printf(" (inf s - try going to a few dimensions up!)\n");
			}
			else{

				build_matrices(dst);
	        	fill_matrix(l);
	        	fill_matrix(u);
	        	decompose_lu(q);
				sanitize(u,0); // upper triangular
	        	sanitize(l,1); // lower triangular
	        	eta = solve_lu();
	        	int t = avg_time(src);
	        	printf("The estimated travel time from %d (blue)", src);
	        	printf(" to %d (red) is currently %d s.\n", dst, t);
	        }
        }
        ~HyperCube(){
        	delete[] q;
        	delete[] l;
        	delete[] u;
        	delete[] eta;
        }
        void build_matrices(int dst){
			#pragma omp parallel for shared(q,l,u) schedule(static)
			for (int i = 0; i < v; ++i){
				q[i] = new double[v-1]();
				l[i] = new double[v-1]();
				u[i] = new double[v-1]();
			}

			#pragma omp parallel for schedule(static)
			for(int i = 0; i < v; ++i){
				for(int s = 0; s < dim; ++s){
					if (i == dst) {continue;}
					q[i][i^(1<<s)] = 1.0/dim;
					// printf("%d, %d\n", i, i^(1<<s));
				}
			}

			#pragma omp parallel for schedule(static)
			for (int idx = 0; idx < v*v; ++idx){
				int r = idx/v;
				int c = idx%v;
				q[r][c] = ((r == c) ? 1.0 : 0.0) - q[r][c];
			}
		}

		void fill_matrix(double** mat, int diag=0){
			std::random_device rd;
			std::default_random_engine gen(rd());
			std::uniform_real_distribution<double> dist(0.1,1.0);

			#pragma omp parallel for schedule(static)
			for (int idx = 0; idx < v*v; ++idx){
				if (diag && idx/v != idx%v) continue;
				mat[idx/v][idx%v] = dist(gen);
			}
		}

		double* fill_vec(double* vec, double val){
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < v; ++i){
				vec[i] = val;
			}
			return vec;
		}

		void decompose_lu(double** mat){
			int r = 0;
			int k = 0;

			omp_lock_t lock;
			omp_init_lock(&lock);

			for (int x = 0; x < v; ++x){
				for (int y = 0; y < v; ++y){
					if (y < x) l[y][x] = 0;
					else {
						l[y][x] = mat[y][x];
						#pragma omp parallel for schedule(static)
						for (int z = 0; z < x; z++){
							l[y][x] = l[y][x] - l[y][z] * u[z][x];
						}
					}
				}
				for (int y = 0; y < v; ++y){
					if (y  < x) u[x][y] = 0;
					if (y == x) u[x][y] = 1;
					else{
						u[x][y] = mat[x][y] / l[x][x];
						#pragma omp parallel for schedule(static)
						for (int z = 0; z < x; z++){
							u[x][y] = u[x][y] - ((l[x][z] * u[z][y])/l[x][x]);
						}
					}

				}
			}
			omp_destroy_lock(&lock);
		}
		void sanitize(double** mat, int wipe_upper){
			#pragma omp parallel for schedule(static)
			for (int idx = 0; idx < v*v; ++idx){
				int r = idx/v;
				int c = idx%v;
				if (wipe_upper == 1 && r < c) mat[r][c] = 0;
				else if (!wipe_upper && r > c) mat[r][c] = 0;
			}
		}

		void solver(double** mat, double* b){
			omp_lock_t lock;
			omp_init_lock(&lock);

			#pragma omp parallel for collapse(3) schedule(static)
			for (int i = 0; i < v; i++){
				for (int k = i; k < v; k++){
					for (int j = 0; j < v; j++){
						mat[k+1][j] -= (mat[k+1][i]/mat[i][i]) * mat[i][j];
					}
				}
			}

			#pragma omp parallel for collapse(2) schedule(static)
			for (int i = v-1; i >= 0; --i){
				for (int j = v-1; j > i; --j){
					b[i] -= mat[i][j] * b[j];
				}
				b[i] /= mat[i][i];
			}
			omp_destroy_lock(&lock);
		}

		arma::mat _dpp2arma(double** mat, int n_rows, int n_cols){
			arma::mat A(n_rows, n_cols, arma::fill::randu);
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < n_rows; ++i){
				for (int j = 0; j < n_cols; ++j){
					A(i,j) = mat[i][j];
				}
			}
			return A;
		}

		double** _arma2dpp(arma::mat A){
			double** a = new double* [A.n_cols];
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < A.n_rows; ++i){
				for (int j = 0; j < A.n_cols; ++j){
					a[i][j] = A(i,j);
				}
			}
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < A.n_cols; ++i){
				a[i] = new double[A.n_rows]();
			}
			return a;
		}

		double* _arma2dp(arma::vec V){
			double* a = new double[V.n_elem]();
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < V.n_elem; ++i){
					a[i] = V(i);
			}
			return a;
		}

		arma::vec _dp2arma(double* V, int n_elem){
			arma::vec w(n_elem);
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < n_elem; ++i){
					w(i) = V[i];
			}
			return w;
		}

		double* solve_lu(){
			arma::mat L = _dpp2arma(l,v,v);
			arma::mat U = _dpp2arma(u,v,v);
			arma::vec y = arma::zeros(v);

			y = arma::solve(trimatl(L), arma::ones(v), \
				arma::solve_opts::equilibrate + arma::solve_opts::allow_ugly + arma::solve_opts::refine);
			return _arma2dp(arma::solve(trimatu(U), y, \
				arma::solve_opts::equilibrate + arma::solve_opts::allow_ugly + arma::solve_opts::refine));

		}
		int avg_time(int src){
			float sum = eta[src];
			#pragma omp parallel for reduction(+:sum) schedule(static)
			for(int s = 0; s < dim; ++s){
				sum += eta[src^(1<<s)];
			}
			return rint((0.95*sum)/(dim+1)); // adjust for float carryover
		}
    private:
		double** q;
    	double** l;
    	double** u;

    };
}

#endif
