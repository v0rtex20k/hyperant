#include <string>
#include <stdio.h>
#include <iostream>
#include "hypercube.h"

using namespace BMC;
using namespace std;


// Utility printing functions that I used for debugging -
// not needed.

/*
void imprint(int** mat, int n_rows, int n_cols){
    for (int i = 0; i < n_rows; ++i){
        cout << " [ ";
        for (int j = 0; j < n_cols; ++j){
            cout << mat[i][j] << ", ";
        }
        cout << "]," << endl;
    }
}

void fmprint(double** mat, int n_rows, int n_cols){
    for (int i = 0; i < n_rows; ++i){
        cout << "[ ";
        for (int j = 0; j < n_cols; ++j){
            cout << mat[i][j] << ",";
        }
        cout << "]," << endl;
    }
}

void fvprint(double* vec, int len){
    cout << "\n";
    for (int i = 0; i < len; ++i){printf("\t%f\n",vec[i]);}
}*/

int main(int argc, char* argv[]) {
    cout<<flush; // flush homebrew output
    if (argc != 5){printf("Invalid arguments received.\n");}
    else {
        HyperCube HC = HyperCube(stoi(argv[1]), stoi(argv[2]), \
                                   stoi(argv[3]), stoi(argv[4]));
    }
}
