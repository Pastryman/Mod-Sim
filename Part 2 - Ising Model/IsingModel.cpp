#include <iostream>
#include "mt19937.h"

#define N 100


// Parameters
double J = 1; // J>0: prefers alignment J<0: prefers anti alignment

// Declaration of variables
double H;
int M;

int lattice[N][N];

void initialize_lattice(bool random = true) {

    if (random) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (dsfmt_genrand()>0.5) {lattice[i][j] = 1;}
                else {lattice[i][j]=-1;}
            }
        }
    }

    if (!random) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                lattice[i][j]=1;
            }
        }
    }

}

double H_term(int i, int j, int l,int r,int u,int d){
    // i belongs to u and d
    // j belongs to l and r

    double Htemp = 0;

    Htemp += lattice[i][j] * lattice[u][j];
    Htemp += lattice[j][j] * lattice[d][j];
    Htemp += lattice[i][j] * lattice[i][l];
    Htemp += lattice[i][j] * lattice[i][r];

    Htemp += lattice[i][j] * lattice[u][r];
    Htemp += lattice[i][j] * lattice[u][l];
    Htemp += lattice[i][j] * lattice[d][r];
    Htemp += lattice[i][j] * lattice[d][l];

    return Htemp;
}


void hamiltonian() {

    int iplus;
    int jplus;
    int imin;
    int jmin;
    H = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            iplus = (i + 1)%N;
            jplus = (j + 1)%N;
            imin = i - 1;
            jmin = j - 1;

            if (imin<0) {imin = N - 1;}
            if (jmin<0) {jmin = N - 1;}

            H+=H_term(i,j,jmin,jplus,imin,iplus);
        }
    }
    // Apply coupling strength (between spins)
    H=-H*J/2;

    //std::cout<<"\nH = " << H;
}

void magnetization(){
    M = 0;

    // Loop to get the total magnetization
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            M+=lattice[i][j];
        }
    }

    // Take the average
    M/=N*N;
}

int main() {
    dsfmt_seed(time(NULL));

    initialize_lattice(false);

    lattice[0][0]=-1;

    hamiltonian();

//    int test = 0;
//    for (int n = 0; n < 10000; ++n) {
//        initialize_lattice(false);
//        hamiltonian();
//        test += H;
//    }

    std::cout<<"\nfinal result = " << H;

    return 0;
}