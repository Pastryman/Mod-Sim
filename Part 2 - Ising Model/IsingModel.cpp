#include <iostream>
#include "mt19937.h"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "mt19937.h"


#define N 100


// Simulation Parameters
int output_data_steps = 100;    // Write data to file and console
int measurements = 1000;        // Amount of measurements we want to do
int initialize_steps = 100000;  // Amount of steps until we start measuring

// Physical Parameters
int J = 1; // J>0: prefers alignment J<0: prefers anti alignment
double Beta = 1000.;



// Declaration of variables
int H; // Hamiltonian
double m; // Average megnetization

int lattice[N][N];

int H_term(int i, int j, int l,int r,int u,int d){
    // i belongs to u and d
    // j belongs to l and r

    int Htemp = 0;

    Htemp += lattice[i][j] * lattice[u][j];
    Htemp += lattice[i][j] * lattice[d][j];
    Htemp += lattice[i][j] * lattice[i][l];
    Htemp += lattice[i][j] * lattice[i][r];

    Htemp += lattice[i][j] * lattice[u][r];
    Htemp += lattice[i][j] * lattice[u][l];
    Htemp += lattice[i][j] * lattice[d][r];
    Htemp += lattice[i][j] * lattice[d][l];

    return -Htemp*J;
}

int hamiltonian_one(int i, int j){
    int iplus;
    int jplus;
    int imin;
    int jmin;
    iplus = (i + 1)%N;
    jplus = (j + 1)%N;
    imin = i - 1;
    jmin = j - 1;

    if (imin<0) {imin = N - 1;}
    if (jmin<0) {jmin = N - 1;}

    return H_term(i,j,jmin,jplus,imin,iplus);
}

int hamiltonian() {
    H = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            H+=hamiltonian_one(i,j);
        }
    }
    // Apply coupling strength (between spins)
    H=H/2;

    //std::cout<<"\nH = " << H;
}

void initialize_lattice(bool random = false) {

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
                lattice[i][j]=-1;
            }
        }
    }
    hamiltonian();
}

void magnetization(){
    m = 0;

    // Loop to get the total magnetization
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m+=lattice[i][j];
        }
    }

    // Take the average
    m=m/double(N*N);
}

int attempt_flip(){
    // Perform flip of random particle
    int x = int(dsfmt_genrand()*N);
    int y = int(dsfmt_genrand()*N);

    // The change in energy
    int dH = -2*hamiltonian_one(x, y);
    // The move if...
    if(dsfmt_genrand()<exp(-Beta*dH)){
        lattice[x][y]*=-1;
        H+=dH;
        //std::cout<<"\nmove is accepted!";
        return 1;
    }
    return 0;
}

void MC_sweep(){
    for (int n = 0; n < N * N; ++n) {
        attempt_flip();
    }
}

int main() {
    dsfmt_seed(time(NULL));


    initialize_lattice();
    printf("\nInitial: \t E=%d \t m = %lf", H, m);

    printf("Initializing done!");



//    magnetization();
//    hamiltonian();
//    std::cout<<"\nHone = "<< -2*hamiltonian_one(0,0);
//    std::cout<<"\nH = "<<H;
//
//    magnetization();
//    hamiltonian();
//    std::cout<<"\nHone = "<< -2*hamiltonian_one(0,0);
//    std::cout<<"\nH = "<<H;

    for (int i = 0; i < 100; ++i) {
        if(attempt_flip()){
            magnetization();
            std::cout<<"\nH = "<<H;
            std::cout<<"\nm = "<<m;
        }
    }

    for (int step = 0; step < initialize_steps + measurements; ++step) {
        MC_sweep();

        if(step%output_data_steps==0){
            magnetization();
            //hamiltonian();
            printf("\nstep = %d \t E=%d \t m = %lf",step, H, m);
        }
    }

    return 0;
}