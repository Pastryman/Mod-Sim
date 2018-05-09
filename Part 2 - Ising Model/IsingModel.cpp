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
int output_data_steps = 10;    // Write data to file and console
int output_console_steps = 50;
int measurements = 1000;        // Amount of measurements we want to do
int initialize_steps = 1000;  // Amount of steps until we start measuring

// Physical Parameters
int J = 1; // J>0: prefers alignment J<0: prefers anti alignment
double T;
double T_i = 5.0;
double T_f = 5.0;
double dT = 0.08;



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
                lattice[i][j]=1;
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
    // Perform flip of random particle4
    int x = int(dsfmt_genrand()*N);
    int y = int(dsfmt_genrand()*N);

    // The change in energy
    int dH = -2*hamiltonian_one(x, y);
    // Accept the move if MC condition is met
    if(dsfmt_genrand()<exp(-dH/T)){
        lattice[x][y]*=-1;
        H+=dH;
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
    magnetization();
    printf("\nInitial: \t E=%d \t m = %lf", H, m);

    printf("\nInitializing done!\n");


    for (double b = T_i; b <= T_f; b+=dT) {
        T = b;
        // Initialize output file
        char buffer[128];
        sprintf(buffer, "IsingModel1_T_%.8f.dat", T);
        FILE *fp = fopen(buffer, "w");
        fprintf(fp, "#Step \t Magnetization \t Energy");

        int meas = 0;
        int step = 0;
        while (meas<measurements) {
            MC_sweep();

            if (step % output_data_steps == 0) {
                magnetization();
                fprintf(fp, "\n %d \t %lf \t %d", step, m, H);
                meas++;
            }

            if (step % output_console_steps ==0){
                printf("\nstep = %d \t m = %lf \t E = %d", step, m, H);
            }
            step++;
        }

        fclose(fp);
    }
    return 0;
}