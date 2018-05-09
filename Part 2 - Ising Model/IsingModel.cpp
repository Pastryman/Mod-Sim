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
int output_data_steps = 10;         // Write data to file and console
int output_console_steps = 50;
int measurements = 1000;            // Amount of measurements we want to do
int initialize_steps = 1000;        // Amount of steps until we start measuring

// Physical Parameters
int J = 1;                          // J>0: prefers alignment J<0: prefers anti alignment
double T;
double T_i = 5.0;
double T_f = 5.0;
double dT = 0.08;

// Declaration of variables
int H;                              // Hamiltonian
double m;                           // Average megnetization

int lattice[N][N];

// Function that what the contribution of energy one lattice site is to the Hamiltionian
// i,j is location of lattice site
// l,r,u,d are left, right, up, down
// Lattice site is compared with all its 8 neighbours
int H_term(int i, int j, int l,int r,int u,int d){
    // i belongs to u and d
    // j belongs to l and r

    int Htemp = 0;

    Htemp += lattice[i][j] * lattice[u][j]; // Up
    Htemp += lattice[i][j] * lattice[d][j]; // Down
    Htemp += lattice[i][j] * lattice[i][l]; // Left
    Htemp += lattice[i][j] * lattice[i][r]; // Right

    Htemp += lattice[i][j] * lattice[u][r]; // Up-Right
    Htemp += lattice[i][j] * lattice[u][l]; // Up-Left
    Htemp += lattice[i][j] * lattice[d][r]; // Down-Right
    Htemp += lattice[i][j] * lattice[d][l]; // Down-Left

    return -Htemp*J;
}

// Function that determines what the neighbours are for a lattice site and calculates its contribution to H
// Ising model makes use of the periodic boundary condition (BC)
int hamiltonian_one(int i, int j){

    int iplus;  // Row down
    int jplus;  // Column right
    int imin;   // Row up
    int jmin;   // Column left

    iplus = (i + 1)%N; // Right BC
    jplus = (j + 1)%N; // Bottom BC
    imin = i - 1;
    jmin = j - 1;

    if (imin<0) {imin = N - 1;} // Left BC
    if (jmin<0) {jmin = N - 1;} // Right BC

    return H_term(i,j,jmin,jplus,imin,iplus);
}

// Function that calculates the total Hamiltionian of the system
int hamiltonian() {
    H = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            H+=hamiltonian_one(i,j);
        }
    }

    // Divide by 2 for overcounting
    H=H/2;

}


// Function that initialize the lattice
// - false: all lattice point point upwards (1)
// - true: random distribution between up (1) and down (-1)
void initialize_lattice(bool random = false) {

    // Random
    if (random) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (dsfmt_genrand()>0.5) {lattice[i][j] = 1;}
                else {lattice[i][j]=-1;}
            }
        }
    }

    // All upwards
    if (!random) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                lattice[i][j]=1;
            }
        }
    }

    hamiltonian();
}

// Function that calculates the average magnitization (Magnitization/#Sites)
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

// Function that tries to perform a flip for a random site
// First calculates the difference in energy for a site flip
// Secondly, accepts the flip according to a Boltzmann distribution
// Lastly, performs flip
// Returns accepted (1) or rejected (0)
int attempt_flip(){

    // Select random particle
    int x = int(dsfmt_genrand()*N);
    int y = int(dsfmt_genrand()*N);

    // The change in energy
    int dH = -2*hamiltonian_one(x, y);

    // Perform the move if Boltzmann distribution is accepted
    if(dsfmt_genrand()<exp(-dH/T)){
        lattice[x][y]*=-1;
        H+=dH;
        return 1; // Accepted
    }

    return 0; // Rejected
}

// Function that performs a 'Sweep'
// Sweep is #Site flip attempts
void MC_sweep(){
    for (int n = 0; n < N * N; ++n) {
        attempt_flip();
    }
}

int main() {

    // Initialize random seed
    dsfmt_seed(time(NULL));

    // Initialize lattice
    initialize_lattice();
    magnetization();
    printf("\nInitial: \t E=%d \t m = %lf", H, m);
    printf("\nInitializing done!\n");

    // Measurements
    // Goes from T_initial to T_final with dT
    for (double b = T_i; b <= T_f; b+=dT) {
        T = b;

        // Initialize output file
        char buffer[128];
        sprintf(buffer, "IsingModel1_T_%.8f.dat", T);
        FILE *fp = fopen(buffer, "w");
        fprintf(fp, "#Step \t Magnetization \t Energy");        // Header: name variables

        int meas = 0;
        int step = 0;

        // Measuring
        while (meas<measurements) {
            MC_sweep();

            // Output data
            if (step % output_data_steps == 0) {
                magnetization();
                fprintf(fp, "\n %d \t %lf \t %d", step, m, H);
                meas++;
            }

            // Write data to console
            if (step % output_console_steps ==0){
                printf("\nstep = %d \t m = %lf \t E = %d", step, m, H);
            }
            step++;
        }

        fclose(fp);
    }
    return 0;
}