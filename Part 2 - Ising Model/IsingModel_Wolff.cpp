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
int output_data_steps = 50;    // Write data to file and console
int output_console_steps = 1;
int measurements = 1000;        // Amount of measurements we want to do
int initialize_steps = 1000;  // Amount of steps until we start measuring

// Physical Parameters
int J = 1; // J>0: prefers alignment J<0: prefers anti alignment
double T;
double T_i = 5.3;
double T_f = 5.3;
double dT = 0.3;



// Declaration of variables
int H; // Hamiltonian
double m; // Average megnetization

int lattice[N][N];
int cluster_lattice[N][N];

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

void calc_cluster(int x, int y, bool debug = false){
    // Get all neighbours of the particle
    int x_nb[3];
    int y_nb[3];
    x_nb[0] = (x + 1)%N; // x-coords neighbour 1
    y_nb[0] = (y + 1)%N; // y-coords neighbour 1
    x_nb[1] = x; // x-coords neighbour 2
    y_nb[1] = y; // y-coords neighbour 2
    x_nb[2] = x - 1; // x-coords neighbour 3
    y_nb[2] = y - 1; // y-coords neighbour 3
    if (x_nb[2]<0) {x_nb[2] = N - 1;} // boundary conds
    if (y_nb[2]<0) {y_nb[2] = N - 1;} // boundary conds

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i==1 & j==1){ continue;} // don't add self
            else if (cluster_lattice[x_nb[i]][y_nb[j]] == 1) { continue;} // don't add sites twice
                // If the neighbour has the same spin, attempt to add it to the cluster
            else if (lattice[x_nb[i]][y_nb[j]] == lattice[x][y]) {
                // Calculate the chance to add it
                double P_add = 1-exp(-2.*double(J)/T);
                if (dsfmt_genrand() < P_add) {
                    //Add the nb to the cluster
                    cluster_lattice[x_nb[i]][y_nb[j]] = 1;
                    //If we have added it, we want to include the neighbours of the neighbours to the cluster
                    calc_cluster(x_nb[i],y_nb[j],debug);
                    if (debug) {
                        std::cout << "\nsite added to cluster: " << x_nb[i] << "\t" << y_nb[j];
                    }
                }
            }
        }
    }
}

void attempt_flip_cluster(bool debug = false){
    for (int k = 0; k < N; ++k) {
        for (int l = 0; l < N; ++l) {
            cluster_lattice[k][l]=0;
        }
    }

    // Take a particle at random
    auto x = int(dsfmt_genrand()*N);
    auto y = int(dsfmt_genrand()*N);

    calc_cluster(x,y,debug);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (cluster_lattice[i][j]){
                lattice[i][j]*=-1;
            }
        }
    }
}

int main() {
    dsfmt_seed(time(NULL));


    initialize_lattice();
    magnetization();
    printf("\nInitial: \t E=%d \t m = %lf", H, m);

    printf("\nInitializing done!\n");

    T = T_i;

    for (double b = T_i; b <= T_f; b+=dT) {
        T = b;
        // Initialize output file
        char buffer[128];
        sprintf(buffer, "IsingModel1_Wolff_T_%.8f.dat", T);
        FILE *fp = fopen(buffer, "w");
        fprintf(fp, "#Step \t Magnetization \t Energy");

        int meas = 0;
        int step = 0;
        while (meas<measurements) {
            attempt_flip_cluster();

            if (step % output_data_steps == 0) {
                hamiltonian();
                magnetization();
                fprintf(fp, "\n %d \t %lf \t %d", step, m, H);
                printf("\nstep = %d \t m = %lf \t E = %d", step, m, H);
                meas++;
            }
            step++;
        }

        fclose(fp);
    }
    return 0;
}