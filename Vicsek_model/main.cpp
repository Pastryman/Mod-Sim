#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#include <iostream>
#include <fstream>
#include <random>
#include "gauss.h"
#include <algorithm>
#include <chrono>  // for high_resolution_clock

using namespace std;


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 2
#define N 2000


//GIT TEST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// SET TO TRUE IF YOU WANT TO PRINT DETAILS (for debugging)
const bool debug = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/* Constants */
const double L_box=1;
const double rcut = 0.05; // WCA: pow(2,(1/6))*diameter; // LJ: 2.5
const double dt = 0.01;

/* Measurement variables */
const int output_steps = 10;
const double equi_time=.1;

/* System variables */
const double eta=0.1;
int n_particles = 100;

/* Simulation variables */
double r[N][NDIM];
double theta[N][NDIM-1];
double box[NDIM];

void initialize_config(void) {
    for (int d = 0; d < NDIM; ++d) {
        box[d] = L_box;
    }

    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d) {
            r[n][d] = dsfmt_genrand() * box[d];
        }
        theta[n][0] = (dsfmt_genrand() * 2 * M_PI);
        if (NDIM == 3) {
            theta[n][1] = (dsfmt_genrand() * M_PI);
        }
    }
}


void update_r() {
    if (NDIM==2){
        double r_new[NDIM];
        for (int n = 0; n < n_particles; ++n) {
            r_new[0]=r[n][0]+dt*cos(theta[n][0]);
            r_new[1]=r[n][1]+dt*sin(theta[n][0]);
            for (int d = 0; d < NDIM; ++d) {
                if(r_new[d]<0){r_new[d]=box[d]-fmod(-r_new[d],box[d]);}
                if(r_new[d]>box[d]){r_new[d]=fmod(r_new[d],box[d]);}
                r[n][d]=r_new[d];
            }
        }
    }
    else if (NDIM==3){
        double r_new[NDIM];
        for (int n = 0; n < n_particles; ++n) {
            r_new[0]=r[n][0]+dt*cos(theta[n][0])*sin(theta[n][1]);
            r_new[1]=r[n][1]+dt*sin(theta[n][0])*sin(theta[n][1]);
            r_new[2]=r[n][2]+dt*cos(theta[n][1]);

            for (int d = 0; d < NDIM; ++d) {
                if(r_new[d]<0){r_new[d]=box[d]-fmod(-r_new[d],box[d]);}
                if(r_new[d]>box[d]){r_new[d]=fmod(r_new[d],box[d]);}
                r[n][d]=r_new[d];
            }
        }
    }
}

void update_theta() {
    double theta_new[N][NDIM - 1];

    for (int n = 0; n < n_particles; ++n) {
        double theta_nb[NDIM - 1] = {0.0};
        int neighbours = 0;

        for (int n_2 = 0; n_2 < n_particles; ++n_2) {

            double dist2 = 0;
            for (int d = 0; d < NDIM; ++d) {
                double dist_temp = abs(r[n][d] - r[n_2][d]);

                if (dist_temp > 0.5 * box[d]) {
                    dist_temp = box[d] - dist_temp;
                }

                dist2 += pow(dist_temp, 2.0);

            }

            if (dist2 <= pow(rcut, 2.0)) {
                for (int k = 0; k < NDIM - 1; ++k) {
                    theta_nb[k] += (theta[n_2][k]-M_PI); // Temp translation to let angles run between -pi and pi (for averages)
                    neighbours++;
                }
            }
        }

        for (int d = 0; d < NDIM - 1; ++d) {
            if (neighbours > 1) {
                theta_nb[d] /= double(neighbours);
                theta_new[n][d] = fmod(theta_nb[d] + (gaussian_rand() * eta) + M_PI, 2*M_PI);
                if (theta_new[n][d]<0){
                    theta_new[n][d]+=2*M_PI;
                }
            }
            else {
                theta_new[n][d]=theta[n][d];
            }
        }
    }

    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM - 1; ++d) {
            theta[n][d] = theta_new[n][d];
        }
    }
}

void write_data(double time){
    char buffer[128];
    sprintf(buffer, "configuration_time%f.dat",time);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    fprintf(fp, "%lf %lf\n", 0.0,1.0);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
        fprintf(fp, "%lf\t", 0.0);
        fprintf(fp, "%lf\n", (box[0]/100.));
    }
    fclose(fp);
}

double measure_order(){
    double order = 0;
    for (int n = 0; n < n_particles; ++n) {
        order += theta[n][0];
    }
    order/=(n_particles*2*M_PI);

    return order;
}

int main(int argc, char* argv[]){

    dsfmt_seed(time(NULL));

    double time;

    initialize_config();

    for (int i = 0; i < 1000; ++i) {
        time=i*dt;
        update_theta();
        update_r();
//        measure_order();
        write_data(time);
    }



    return 0;
}