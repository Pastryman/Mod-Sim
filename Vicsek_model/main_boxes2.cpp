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
using namespace std::chrono;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 2
#define N 4000


//GIT TEST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// SET TO TRUE IF YOU WANT TO PRINT DETAILS (for debugging)
const bool debug = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*System variables*/
int n_particles;
//n_particles = [100,500,1000,1500,2000]

/*System constants*/
const double eta = 0.63;
const double L_box = 20;
const double rcut = 1.0; // WCA: pow(2,(1/6))*diameter; // LJ: 2.5

/* Measurement constants */
const int output_steps = 10;
const int nr_measurements = 10000;
const double equi_time=50;

/* System variables */
double dt = 0.01;
double density=4;
double d_density=2;
double max_density=11;

/* Simulation variables */
double r[N][NDIM];
double theta[N][NDIM-1];
double box[NDIM];
const int no_of_cells = int(ceil(L_box));
int grid[no_of_cells][no_of_cells][N];
int cell_count[no_of_cells][no_of_cells];

double order;

void fill_matrix() {
    int cel[2];
    int cnt;
    memset(grid, 0, sizeof(grid));
    memset(cell_count, 0, sizeof(cell_count));

//    cout << n_particles << "\n";
    for (int n = 0; n < n_particles; n++) {
//        cout << "\nParticle no. " << n << "\n";
        for (int d = 0; d < NDIM; ++d) {
            cel[d] = int(r[n][d] / rcut);
        }
//        cout << "Particle coords: " << r[n][0] << "\t" << r[n][1] << "\n";
//        cout << "Added to matrix: " << cel[0] << "\t" << cel[1] << "\n";
//        cout << "Cell count: " << cell_count[cel[0]][cel[1]] << "\n";
//        cout << "Grid: " <<  grid[cel[0]][cel[1]][cell_count[cel[0]][cel[1]]] << "\n";
        cnt = cell_count[cel[0]][cel[1]];
        cell_count[cel[0]][cel[1]]+=1;
        grid[cel[0]][cel[1]][cnt]=n;
    }
}

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
    fill_matrix();
}


void update_r() {
//    int cel[2];

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
//            for (int d = 0; d < NDIM; ++d) {
//                cel[d] = int(r[n][d] / rcut);
//            }
//            matrix[cel[0]][cel[1]].push_back(n);
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

    fill_matrix();
}

double calc_dist2(int n1, int n2){
    double dist2 = 0;
    for (int d = 0; d < NDIM; ++d) {
        double dist_temp = abs(r[n1][d] - r[n2][d]);

        if (dist_temp > 0.5 * box[d]) {
            dist_temp = box[d] - dist_temp;
        }

        dist2 += pow(dist_temp, 2.0);

    }

    return dist2;
}

void update_theta() {
    double theta_new[N][NDIM - 1];

    for (int n = 0; n < n_particles; ++n) {
        auto cel_x = int(r[n][0] / rcut);
        auto cel_y = int(r[n][1] / rcut);
        int neighbours = 0;
        double theta_nb[NDIM - 1][2] = {0.0};

        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                if (cel_x + i < 0 | cel_y + j < 0 | cel_x + i >= no_of_cells | cel_y + j >= no_of_cells) {
                    continue;
                }

                for (int k = 0; k < cell_count[cel_x + i][cel_y + j]; ++k) {
//                    cout << "Particle: " << n << "\tNeighbouring partcile: " << grid[cel_x + i][cel_y + j][k] << "\n";

                    int n2 = grid[cel_x + i][cel_y + j][k];

                    double dist2 = calc_dist2(n, n2);
                    if (dist2 <= pow(rcut, 2.0)) {
                        for (int l = 0; l < NDIM - 1; ++l) {
                            theta_nb[l][0] += sin(
                                    theta[n2][l]); // Temp translation to let angles run between -pi and pi (for averages)
                            theta_nb[l][1] += cos(
                                    theta[n2][l]); // Temp translation to let angles run between -pi and pi (for averages)
                            neighbours++;
                        }
                    }
                }

            }
        }

        for (int d = 0; d < NDIM - 1; ++d) {
            theta_nb[d][0] /= double(neighbours);
            theta_nb[d][1] /= double(neighbours);
            double mean_angle = atan2(theta_nb[d][0], theta_nb[d][1]);

            theta_new[n][d] = fmod(mean_angle + (gaussian_rand() * eta), 2 * M_PI);
            if (theta_new[n][d] < 0) {
                theta_new[n][d] += 2 * M_PI;
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
    double order_x = 0;
    double order_y = 0;
    for (int n = 0; n < n_particles; ++n) {
        order_x += cos(theta[n][0]);
        order_y += sin(theta[n][0]);
    }
    order = sqrt(order_x*order_x+order_y*order_y);
    order /= n_particles;
}

int main(int argc, char* argv[]){
    cout << "rho = " << density << "\n";
    cout << "Cell size = " << no_of_cells << "\n";
    cout << "N = " << n_particles << "\n";
    cout << "L_box = " << L_box << "\n";
    cout << "rcut = " << rcut << "\n";
    cout << "eta = " << eta << "\n";

    dsfmt_seed(time(NULL));
    initialize_config();

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    while(density<=max_density){
        n_particles = int(L_box*L_box*density);
        initialize_config();

        cout << "\nDensity = "<< density;
        cout << "\nN = " << n_particles << "\n";




        char buffer[128];
        sprintf(buffer, "Measurement_eta%.2f_rho%.3f_dens.dat",eta,density);
        FILE* fp = fopen(buffer, "w");
        fprintf(fp, "# N = %lf \n",float(n_particles));
        fprintf(fp, "# eta = %lf \n",float(eta));
        fprintf(fp, "# rho = %lf \n",float(density));
        fprintf(fp, "# time \t order\n");

        cout << "\nConfiguration initialized\n";

        int measure = 0;
        double time = 0;
        int step=0;

        while (measure<nr_measurements) {
            update_theta();
            update_r();
            measure_order();

            if(step % output_steps==0){
                printf("\nTime: %.3f \t Order: %.5f",time,order);
            }
            if(time>equi_time){
                fprintf(fp, "%.5f \t %lf\n",time,order);
                measure++;
            }
            time+=dt;
            step++;

        }

        fclose(fp);
        density+=d_density;
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "\nDuration was: " << duration;

    return 0;
}