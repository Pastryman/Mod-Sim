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

#define NDIM 3
#define N 1000


//GIT TEST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// SET TO TRUE IF YOU WANT TO PRINT DETAILS (for debugging)
const bool debug = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/* Initialization variables */
const int output_steps = 10;
const double packing_fraction = 0.74; //0.7
const double diameter = 1.0;
const double dt = 0.0001;
const double Gamma = 1.0;
const double rcut = pow(2,(1/6))*diameter; // WCA: pow(2,(1/6))*diameter; // LJ: 2.5
const double equi_time=.1;

/* Temperature variables */
const double temp = 2.0;              // kT of the system
const double freq = 0.5;            // freq*∆t the frequency of stochastic collisions which determines the coupling strength to the heat bath

/* Global variables */
const char* init_filename = "fcc.dat";
double ecut;
bool measuring = false;

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
float r[N][NDIM];
float r_cum[N][NDIM];
float v[N][NDIM];
float F[N][NDIM];;
double box[NDIM];
double PotE;
double KinE;

/*
 * Method: read_data
 * Description:
 *  Reads data from a file and stores the results in specified global variables
 *  Processing of the file goes like this:
 *  Line 1:      number of particles in the file
 *  Line 2-4:    dimensions of the box
 *  Line 5-end:  x, y, z, R     coordinates and size of particle (R is never used)
 */
void read_data(void){
    // Open the file and print if the file cannot be opened
    FILE* dataFile;
    if((dataFile = fopen(init_filename,"r"))==NULL)
        printf("file could not be opened\n");

    // Read line 1: Number of particles
    fscanf(dataFile, "%d", &n_particles);
    if (debug){std::cout << "nParticles = " << n_particles << "\n";}

    // Read lines 2-4: dimensions of the box
    float boxDim1, boxDim2;
    int i = 0;
    while (i < NDIM)
    {
        fscanf(dataFile, "\n%f %f",&boxDim1, &boxDim2);
        box[i]=abs(boxDim2 - boxDim1);
        if (debug){std::cout << "boxdim" << i << " =" << box[i] << "\n";}
        i++;
    }

    // Line 5-end: reading coordinates of the particles
    // A "dummy_r" is created to store the value of the third component of each line
    float dummy_r;
    i=0;
    // The number of lines to be read has to be smaller that n_particles or N
    while (i < n_particles && i < N)
    {
        // Read the values on the line, the last value is ignored
        fscanf(dataFile, "%f %f %f %*f",&r[i][0], &r[i][1],&dummy_r);
        // Print the coordinates of the read particle (for debugging)
        if (debug == 1){std::cout << "\n" << i+1 << ": x=" << r[i][0] << ", y=" << r[i][1];}
        if (NDIM == 3){
            r[i][2]=dummy_r;
            if (debug){std::cout << ", z=" << r[i][2];}
        }
        i++;
    }


}

void write_data(double time, float list[][NDIM], char name){
    char buffer[128];
    sprintf(buffer, "%c_time%f.dat", name,time);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", list[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];
    std::cout<<"\nold system volume = " <<volume << "\n";

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;

}

void update_forces(){

    // Resetting Forces and Energy for this iteration
    PotE=0;
    for(int i = 0; i < n_particles; i++){
        for (int d = 0; d<NDIM; d++) {
            F[i][d]=0.0;
        }
    }

    for(int i = 0; i < n_particles; i++){
        for (int j = 0; j < i; j++){
            double dist2=0;
            double dist[NDIM];
            for (int d = 0; d<NDIM; d++) {
                // Distance between particles
                dist[d] = r[i][d] - r[j][d];

                // Nearest image convention
                if(dist[d]>0.5*box[d]){
                    dist[d] = -(box[d]-dist[d]);
                }
                if(dist[d]< -(0.5*box[d])){
                    dist[d] = box[d]+dist[d];
                }

                // Squared distance
                dist2 += dist[d]*dist[d];
            }

            // Cut off distance
            if (dist2 < rcut*rcut){

                // Calculate the force on the particles using the LJ potential
                double r2i=1/dist2;
                double r6i=pow(r2i,3.);
                double ff = 48 * r2i * r6i*(r6i-0.5);

                // Apply the Force to the particles (in opposite direction)
                for(int d=0; d<NDIM;d++){
                    F[i][d] += ff*dist[d];
                    F[j][d] -= ff*dist[d];
                }

                // Update the potential energy
                PotE += 4*r6i*(r6i-1)-ecut;
            }
        }
    }
}

void update_kinematics(){

    // Update forces F(t)->F(t+dt)
    update_forces();

    // Update positions r(t)->r(t+dt)
    for(int n = 0; n < n_particles; ++n){
        for(int d = 0; d < NDIM; ++d){

            // Brownian motion numeric integration
            double r_new=r[n][d]+sqrt(2.0*temp*dt/Gamma)*gaussian_rand()+F[n][d]*dt/Gamma;

            // When measurement start it saves the cumulative distance
            if(measuring){
                r_cum[n][d]+=float(r_new)-r[n][d];
            }

            // Nearest image convention
            if(r_new<0){r_new=box[d]-fmod(-r_new,box[d]);}
            if(r_new>box[d]){r_new=fmod(r_new,box[d]);}
            r[n][d]=float(r_new);
        }
    }

}

int main(int argc, char* argv[]){

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    // Cut off energy, set potential energy to zero at this distance
    ecut = 4.*(pow(rcut,-12.)-pow(rcut,-6.));

    set_packing_fraction();

    // Initializing r cumulative
    for(int i = 0; i < n_particles; i++){
        for (int d = 0; d<NDIM; d++) {
            r_cum[i][d]=0.0;
        }
    }

    // Randomizer
    dsfmt_seed(time(NULL));

    std::cout << "\nStart initializing ";
    std::cout << "\nNumber of particles = " << n_particles;
    std::cout << "\nInitializing configuration done\n";

    // Initialize local variables
    double time=0;
    int step = 0;
    double TotE=0;
    bool stopMeasure = false;
    int nr_of_measurements=0;

    // Initialize output file
    char buffer[128];
    sprintf(buffer, "brownianmotion_fcc_pf%.7f_T%.2f.dat", packing_fraction, temp);
    FILE* fp = fopen(buffer, "w");
    fprintf(fp, "#Time \t <r2_cum (3D)> \t <r_cum>");

    while(!stopMeasure)
    {

        update_kinematics();

        // Time to get to get the system to equilibrium
        if (time>=equi_time && !measuring)
        {
            measuring = true;
            std::cout << "\nMeasurements have started!!!";
        }

        if ((step%output_steps) == 0)
        {
            // Mean-square-displacement at time t
            double r_sum=0;

            // Displacement at time t in each direction
            double x_mean=0;
            double y_mean=0;
            double z_mean=0;

            // Measurements
            if (measuring) {

                // Mean-square-displacement at time t
                for (int n = 0; n < n_particles; ++n) {
                    for (int d = 0; d < NDIM; ++d) {
                        r_sum+=pow(r_cum[n][d],2.0)/n_particles;
                    }

                    // Displacement at time t in each direction
                    x_mean=r_cum[n][0];
                    y_mean=r_cum[n][1];
                    z_mean=r_cum[n][2];
                }

                fprintf(fp,"\n %.4f \t %lf",time,r_sum);

                // Number of measurements
                if (nr_of_measurements>500){stopMeasure = true;}
                nr_of_measurements++;
                }

            // Print in console
            printf("\nTime = %.4f \t PotE = %.2f \t r_cum = %.6f \t x_mean = %.6f \t y_mean = %.6f \t z_mean = %.6f",time,PotE,r_sum,x_mean,y_mean,z_mean);
        }

        // Next iteration
        time+=dt;
        step++;
    }

    return 0;
}