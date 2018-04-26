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
#define measurements 400

//GIT TEST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// SET TO TRUE IF YOU WANT TO PRINT DETAILS (for debugging)
const bool debug = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  rhosigma^3 = 1
//  (pi/6)*1.2= pf = 0.63
// google lennard jones phase diagram

/* Initialization variables */
int output_console_steps = 10;
int output_steps = 1;
int initialize_steps = 30000;
const double packing_fraction = 0.6; //0.7
const double diameter = 1.0;
const double dt = 0.001;


/* Temperature variables */
const double temp = 2.;                     // kT of the system
double kT;                                  // Global, used to calculate kinetic temperature of system
const double freq = 100.;                   // freq*∆t=0.1, the frequency collisions,determines the strength to the heat bath
const bool NVT = true;                      // NVT=True, NVE=False

/* Reduced pressure \beta P */
const double betaP = 5;
const char* init_filename = "fcc.dat";
const double rcut = 2.5;
double ecut;


/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
float r[N][NDIM];
float v[N][NDIM];
float v_t[N][NDIM][measurements];
float F[N][NDIM];
float F_old[N][NDIM];
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


    // Resetting Forces and Energy
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
    // Velocity Verlet intergration

    // Update positions r(t)->r(t+dt)
    for(int n = 0; n < n_particles; ++n){
        for(int d = 0; d < NDIM; ++d){

            // Calculating new position
            double r_new = r[n][d] + v[n][d]*dt + F[n][d]*dt*dt/(2.);

            // Nearest image convention (modulus)
            if (r_new<0) {r_new = box[d]-fmod(-r_new,box[d]);}
            if (r_new>box[d]) {r_new = fmod(r_new,box[d]);}

            // r(t)->r(t+dt)
            r[n][d] = float(r_new);

            // Save old force for calculating velocity
            F_old[n][d] = F[n][d];
        }
    }

    // Update forces F(t)->F(t+dt)
    update_forces();

    // Update velocities v(t)->v(t+dt)
    // Calculating temperature: Ekin=(1/2)*m*Σv^2=(d/2)*N*kT // d=3 in 3D

    kT=0;               // Temperature from kinetic energy
    KinE=0;             // Initilize kinetic energy

    for(int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d) {

            // Update velocity v(t)->v(t+dt)
            v[n][d] = v[n][d] + (F[n][d] + F_old[n][d]) * dt / 2.;

            // Ekin=(1/2)*m*Σv^2
            KinE+=v[n][d]*v[n][d];
        }
    }

    kT=KinE/(NDIM*n_particles);
    KinE=KinE/2.0;

    // The Andersen Thermostat
    if(NVT){
        double sigma = sqrt(temp);
        for(int n = 0; n < n_particles; ++n) {
            if(freq*dt>dsfmt_genrand()){
                for (int d = 0; d < NDIM; ++d) {
                    double v_temp = sigma*gaussian_rand(); // Gaussian Distribution
                    v[n][d] = float(v_temp);
                }
            }
        }
    }

}

void initialize_config() {
    // Giving initial velocity to the particles
    // Conditions:
    // - Total momentum = 0
    // - Ekin=(1/2)*m*Σv^2=(d/2)*N*kT // d=3 in 3D

    float sumv = 0;
    float sumv2 = 0;

    for (int d = 0; d < NDIM; d++){
        sumv = 0;
        for (int n = 0; n < n_particles; n++)
        {
            // Giving velocity from [-1,1]
            float v_temp = dsfmt_genrand()*2.-1.;
            v[n][d] = v_temp;
            sumv += v_temp;
            sumv2 += v_temp*v_temp;
        }
        sumv2/=n_particles;

        // Total energy corresponds to temperature
        float fs = sqrt(3.*temp/sumv2);

        //Total momentum (in each dimension) must be zero
        for (int n = 0; n < n_particles; n++){
            v[n][d]=(v[n][d]-sumv/double(n_particles))*fs;
        }

    }

}

int main(int argc, char* argv[]){
    dsfmt_seed(time(NULL));

    // Radius particle
    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    // Initializing configuration
    read_data();
    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }
    set_packing_fraction();

    // Calculate cut off distance, potential is zero at this distance
    ecut = 4.*(pow(rcut,-12.)-pow(rcut,-6.));

    // Initializing configuration, for first step
    std::cout << "\nStart initializing ";
    std::cout << "\nNumber of particles = " << n_particles;
    initialize_config();
    std::cout << "\nInitializing configuration done\nStarting to update forces (first time)";
    update_forces();
    std::cout << "\nUpdating forces done";

    // Initializing local variables for measuring
    bool measuring = false;
    int current_measurements = 0;
    double time=0;
    int step = 0;
    double TotE=0;
    double start_time;

    while(current_measurements < measurements)
    {
        // Update the kinematics (the actual work)
        update_kinematics();

        // Calculate the total energy
        TotE=KinE+PotE;

        // After some steps, we will start measuring
        if (step == initialize_steps)
        {
            measuring = true;
            std::cout << "\nMeasurements have started!!!";
        }

        // Outputting velocity-velocity autocorrelation
        if ((step%output_steps) == 0 && measuring) {
            for (int n = 0; n < n_particles; ++n) {
                for (int dim = 0; dim < NDIM; ++dim) {
                    v_t[n][dim][current_measurements] = v[n][dim];
                }
            }
            current_measurements++;
        }

        if ((step%output_console_steps) == 0) {
            printf("\nTime = %.4f \t PotE = %.2f \t KinE = %.2f \t TotE = %.2f \t kT = %.4f", time, PotE,
                   KinE, TotE, kT);
        }

        time+=dt;
        step++;
    }

    // Initializing .dat file
    char buffer[128];
    sprintf(buffer, "v(0)v(t)_pf%.3f_T%.2f.dat", packing_fraction, temp);
    FILE* fp = fopen(buffer, "w");
    fprintf(fp, "#packing_fraction = %.2f",packing_fraction);
    fprintf(fp, "\n#temperature = %.2f",temp);
    fprintf(fp, "\n#dt = %lf",dt);
    fprintf(fp, "\n#measurements = %i",measurements);

    float v0_vt;

    // Export velocity-velocity autocorrelation
    for (int dt_step = 1; dt_step <= measurements - 100; ++dt_step) {
        v0_vt = 0;
        for (int t0 = 0; t0 < measurements - dt; ++t0) {
            for (int dim = 0; dim < NDIM; ++dim) {
                for (int n = 0; n < n_particles; ++n) {
                    v0_vt += v_t[n][dim][t0]*v_t[n][dim][t0+dt_step];
                }
            }
        }
        v0_vt = v0_vt/(n_particles*(measurements-dt_step));
        printf("\n%lf \t %lf",v0_vt, dt_step*dt);
        fprintf(fp,"\n%lf \t %lf",v0_vt, dt_step*dt);
    }

    return 0;
}