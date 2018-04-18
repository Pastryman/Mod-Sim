#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#include <iostream>
#include <fstream>
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
const int mc_steps = 200000;
int output_steps = 1;
int min_output_steps = 10000;
const double packing_fraction = 0.2; //0.7
const double diameter = 1.0;
const double dt = 0.00001;
const int maxMeasure = 500;
const double sim_time = 0.2;

/* Volume change -deltaV, delta V */
double delta  = 0.05;
double deltaV = 0.05;

/* Reduced pressure \beta P */
const double temp = 1;
const double betaP = 5;
const char* init_filename = "fcc.dat";
const double rcut = 2.5;
double ecut;
bool liquid = 1;



/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
float r[N][NDIM];
float v[N][NDIM];
float F[N][NDIM];
float F_old[N][NDIM];
double box[NDIM];
double PotE;

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
    std::cout<<"\nold system volume = " <<volume;

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;

}

void update_forces(){
    for(int i = 0; i < n_particles; i++){
        for (int j = 0; j < i; j++){

            double dist2=0;
            double dist[NDIM];
            for (int d = 0; d<NDIM; d++) {
                // Distance between particles
                dist[d] = r[i][d] - r[j][d];
                dist[d] = dist[d] - int(2.0 * dist[d] / box[d]) * box[d];
                dist2 += dist[d]*dist[d];
            }

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
    for(int n = 0; n < n_particles; ++n){
        for(int d = 0; d < NDIM; ++d){
            double r_new = r[n][d] + v[n][d]*dt + F[n][d]*dt*dt/(2.);
            if (r_new<0) {r_new = box[d]-fmod(-r_new,box[d]);}
            if (r_new>box[d]) {r_new = fmod(r_new,box[d]);}
            r[n][d] = float(r_new);
            F_old[n][d] = F[n][d];
        }
    }
    update_forces();
    for(int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d) {
            v[n][d] = v[n][d] + (F[n][d] + F_old[n][d]) * dt / 2.;
        }
    }

}

void initialize_config() {
    float sumv = 0;
    float sumv2 = 0;

    for (int d = 0; d < NDIM; d++){
        sumv = 0;
        for (int n = 0; n < n_particles; n++)
        {
            float v_temp = dsfmt_genrand()*2.-1.;
            v[n][d] = v_temp;
            sumv += v_temp;
            sumv2 += v_temp*v_temp;
        }
        sumv2/=n_particles;
        float fs = sqrt(3.*temp/sumv2);
        //Total momentum (in each dimension) must be zero
        for (int n = 0; n < n_particles; n++){
            v[n][d]=(v[n][d]-sumv/double(n_particles))*fs;
        }

    }

}

int main(int argc, char* argv[]){
    std::cout << "\nBetaP = " << betaP << "\n";

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

    ecut = 4.*(pow(rcut,-12.)-pow(rcut,-6.));

    set_packing_fraction();

    dsfmt_seed(time(NULL));

    std::cout << "\nStart initializing ";
    initialize_config();
    std::cout << "\nInitializing configuration done\nStarting to update forces (first time)";
    update_forces();
    std::cout << "\nUpdating forces done" << std::flush;

    double time=0;
    int step = 0;
    while(time<sim_time)
    {
        update_kinematics();

        if ((step%output_steps) == 0 && step > min_output_steps )
        {
            std::cout << "\ntime = " << time;
            std::cout << "\nPotential Energy = " << PotE;
            write_data(time,r,'r');
        }

        time+=dt;
        step++;
    }


    /*char buffer[128];
    if (liquid) {sprintf(buffer, "Exercise2_P%.0f_liquid.dat",betaP);}
    else {sprintf(buffer, "Exercise2_P_fcc%.0f.dat",betaP);}
    FILE* fp = fopen(buffer, "w");

    printf("\n#Step \t Volume \t Move-acceptance\t Volume-acceptance");
    fprintf(fp,"#Step \t Volume \t Move-acceptance\t Volume-acceptance");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
    double volume_old = liquid?100000:0;
    bool measure = 0;
    int measurements = 0;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            move_accepted += move_particle();
        }
        vol_accepted += change_volume();

        double moveRatio = double(move_accepted) / (double(n_particles) * double(output_steps));
        double volumeRatio = double(vol_accepted) /  double(output_steps);
        if(step % output_steps == 0){
            std::cout << std::flush;
            printf("\n%d \t %lf \t %lf \t %lf",
                   step, box[0] * box[1] * box[2],
                   moveRatio,
                   volumeRatio);
            if(step % 10000 == 0 && !measure)
            {
                double volume_current = box[0]*box[1]*box[2];
                if((volume_current < volume_old && !liquid)||(volume_current > volume_old && liquid))
                {
                    measure = 1;
                    std::cout << "\nMetingen zijn begonnen!\n";
                }
                volume_old = volume_current;
            }
            if(measure)
            {
                fprintf(fp,"\n%d \t %lf \t %lf \t %lf",
                        step, box[0] * box[1] * box[2],
                        moveRatio,
                        volumeRatio);
                measurements++;
                if (measurements>maxMeasure){fclose(fp);return 1;}
            }
            if(moveRatio < 0.35) {delta *= 0.8; std::cout << "\ndelta changed to: " << delta;}
            if(moveRatio > 0.65) {delta *= 1.2; std::cout << "\ndelta changed to: " << delta;}
            if(volumeRatio < 0.35) {deltaV *= 0.8; std::cout << "\ndeltaV changed to: " << deltaV;}
            if(volumeRatio > 0.65) {deltaV *= 1.2; std::cout << "\ndeltaV changed to: " << deltaV;}
            //std::cout << "Vol accepted: " << vol_accepted<< "\n";
            move_accepted = 0;
            vol_accepted = 0;
            //write_data(step);
        }
    }
    fclose(fp);*/

    return 0;
}