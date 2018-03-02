#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

bool debug = 0;

/* Initialization variables */
const int configs = 500;
const int output_steps = 100;
const int initializeSteps = 2000;
const double packing_fraction = 0.05;
const double diameter = 1.0;
double delta = 0.3;
const char* init_filename = "fcc.dat";
const double dr = 0.1;

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
float r[N][NDIM];
double box[NDIM];

// Make the Histogram array
float hist[N];

/* Functions */
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
    while (i < n_particles && i < N) {
        // Read the values on the line, the last value is ignored
        fscanf(dataFile, "%f %f %f %*f", &r[i][0], &r[i][1], &dummy_r);
        // Print the coordinates of the read particle (for debugging)
        if (debug == 1) { std::cout << "\n" << i + 1 << ": x=" << r[i][0] << ", y=" << r[i][1]; }
        if (NDIM == 3) {
            r[i][2] = dummy_r;
            if (debug) { std::cout << ", z=" << r[i][2]; }
        }
        i++;
    }
}

int move_particle(void){
    int particle =(int)(dsfmt_genrand()*n_particles);
    //std::cout << "\n Particle: " << particle + 1;
    //std::cout << "\n From: " << "x="<< r[particle][0] << " ,y="<< r[particle][1] << " ,z="<< r[particle][2];

    // Generate a random new position of the particle, the new position is (temporarily) stored in another variable
    // This is done in each dimension separately
    double positionAttempt[3];
    for(int i=0;i<NDIM;i++)
    {
        // Generate a random displacement
        double disp =  dsfmt_genrand()*(2.*delta)-delta;
        // Calculate the new position
        positionAttempt[i] = r[particle][i] + disp;
        // Apply boundary conditions
        if(positionAttempt[i] > box[i]) { positionAttempt[i]-=box[i]; }
        if(positionAttempt[i] < 0) { positionAttempt[i]+=box[i]; }
    }

    //std::cout << "\n To: " << "x="<< positionAttempt[0] << " ,y="<< positionAttempt[1] << " ,z="<< positionAttempt[2];

    // Now we check if the particle (in it's new position) doesn't overlap with another particle
    // We do this by checking with every other particle in a for-loop
    for(int i = 0; i<n_particles; i++)
    {
        // We don't check if the particle overlaps with itself
        if(i == particle){ continue;}

        // The displacement between the particles (in x and y)
        double dx = abs(r[i][0]-positionAttempt[0]);
        double dy = abs(r[i][1]-positionAttempt[1]);

        // Apply nearest image convention
        if (dx>.5*box[0]) {dx = box[0] - dx;}
        if (dy>.5*box[1]) {dy = box[1] - dy;}

        // We also want to do this for the z-coordinate, if NDIM == 3
        double dz;
        if (NDIM == 3)
        {
            dz = abs(r[i][2] - positionAttempt[2]);
            if (dz > .5 * box[2]) { dz = box[2] - dz; }
        }

        // We quadraticly sum the coordinate displacements to get the actual displacement
        double dr = pow(dx,2.) + pow(dy,2.);
        if (NDIM == 3) {dr += pow(dz,2.);}

        // If the displacement is smaller than one, we reject the move (by returning)
        if (dr < 1)
        {
            //std::cout << "\nRejected! Overlap with: " << i+1 << "\tdr = " << dr;
            return 0;
        }

    }

    //std::cout<< "\nSuccess!\n";

    // If we get out of the for loop it means the particle does not overlap with any other particle!
    // Perform the displacement and return
    r[particle][0]=float(positionAttempt[0]);
    r[particle][1]=float(positionAttempt[1]);
    r[particle][2]=float(positionAttempt[2]);
    return 1;
}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d_volume%i.dat", step,int(box[0]*box[1]*box[2]));
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

void get_distances(void)
{
    double dist;
    int bin;

    for (int i = 0; i < n_particles; i++)
    {
        for (int j = i ; j < n_particles; j++)
        {
            if (j == i){ continue;}

            // The displacement between the particles (in x and y)
            double dx = abs(r[i][0]-r[j][0]);
            double dy = abs(r[i][1]-r[j][1]);
            double dz = abs(r[i][2] - r[j][2]);
            // Apply nearest image convention
            if (dx>.5*box[0]) {dx = box[0] - dx;}
            if (dy>.5*box[1]) {dy = box[1] - dy;}
            if (dz > .5 * box[2]) { dz = box[2] - dz; }

            dist = sqrt(pow(dx,2.) + pow(dy,2.) + pow(dz,2.));

            bin = int(dist/dr);
            hist[bin] += 2;

        }
    }
}

int main(int argc, char* argv[]) {


    assert(packing_fraction > 0.0 && packing_fraction < 1.0);
    assert(diameter > 0.0);
    assert(delta > 0.0);
    radius = 0.5 * diameter;
    if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if (NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else {
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();
    if (n_particles == 0) {
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }
    set_packing_fraction();
    double volume = box[0] * box[1] * box[2];
    std::cout<< "Volume of system: "<< volume << "\n";

    double rmax = sqrt(3.) * pow(volume, 1. / 3.) / 2.;
    std::cout<< "Maximum distance: "<< rmax << "\n";


    int nbins = int(rmax / dr) + 1;


    dsfmt_seed(time(NULL));
    int accepted = 0;
    int step, n;
    int measurement = 1;
    for (step = 0; step <= output_steps*configs + initializeSteps; step++) {
        for (n = 0; n < n_particles; ++n)
        {
            accepted += move_particle();
        }
        double moveRatio;
        if (step % output_steps == 0)
        {
            moveRatio = double(accepted) / (double(n_particles) * double(output_steps));
            printf("Step %d. Move acceptance: %lf.\n", step, moveRatio);
            accepted = 0;
            if (step>initializeSteps)
            {
                std::cout << "Measuring configuration: " << measurement << "\n";
                get_distances();
                measurement++;
            }
        }
    }

    // Initialize file
    char buffer[128];
    sprintf(buffer, "Exercise2_gList_Volume%i_configs%i.dat", int(volume), configs);
    FILE *fp = fopen(buffer, "w");
    fprintf(fp, "##Volume:\t%lf\n", volume);
    fprintf(fp, "#bin\tr\tgr\n");

    // Scale bincounts
    double nID[nbins];
    double R = 0;
    double g;

    for (int i = 0; i < nbins; i++) {
        hist[i] /= float(configs);
        nID[i] = 4. * M_PI * (500. / volume) * (pow(R + dr, 3.) - pow(R, 3.)) / 3;
        g = hist[i] / (nID[i] * 500.);

        fprintf(fp, "%i \t %lf \t %lf\n",i,R + 0.5 * dr, g);

        R += dr;
    }

    fclose(fp);
    return 1;
}