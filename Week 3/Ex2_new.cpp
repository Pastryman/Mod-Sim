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


//GIT TEST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// SET TO TRUE IF YOU WANT TO PRINT DETAILS (for debugging)
const bool debug = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/* Initialization variables */
const int mc_steps = 200000;
int output_steps = 100;
const double packing_fraction = 0.74; //0.7
const double diameter = 1.0;
double delta  = 0.05;

const int step_equi=100000;

/* Volume change -deltaV, delta V */
double deltaV = 0.05;

/* Reduced pressure \beta P */
const double betaP = 30.0;
const char* init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
float r[N][NDIM];
double box[NDIM];


/* Functions */
int change_volume(void){

    // Initialize new box parameters, and new particle positions
    double box_new[NDIM];
    float r_new[N][NDIM];

    // Calculate new volume
    double V=box[0]*box[1]*box[2];
    double dV= dsfmt_genrand()*(2.*deltaV)-deltaV;
    double V_new=V+dV;

    //std::cout<< "\ndV: "<<dV;
    //std::cout<<"V_new: "<<V_new<<"\n";
    //std::cout<<"V_old: "<<V<<"\n";

    // Ratio of change for every box dimension
    double ratio=pow(V_new/V,(1.0/3.0));

    // Calculate new particle positions and box dimension after volume change
    for(int i=0;i<NDIM;i++)
    {
        for(int n=0; n<n_particles;n++)
        {
            r_new[n][i] = float(r[n][i]*ratio);
        }
        box_new[i]=box[i]*ratio;
    }

    // Compare every new particle in the new box with all the other new particles in the new box
    // j - particle selected which is compared with the others
    // i - to which particle the selected particle is compared

    if(dV<0)
    {
        for(int j=1;j<n_particles;j++)
        {
            for(int i = 0; i<j; i++)
            {
                // We don't check if the particle overlaps with itself
                if(i == j){ continue;}

                // The displacement between the particles (in x and y)
                double dx = abs(r_new[i][0]-r_new[j][0]);
                double dy = abs(r_new[i][1]-r_new[j][1]);

                // Apply nearest image convention
                if (dx>2.*box_new[0]) {dx = box_new[0] - dx;}
                if (dy>2.*box_new[1]) {dy = box_new[1] - dy;}

                // We also want to do this for the z-coordinate, if NDIM == 3
                double dz;
                if (NDIM == 3)
                {
                    dz = abs(r_new[i][2] - r_new[j][2]);
                    if (dz > 2. * box[2]) { dz = box_new[2] - dz; }
                }

                // We quadraticly sum the coordinate displacements to get the actual displacement
                double dr = pow(dx,2.) + pow(dy,2.);
                if (NDIM == 3) {dr += pow(dz,2.);}

                // If the displacement is smaller than one, we reject the move (by returning)
                if (dr < 1)
                {
                    //std::cout << "\nRejected! Overlap with: " << i+1 << "\tdr = " << dr;
                    //std::cout<< "\nVolume change (rejected,overlap):";
                    return 0;
                }

            }
        }
    }


    double acc_value = exp(-betaP * dV)*pow((V_new/V),n_particles);
    //std::cout<< "\nAcc value: " << acc_value;

    if(!(dsfmt_genrand()<acc_value))
    {
        return 0;
    }

    // Change all the particle positions and volume
    for(int i=0;i<NDIM;i++)
    {
        for(int n=0; n<n_particles;n++)
        {
            r[n][i]=r_new[n][i];
        }
        box[i]=box_new[i];
    }

    //std::cout<< "\nVolume change (accepted): " << dV << "\t" << "Acc value: " << acc_value;
    return 1;

}


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

/*
 * Method: move_particle
 * Description:
 *  Attempt to move a random particle in the box
 *  If the particle does not overlap with another particle the move is accepted
 *  If the particle does overlap with another particle the move is declined, the particle stays in it's original
 *      position
 */
int move_particle(void)
{
    int particle =(int)(dsfmt_genrand()*n_particles);
    //std::cout << "\n Particle: " << particle + 1;
    //std::cout << "\n From: " << "x="<< r[particle][0] << " ,y="<< r[particle][1] << " ,z="<< r[particle][2];

    // Generate a random new position of the particle, the new position is (temporarily) stored in another variable
    // This is done in each dimension separately
    float positionAttempt[3];
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
        if (dx>2.*box[0]) {dx = box[0] - dx;}
        if (dy>2.*box[1]) {dy = box[1] - dy;}

        // We also want to do this for the z-coordinate, if NDIM == 3
        double dz;
        if (NDIM == 3)
        {
            dz = abs(r[i][2] - positionAttempt[2]);
            if (dz > 2. * box[2]) { dz = box[2] - dz; }
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
    r[particle][0]=positionAttempt[0];
    r[particle][1]=positionAttempt[1];
    r[particle][2]=positionAttempt[2];
    return 1;
}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
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

    set_packing_fraction();

    dsfmt_seed(time(NULL));

    char buffer[128];
    sprintf(buffer, "Exercise2_P%.0f_pf%.2f_dV%.3f_dstep%.2f.dat",betaP,packing_fraction,deltaV,delta);
    FILE* fp = fopen(buffer, "w");

    printf("\n#Step \t Volume \t Move-acceptance\t Volume-acceptance");
    fprintf(fp,"#Step \t Volume \t Move-acceptance\t Volume-acceptance");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
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
            if(step>=step_equi)
            {
                //output_steps=10;
                fprintf(fp,"\n%d \t %lf \t %lf \t %lf",
                        step, box[0] * box[1] * box[2],
                        moveRatio,
                        volumeRatio);
            }
            if(moveRatio < 0.35) {delta *= 0.2; std::cout << "\ndelta changed to: " << delta;}
            if(moveRatio > 0.65) {delta *= 1.2; std::cout << "\ndelta changed to: " << delta;}
            if(volumeRatio < 0.35) {deltaV *= 0.2; std::cout << "\ndeltaV changed to: " << deltaV;}
            if(volumeRatio > 0.65) {deltaV *= 1.2; std::cout << "\ndeltaV changed to: " << deltaV;}
            //std::cout << "Vol accepted: " << vol_accepted<< "\n";
            move_accepted = 0;
            vol_accepted = 0;
            //write_data(step);

        }
    }
    fclose(fp);

    return 0;
}