#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512

/* Initialization variables */
const int    mc_steps      = 10000;
const int    output_steps  = 100;
const double density       = 0.8;
const double delta         = 0.1;
const double r_cut         = 2.5;
const double beta          = 0.5;
const char*  init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double e_cut;
double r[N][NDIM];
double box[NDIM];

double energy = 0.0;
double virial = 0.0;

//Struct containing the contribution to the energy and to the virial of a given particle
typedef struct{
    double energy;
    double virial;
}particle_info_t;

//Struct containing the measurements of pressure and excess chemical potential for the whole system
//(to be calculated in the routine measure)
typedef struct{
    double average_pressure;
    double mu_excess;
}measurement_t;

//Function for the calculation of the contribution to energy and to virial of a given particle
particle_info_t particle_energy_and_virial(int);

/* Functions */
measurement_t measure(void){
    measurement_t result;
    /*--------- Your code goes here -----------*/

    return result;
}

//Function that calculates the contribution of the particle to the energy and the virial
//When you use it, be carefull not to make double countings
particle_info_t particle_energy_and_virial(int pid){
    particle_info_t info;
    info.energy = 0.0;
    info.virial = 0.0;
    int n, d;
    for(n = 0; n < n_particles; ++n){ 
        if(n == pid) continue; //For each OTHER particle
	//Its relative distance with pid is calculated (using the MINIMUM IMAGE convention)
        double dist2 = 0.0;
        for(d = 0; d < NDIM; ++d){
            double min_d = r[pid][d] - r[n][d];
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            dist2 += min_d * min_d;
        }

	//And the contribution to energy and virial of the pid-n pair of particles 
        if(dist2 <= r_cut * r_cut){
            double temp = 1.0 / (dist2 * dist2 * dist2);
            info.energy += 4.0 * temp * (temp - 1.0) - e_cut;
            info.virial += 24.0 * temp * (2.0 * temp - 1.0);
        }
    }

    return info;
}

void read_data(void){
    FILE* fp;

    if((fp = fopen(init_filename, "r")) == NULL) {
	printf("Error opening the file \"%s\", program will be arrested.\n", init_filename);
	exit(EXIT_FAILURE);
    }
    int n, d;
    double dmin,dmax;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = abs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
    }
    fclose(fp);
}

int move_particle(void){
    int rpid = n_particles * dsfmt_genrand();


    particle_info_t info = particle_energy_and_virial(rpid);

    double old_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    particle_info_t new_info = particle_energy_and_virial(rpid);

    double dE = new_info.energy - info.energy;
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        virial += new_info.virial - info.virial;
        return 1;
    }

    for(d = 0; d < NDIM; ++d) r[rpid][d] = old_pos[d];

    return 0;
}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void set_density(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}




/* ############################################ MAIN ############################################ */

int main(int argc, char* argv[]){
    /* INITIALIZATION */


    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));


    //The starting configuration is read and check
    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    //Density is set
    set_density();


    //The box size is checked to be big enough given the interaction's critical radius (for MINIMUM IMAGE convention)
    int d;
    for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box[d]);

    //The total initial energy and virial are calculated
    int step, n;
    for(n = 0; n < n_particles; ++n){
        particle_info_t info = particle_energy_and_virial(n);
        energy += info.energy;
        virial += info.virial;
    }
    energy *= 0.5;
    virial *= 0.5;

    //Random number generation is initialized
    size_t seed = time(nullptr);
    dsfmt_seed(seed);

    //Volume is calculated
    double volume = 1.0;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    //And the starting state is printed out
    printf("Starting volume: %f\n", volume);
    printf("Starting energy: %f\n", energy);
    printf("Starting virial: %f\n", virial);
    printf("Starting seed: %lu\n", seed);

    FILE* fp = fopen("measurements.dat", "w");





    /* SIMULATION */


    int accepted = 0;
    for(step = 0; step < mc_steps; ++step){ //For each Monte Carlo cycle
        for(n = 0; n < n_particles; ++n){ //n_particles particle displacements are attempted
            accepted += move_particle();
        }

	//And the observable are measured
        measurement_t ms = measure();
	//And printed out
        fprintf(fp, "%d\t%f\t%f\n", step, ms.average_pressure, ms.mu_excess);

	//Every output_steps the state of the simulation is stored
        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %f.\n",
                step, (double)accepted / (n_particles * output_steps)
            );
            accepted = 0;
            write_data(step);
        }
    }

    fclose(fp);

    return 0;
}
