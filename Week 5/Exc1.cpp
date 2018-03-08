#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "sinfft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define points 2048

// Variables
const double diam=1.0;
const double N=500;
const double packing_fraction=0.3;
const double TOL=0.0001;
const double runs=100000;


double dr=0.005; // Note: grid cut-off must be at diam, dr*(int)=1.0;
double gamm[points] = {0.0}; // [Step I]
double gamm_old[points];
double c[points];

void print_file(double dr,double dq, double rho, double pf, int run)
{
    // Initialize file and printing first parameters
    char buffer[128];
    sprintf(buffer, "gList_pf%.2f_rho%.4f.dat", pf, rho);
    FILE *fp = fopen(buffer, "w");
    fprintf(fp, "# Number of particle: \t%d\n", N);
    fprintf(fp, "# dr: \t%f\n", dr);
    fprintf(fp, "# TOL: \t%f\n", TOL);
    fprintf(fp, "# Packing fraction: \t%f\n", packing_fraction);
    fprintf(fp, "# Stopped after number of runs: \t%d\n", run);
    fprintf(fp, "# i \t r \t q \t g[r] \t c[r] \t S[q] \n");

    double S[points]; // For calculation of it see lecture notes equation (7.8)
    double g[points]; // γ(r) = h(r) − c(r) and h(r) = g(r) - 1
    double r;
    double q;


    for(int i=0; i<points; i++)
    {
        r=(i+1)*dr;
        q=(i+1)*dq;

        // NOTE: c is until know in q-space
        // Calculate S[q] with c[q]
        S[i]=1/(1-rho*c[i]);

        // Calculate final c(r) for final gamma(r)
        if(r<diam)
        {
            c[i]=-gamm[i]-1.0;
        }
        else
        {
            c[i]=0.0;
        }
        // NOTE: c is know in r-space

        // Calculate g(r) // γ(r) = h(r) − c(r) and h(r) = g(r) - 1
        g[i]=gamm[i]+c[i]+1;

        // Output step, i, r, q, g(r), c(r), S(q)
        fprintf(fp, "%i \t %lf \t %lf \t %lf \t %lf \t %lf \n",(i+1),r, q, g[i], c[i],S[i]);
    }
    fclose(fp);

}


int main() {

    double V=(M_PI*pow(diam,3.0)*N)/(6*packing_fraction);
    double rho=N/V;
    double r;
    double q;
    double dq=M_PI/(points*dr); // drdq = π/L (see lecture notes)

    for(int run=0; run<runs;run++)
    {
        // [Step II]
        //Calculating c(r)
        for(int i=0;i<points;i++){

            // Saving old value of gamma
            gamm_old[i]=gamm[i];

            r=(i+1)*dr;         // Point on the right border of the grid
            if(r<diam)
            {
                c[i]=-gamm[i]-1;
            }
            else
            {
                c[i]=0.0;
            }

            //std::cout<< "c_r[" <<i << "," << r << "] = \t" << c[i] << "\n";
        }

        // [Step III]
        // Fourier transform c(r) to c(q)

        sinft(c,points);
        for(int i=0;i<points;i++)
        {
            // Important that the pre-factors of the Fourier transform are added
            q=(i+1)*dq;   // Point on the right border of the grid
            c[i]=(4*M_PI*dr)/(q) * c[i];

            //std::cout<< "c_q[" <<i << "," << q << "] = \t" << c[i] << "\n";
        }


        // [Step IV]
        // Calculate gamma^hat = gamma(q)
        for(int i=0;i<points;i++)
        {
            gamm[i]=rho*pow(c[i],2.0)/(1-rho*c[i]);
        }

        // [Step V]
        // Invert gamma(q) to gamma(r) using same Fourier transform
        sinft(gamm,points);
        for(int i=0;i<points;i++)
        {
            r=(i+1)*dr;
            gamm[i]=dq/(2*pow(M_PI,2.0)*r)*gamm[i];

            //std::cout<< "gamma_new[" <<i << "," << r << "] = \t" << gamm[i] << "\n";
        }


        // Next step is to calculate if we have a new solution gamma
        double gamma_difference = 0;
        for(int i=0; i<points; i++)
        {
            gamma_difference+=abs(gamm[i]-gamm_old[i]);
        }
        if(gamma_difference<=TOL)
        {

            print_file(dr,dq,rho,packing_fraction,run);

            std::cout<< "\nEnd solution found: Sum of absolute gamma difference = " << gamma_difference << "\n";
            std::cout<< "In "<< run <<" runs." << "\n";

            return 2;
        }
        std::cout<< "Sum of absolute gamma difference = " << gamma_difference << "\n";

    }

    std::cout<< "\nDid NOT found solution with gamma difference smaller than TOL! Increase runs." << "\n";

    return 1;



}