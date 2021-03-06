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

// Constant variables
const double diam=1.0;
const int N=500;
const double TOL=0.0001;
const double runs=1000;
const double alpha=0.5;

// Variables
double packing_fraction=0.0;
double dr=0.005; // 0.005           // Note: grid cut-off must be at diam, dr*(a integer)=1.0;
double gamm[points] = {0.0};        // [Step I]: Initial guess for the gamma[r]
double gamm_old[points];
double c[points];


// Function that export: g(r), c(r), S(q) over all used r and q space
void print_file(double dr,double dq, double rho, double pf, int run)
{
    // Initialize file and printing first parameters
    char buffer[128];
    sprintf(buffer, "gList_pf%.2f.dat", pf);
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
        r=(i)*dr;
        q=(i)*dq;

        // NOTE: c is until now in q-space
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
        // NOTE: c is now in r-space
        // Calculate g(r) // γ(r) = h(r) − c(r) and h(r) = g(r) - 1
        g[i]=gamm[i]+c[i]+1;

        // Output step, i, r, q, g(r), c(r), S(q)
        fprintf(fp, "%i \t %lf \t %lf \t %lf \t %lf \t %lf \n",(i+1),r, q, g[i], c[i],S[i]);
    }
    fclose(fp);

}


int main() {

    // It loops over all packingfraction that are asked (0.1,0.2,0.3,0.4)
    // 0.1 uses the guess of gamma(r) that is given above in the variables
    // The others use the gamma(r) that is generated in the previous packing fraction
    for(packing_fraction=0.1; packing_fraction<=0.4; packing_fraction+=0.1)
    {

        double V=(M_PI*pow(diam,3.0)*N)/(6*packing_fraction);   // Packing fraction defines volume of system
        double rho=N/V;                                         // rho defined by packing fraction
        double r;
        double q;
        double dq=M_PI/(points*dr);                             // drdq = π/L (see lecture notes)
        double cr[points];
        double gammq[points];

        std::cout<< "\n----------------";
        std::cout<< "\nRun for packing faction: "<< packing_fraction << " started.";

        for(int run=0; run<runs;run++)
        {

            // [Step II]
            // Calculating c(r)
            for(int i=0;i<points;i++){

                // Saving old value of gamma
                gamm_old[i]=gamm[i];

                r=(i)*dr;                                       // Point on the left border of the grid
                if(r<diam)
                {
                    c[i]=-gamm[i]-1;
                }
                else if(r>=diam)
                {
                    c[i]=0.0;
                }

                //std::cout<< "c_r[" <<i << "," << r << "] = \t" << c[i] << "\n";
            }

            // [Step III]
            // Fourier transform c(r) to c(q)

            for(int i=0;i<points;i++)
            {
                r=(i)*dr;
                cr[i]=c[i]*r;
            }

            // Sinfft.h package that computes the Fourier Transform
            sinft(cr,points);


            for(int i=0;i<points;i++)
            {
                // Important that the pre-factors of the Fourier transform are added
                q=(i)*dq;                                           // Point on the right border of the grid

                if(i!=0){c[i]=(4*M_PI*dr)/(q) * cr[i];}
                else{c[i]=0;}                                       // First value of the fourier transform is 0

            }


            // [Step IV]
            // Calculate gamma^hat = gamma(q)
            for(int i=0;i<points;i++)
            {
                gamm[i]=rho*pow(c[i],2.0)/(1-rho*c[i]);
            }

            // [Step V]
            // Invert gamma(q) to gamma(r) using same Fourier transform
            for(int i=0;i<points;i++)
            {
                q=(i)*dq;
                gammq[i]=gamm[i]*q;
            }

            // Sinfft.h package that computes the Fourier Transform
            sinft(gammq,points);

            for(int i=0;i<points;i++)
            {
                // Important that the pre-factors of the Fourier transform are added
                r=(i)*dr;

                if(i!=0){gamm[i]=dq/(2*pow(M_PI,2.0)*r)*gammq[i];}
                else{gamm[i]=0;}                                    // First value of the fourier transform is 0

                //std::cout<< "gamma_new[" <<i << "," << r << "] = \t" << gamm[i] << "\n";
            }


            // Next step is to calculate if we have a new solution gamma
            double gamma_difference = 0;
            for(int i=0; i<points; i++)
            {
                // Broyles mixing
                gamm_old[i]=alpha*gamm[i]+(1-alpha)*gamm_old[i];

                // Calculating difference
                gamma_difference+=abs(gamm[i]-gamm_old[i]);

            }

            // If the absolute sum of the gamma(r) is smaller than TOL, than this solution for gamma(r) is accepted and exported
            if(gamma_difference<=TOL)
            {

                // Export g(r), c(r), S(q) over all used r and q space
                print_file(dr,dq,rho,packing_fraction,run);

                std::cout<< "\nSolution found!";
                std::cout<<"\nSum of absolute gamma difference = " << gamma_difference;
                std::cout<< "\nIn "<< run <<" runs." << "\n";

                break;
            }

            // Look if it is not diverging, so print after 1000 runs what the gamma_difference is
            if(run % 1000==0)
            {
                if(run!=0)
                {
                    std::cout<< "\nRun: "<< run << ", Gamma difference = " << gamma_difference;
                }
            }

        }

        //std::cout<< "\nDid NOT found solution with gamma difference smaller than TOL! Increase runs." << "\n";

        //return 1;
    }
    std::cout<< "\nProgram stopped" << "\n";
    return 1;


}