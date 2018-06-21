#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <chrono>  // for high_resolution_clock

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Measurements constants
const double eta = pow(5,-2); //[10^-4,0.1]     // Packing fraction
const double sigma = 0.1; //[0.001,0.10]        // Charge density (nm^-2)
const int M = 1000000;                          // Cells

// Initializing measurements for while loop
double ka = 0.1; //[0.1,10]                     // k^-1= debye length (nm)
double ka_max = 10;
double d_ka = 0.005;

// Measurement variables
double phi0_start = 0;                          // First guess phi0
double phi0_max = 1000000;                      // Max values for phi0 in guessing
double delta_start = 50;                        // dphi0 steps, firstly
double delta_min = 0.0000000000001;             // To prevent to get stuck in while loop
const double d_phi_max=0.000000001;             // Acceptance value for phi(M)-phi(M-1)<delta_min

// Variable declaration
double phi_0;
double phi_1;
double phi_0_temp;
double phi_1_temp;

// System constants
double a;   // Colloid size (nm)
double R;   // Cell radius (nm)
double k;   // k^-1 = Debye lentgh (nm)
double Z;   // charge
double lb;  // Bjerumm length (nm)
double h;   // Cell size in (nm)

double phi(){

    // Intializing local variables
    phi_0_temp=phi_0;
    phi_1_temp=phi_1;

    double rj;
    double term1, term2, term3;

    for (int j = 2; j <= M; ++j) {

        // Define rj, and terms of the recursive formula
        rj = a + double(j)*h;
        term1 = 2*rj/(rj+h)*phi_1_temp;
        term2 = ((h-rj)/(h+rj))*phi_0_temp;
        term3 = pow(k*h,2.)*(rj/(rj+h))*sinh(phi_1_temp);

        // Recursion for next iteration
        phi_0_temp=phi_1_temp;
        phi_1_temp=term1+term2+term3;

        // This is build in to find the optimal solution for phi0
        // If the solution is -inf, we have a to low phi0 guess
        // If the solution is +inf, we have a to high phi0 guess
        if (phi_1_temp==INFINITY){
            return 1;
        }
        else if (phi_1_temp==-INFINITY){
            return -1;
        }
    }

    // If not infinity than for every step it has a solution to phi(ri)
    return 0;
}

int main() {

    // Intializing file for exporting data
    char buffer[128];
    sprintf(buffer, "phi_eta%.5f_kaM_sigma%.3f.dat", eta, sigma);
    FILE *fp = fopen(buffer, "w");

    // Variables in the file
    fprintf(fp, "#phi(0)\tphi(R)\teta\tka\tsigma\tdelta");

    // For this measurement loop over all the wanted values for ka and guess what phi0 has to be
    while (ka < ka_max) {

        // Give values to variables for every new loop to look for phi0
        a = 1;                                  // Colloid size (nm)
        R = a * pow(eta, -1. / 3.);             // Cell radius (nm)
        k = ka / a;                             // (Debye length)^-1 (nm)
        Z = sigma * 4 * M_PI * pow(a, 2.);      // Charge, sigma (nm^-2)
        lb = 0.72;                              // Bjerrum length (nm)
        h = (R - a) / M;                        // Cell size (nm)

        // Start guess to phi0 and phi1
        phi_0 = phi0_start;
        phi_1 = phi_0 - (Z * lb * h / pow(a, 2));


        int sign = -1;                  // Added to know if it is in -inf of +inf singularity
        double phi_diff = INFINITY;     // phi_diff=phi1-phi0 for knowing what phi(M)-phi(M-1)=0
        int counter = 0;                // Added to count how much iteration we have for finding phi0
        double delta = delta_start;     // first dphi0 for searching the phi0 for which phi(M)-phi(M-1)=0

        // If phi_diff > d_phi_max we have found our phi0
        while (abs(phi_diff) > d_phi_max) {

            // Result for this phiM, phi(M-1) for this phi0
            double result = phi();

            // If it has a solution for the whole space of R
            if (result == 0) {
                phi_diff = phi_1_temp - phi_0_temp;

                // If solution is smaller than acceptance, a solution is found
                if (abs(phi_diff) < d_phi_max) {
                    std::cout << "\n\n Found the optimal phi(0)!!!!!";
                    printf("\n phi(0) = %.15f", phi_0);
                    printf("\n phi(M)-phi(M-1) = %.15f", phi_diff);
                    break; // Solution found -> break while loop

                    // If no solution is found look closer around phi0 to find solution
                } else if (phi_diff > 0) {
                    result = 1;
                } else if (phi_diff < 0) {
                    result = -1;
                }
            }

            // If in -inf search for higher phi0
            if (result == -1) {
                if (sign != result) {
                    sign *= -1;
                    delta *= 0.5;
                    printf("\nSIGN FLIP. delta is now %.15f", delta);
                }
                phi_0 = phi_0 + delta;                      // New guess for phi0
                phi_1 = phi_0 - (Z * lb * h / pow(a, 2));   // New guess for phi1
                // If we look to close but cannot find solution, cut the while loop off and give phi0=-1, which we now is incorrect
                if (delta<delta_min) {
                    phi_0 = -1;
                    break; // no solution found -> break while loop
                }
                continue;
                // If in +inf search for smaller phi0
            } else if (result == 1) {
                if (sign != result) {
                    sign *= -1;
                    delta *= 0.5;
                    printf("\nSIGN FLIP. delta is now %.15f", delta);
                }
                phi_0 = phi_0 - delta;                      // New guess for phi0
                phi_1 = phi_0 - (Z * lb * h / pow(a, 2));   // New guess for phi0
                // If we look to close but cannot find solution, cut the while loop off and give phi0=-1, which we now is incorr
                if (delta<delta_min) {
                    phi_0 = -1;
                    break; // no solution found -> break while loop
                }
                continue;
            }
            counter += 1;

        }

        // Print some values
        std::cout << "\nh = " << h;
        std::cout << "\nM = " << M;
        std::cout << "\nka = " << ka;

        // Print values to file
        fprintf(fp, "\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf", phi_0, phi_1_temp, eta, ka, sigma, delta);

        // Search for phi0 for next ka
        ka+=d_ka;

    }
    fclose(fp);
    return 1;
}