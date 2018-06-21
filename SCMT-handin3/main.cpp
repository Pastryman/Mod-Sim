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



const double eta = pow(5,-2); //[10^-4,0.1] // packing fraction
const double sigma = 0.1; //[0.001,0.10] charge density
const int M = 10000;


double ka = 0.1; //[0.1,10] k^-1= debye length
double ka_max = 10;
double d_ka = 0.005;

// Measurement variables
double phi0_start = 0;
double phi0_max = 1000000;
double delta_start = 50;
double delta_min = 0.0000000000001;
const double d_phi_max=0.000000001;

// Variable declaration
double phi_0;
double phi_1;
double phi_0_temp;
double phi_1_temp;

// System constants
double a; // Colloid size (nm)
double R;
double k;
double Z;
double lb;
double h;

double phi(){


    phi_0_temp=phi_0;
    phi_1_temp=phi_1;

    double rj;
    double term1, term2, term3;

    for (int j = 2; j <= M; ++j) {
//        std::cout << "\nj = " << j;

        rj = a + double(j)*h;
        term1 = 2*rj/(rj+h)*phi_1_temp;
        term2 = ((h-rj)/(h+rj))*phi_0_temp;
        term3 = pow(k*h,2.)*(rj/(rj+h))*sinh(phi_1_temp);

//        std::cout << "\nterm1,term2,term3 = " << term1 << "," << term2 << "," << term3;

        phi_0_temp=phi_1_temp;
        phi_1_temp=term1+term2+term3;


//        std::cout << "\nphi new = " << phi_1_temp;
//        std::cout << "\nrj = " << rj << "\n";

        if (phi_1_temp==INFINITY){
            return 1;
        }
        else if (phi_1_temp==-INFINITY){
            return -1;
        }
    }
    return 0;
}

int main() {

    char buffer[128];
//    sprintf(buffer, "phi_eta%.5f_ka%.1f_sigma%.3f.dat", eta, ka, sigma);
    sprintf(buffer, "phi_eta%.5f_kaM_sigma%.3f.dat", eta, sigma);
    FILE *fp = fopen(buffer, "w");
    fprintf(fp, "#phi(0)\tphi(R)\teta\tka\tsigma\tdelta");

    while (ka < ka_max) {

        // System constants
        a = 1; // Colloid size (nm)
        R = a * pow(eta, -1. / 3.);
        k = ka / a;
        Z = sigma * 4 * M_PI * pow(a, 2.);
        lb = 0.72;
        h = (R - a) / M;

        phi_0 = phi0_start;
        phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
        int sign = -1;
        double phi_diff = INFINITY;
        int counter = 0;
        double delta = delta_start;

        while (abs(phi_diff) > d_phi_max) {

            //        std::cout << "\nphi0 = " << phi_0 << "\tphi1 = " << phi_1;

            double result = phi();

            //        std::cout << "\nphi function returns " << result << "\n";

            if (result == 0) {
                phi_diff = phi_1_temp - phi_0_temp;
                if (abs(phi_diff) < d_phi_max) {
                    std::cout << "\n\n Found the optimal phi(0)!!!!!";
                    printf("\n phi(0) = %.15f", phi_0);
                    printf("\n phi(M)-phi(M-1) = %.15f", phi_diff);
                    break;
                } else if (phi_diff > 0) {
                    result = 1;
                } else if (phi_diff < 0) {
                    result = -1;
                }
            }

            if (result == -1) {
                if (sign != result) {
                    sign *= -1;
                    delta *= 0.5;
                    printf("\nSIGN FLIP. delta is now %.15f", delta);
                }
                phi_0 = phi_0 + delta;
                phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
                if (delta<delta_min) {
                    phi_0 = -1;
                    break;
                }
                continue;
            } else if (result == 1) {
                if (sign != result) {
                    sign *= -1;
                    delta *= 0.5;
                    printf("\nSIGN FLIP. delta is now %.15f", delta);
                }
                phi_0 = phi_0 - delta;
                phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
                if (delta<delta_min) {
                    phi_0 = -1;
                    break;
                }
                continue;
            }
            counter += 1;

        }

        std::cout << "\nh = " << h;
        std::cout << "\nM = " << M;

        fprintf(fp, "\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf", phi_0, phi_1_temp, eta, ka, sigma, delta);

        std::cout << "\nka = " << ka;

        ka+=d_ka;

    }
    fclose(fp);
    return 1;
}