#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <chrono>  // for high_resolution_clock

// Measurement variables
double phi0_start = 0;
double phi0_max = 100;
double delta = 50;
const double d_phi_max=0.0000001;

// System constants
const double R = 1;
const double a = 0.01;
const int M = 100;
const double k = 1;
const int Z = 1;
const double lb = 0.72;
const double h = (R-a)/M;

// Variable declaration
double phi_0;
double phi_1;
double phi_0_temp;
double phi_1_temp;

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

    std::cout << "\nh = " << h;

    //phi_0 = 84.94;
    phi_0 = phi0_start;
    phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
    int sign = -1;
    double phi_diff = INFINITY;

    while(abs(phi_diff)>d_phi_max) {

//        std::cout << "\nphi0 = " << phi_0 << "\tphi1 = " << phi_1;

        double result = phi();

//        std::cout << "\nphi function returns " << result << "\n";

        if (result == 0){
            phi_diff = phi_1_temp - phi_0_temp;
            if (abs(phi_diff)<d_phi_max){
                std::cout << "\n\n Found the optimal phi(0)!!!!!";
                printf("\n phi(0) = %.10f",phi_0);
                printf("\n phi(M)-phi(M-1) = %.10f",phi_diff);
                return 0;
            }
            else if (phi_diff>0){
                result = 1;
            }
            else if (phi_diff<0){
                result = -1;
            }
        }

        if (result==-1){
            if (sign!=result){
                sign*=-1;
                delta*=0.5;
                printf("\nSIGN FLIP. delta is now %.10f",delta);
            }
            phi_0= phi_0 + delta;
            phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
            continue;
        }

        else if (result == 1){
            if (sign!=result){
                sign*=-1;
                delta*=0.5;
                printf("\nSIGN FLIP. delta is now %.10f",delta);
            }
            phi_0= phi_0 - delta;
            phi_1 = phi_0 - (Z * lb * h / pow(a, 2));
            continue;
        }
    }
    return 0;
}