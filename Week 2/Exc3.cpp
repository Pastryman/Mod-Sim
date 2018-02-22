/*  Mod&Sim 2018 week 2, Exercise 3
 *
 *  Authors: Tom Niessen & Floris Jonkman
 * */

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

ofstream myfile;

/*Simulation variables*/
const float latticeSize = 2;    // Lattice spacing size of the unit cell
                                // NOTE: latticeSize must be larger than sqrt(2).
const float particleSize = 1;   // Always take the maximum particle size
const int NParticles = 1000;    // The amount of  particles we want to generate
                                // NOTE: this number is later corrected to perfectly fit in a cube

// Method: fcc
// Input: float dx,dy,dz: displacement of the lattice cell with respect to the origin
//                          if dx=6 it means we are in the 6th cell in the x direction
// Description: This method generates the coordinates of the particles in a cell from the position of the cell
//              This is fcc, so we have 4 particle per cell
//              The coords of the position of each 4 particles are then written to a file
void fcc(float dx, float dy, float dz)
{
    // Generate the coordinates of the particle in the origin of the cell
    float x = dx*latticeSize;
    float y = dy*latticeSize;
    float z = dz*latticeSize;
    // The space between the other particles and the origin of the cell
    float spacing = float(0.5)*latticeSize;
    // Write the positions of the particles to a file
    myfile << x << " " << y << " "<< z << " " << particleSize << "\n";
    myfile << x + spacing << " " << y + spacing << " "<< z << " " << particleSize << "\n";
    myfile << x + spacing << " " << y << " "<< z + spacing << " " << particleSize << "\n";
    myfile << x << " " << y + spacing << " "<< z+spacing << " " << particleSize << "\n";
}

int main() {
    // The number of cells we will draw
    // This number is rounded down to the first int, so we have only cells with 4 particles
    int NCells = NParticles/4;
    // The amount of cells in each direction
    // The number is rounded down to get the highest number cells that lead to a square box
    int Ndim = int(pow(NCells,1./3.));

    // Open the file
    myfile.open("xyz.dat");

    // Line 1 of the file: The total amount of particles
    myfile << pow(Ndim,3)*4 << "\n";
    // Lines 2-4: The dimensions of the (square) box
    myfile << 0 << " " << Ndim*latticeSize << "\n";
    myfile << 0 << " " << Ndim*latticeSize << "\n";
    myfile << 0 << " " << Ndim*latticeSize << "\n";

    // loop of all cells in the x direction
    for (int i = 0; i<Ndim; i++)
    {
        // loop of all cells in the y direction
        for (int j = 0; j<Ndim; j++)
        {
            // loop of all cells in the z direction
            for (int k = 0; k<Ndim; k++)
            {
                // Generate the particles in the current cell (and write to file)
                fcc(i,j,k);
            }
        }
    }

    myfile.close();
    return 0;
}