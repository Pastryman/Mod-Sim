/*  Mod&Sim 2018 week 2, Exercise 1
 *
 *  Authors: Tom Niessen & Floris Jonkman
 * */

#include <iostream>
#include <fstream>
using namespace std;

// Initialize a stream to our file
ofstream myfile;

/*Simulation variables*/
const float latticeSize = 1;  // Lattice spacing (size of the unit cell)
const float particleSize = 1; // Size of a particle
const int Nx = 10; // Amount of particles in the x direction
const int Ny = 10; // Amount of particles in the y direction
const int Nz = 10; // Amount of particles in the z direction

// Method: Cubic
// Input: float dx,dy,dz: displacement of the lattice cell with respect to the origin
// Description: This method generates the coordinates of the particle from the coordinates of the cell
//              This is cubic, so we only have 1 particle per cell
//              The coords are then written to a file
void Cubic(float dx, float dy, float dz)
{
    // Calculation of the coordinates of the particle in the cell
    float x = dx*latticeSize;
    float y = dy*latticeSize;
    float z = dz*latticeSize;
    // Writing the coordinates to our file (space separated)
    myfile << x << " " << y << " " << z << " " << particleSize << "\n";
}

int main() {
    // Open our file
    myfile.open("xyz.dat");
    // Line 1 of the file: The total amount of particles
    myfile << Nx*Ny*Nz << "\n";
    // Lines 2-4: The dimensions of the box
    myfile << 0 << " " << Nx*latticeSize << "\n";
    myfile << 0 << " " << Ny*latticeSize << "\n";
    myfile << 0 << " " << Nz*latticeSize << "\n";

    // Loop over every cell in our box
    // The x direction
    for (int i = 0; i<Nx; i++)
    {
        //The y direction
        for (int j = 0; j<Ny; j++)
        {
            // The z direction
            for (int k = 0; k<Nz; k++)
            {
                // Generate the particle and write to file
                Cubic(i,j,k);
            }
        }
    }

    // Close the file
    myfile.close();
    return 0;
    }