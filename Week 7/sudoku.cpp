#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "mt19937.h"
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

/* Initializing Sudoku */
const int size      = 9;
int blocksize;
int sudoku[size][size];
int sudoku_initial[size][size];
int current_score;
const char*  init_filename = "sudoku_extreme.dat";

/* Parameters */
double T                = 0.5;
double deltaMin         = 0.99999;
double deltaPlus        = 1.001;
int stuck_value         = 10;
const int output_steps  = 1000;
const int mc_steps      = 1000;
const int max_steps     = 400000;


void read_sudoku(){
    // Method that reads the Sudoku
    // File must be a .dat file seprated by '\t'
    // Vacant spots must be filled with 0

    FILE* fp = fopen(init_filename, "r");
    int r,c;
    for(r = 0; r < size; ++r){
        for(c = 0; c < size; ++c){
            if(c==size-1){
                fscanf(fp, "%d\n", &sudoku[r][c]);
            }
            else{
                fscanf(fp, "%d\t", &sudoku[r][c]);
            }
            // Also saves the the location of the initial (nonzero) values  in the Sudoku
            // These spots are not allowed to be swapped
            if (sudoku[r][c]!=0) {sudoku_initial[r][c]=1;}
        }
    }
    fclose(fp);
}

void print_sudoku(int sud[size][size]){
    // Method that prints the Sudoku
    // Code optimilized for a 9x9 Sudoku

    for(int r = 0; r < size; ++r){
        if(r*r % size == 0){std::cout << "---------------------\n"; }
        for(int c = 0; c < size; ++c){
            if(c*c % size == 0){std::cout << "|";}
            if(c==size-1){
                std::cout << sud[r][c] << "|\n";;
            }
            else{
                std::cout << sud[r][c] << " ";
            }
        }
    }
    std::cout << "---------------------\n";
}

int get_score(){
    // High score is bad (less minus)
    // Low score is good (more minus)
    int score = 0;
    int count;

    // Count all the unique values in the columns
    for (int c = 0; c<size; c++){
        for (int val = 1; val<size+1; val++){
            count = 0;
            for (int r = 0; r<size; r++){
                if (sudoku[r][c] == val){count++;}
            }
            if (count == 1) {score--;}
        }
    }

    // Count all the unique values in the boxes
    for (int c = 0; c < blocksize; c++) {
        for (int r = 0; r < blocksize; r++) {
            for (int val = 1; val < size + 1; val++) {
                count = 0;
                for (int i = 0; i < blocksize; i++) {
                    for (int j = 0; j < blocksize; j++) {
                        if (sudoku[(r * blocksize) + i][(c * blocksize) + j] == val) { count++; }
                    }
                }
                if (count == 1) { score--; }
            }
        }
    }
    return score;
}

int attempt_change()
{
    // Method that tries two swap two values
    // returns 1 (accepted) or 0 (rejected)
    //
    // 1 - Pick a random row
    // 2 - Pick two random integers
    // 3 - Swap integers in rwo
    // 4 - Calculate score of new solution
    // 5 - If better -> accepted
    // 6 - If worse -> Calculate Boltzmann factor -> accepted or rejected
    // (7) - Accept new solution

    // Pick a random row
    int row = int(dsfmt_genrand()*9.);
    // Pick two random positions (columns in that row)
    int pos1=0;
    int pos2=0;
    bool accept = false;

    // If one of the two selected positions is on a spot of a nonzero value in the initial Sudoku
    // Then this another spot
    while(!accept)
    {
        pos1 = int(dsfmt_genrand()*9.);
        if (!sudoku_initial[row][pos1]) {accept = true;}
    }
    accept = false;
    while(!accept)
    {
        pos2 = int(dsfmt_genrand()*9.);
        if (!sudoku_initial[row][pos2]) {accept = true;}
    }

    // Values of the two selected spots
    int nr1 = sudoku[row][pos1];
    int nr2 = sudoku[row][pos2];

    // Swap the numbers and calculate the new score
    sudoku[row][pos1] = nr2;
    sudoku[row][pos2] = nr1;
    int new_score = get_score();

    // If the new score is better, we accept immediately
    if (new_score < current_score)
    {
        current_score = new_score;
        return 1;
    }

    // Calculate the boltzmann factor
    double boltz = exp(double(current_score-new_score)/T);
    double rand = dsfmt_genrand();
    if (rand < boltz)
    {
        current_score = new_score;
        return 1;
    }

    // If we reach this code, the move is rejected
    // Change back the swap we made
    sudoku[row][pos1] = nr1;
    sudoku[row][pos2] = nr2;
    return 0;

}

void first_fill(){
    // Fill first instance of Sudoku
    //
    // Method
    // The sudoku is filled stating that every row will contain the values from 1-size (9), so we loop over
    // every row in order to get a 'first fill' for the Sudoku
    //
    // 1 - Loop over all the rows
    // 2 - Pick an integer 'n' ranging [1,size]
    // 3 - If 'n' in row                -> n+1
    // 4 - If 'n' not in row            -> pick a random spot in the row and check whether it is filled or not
    // 5 - If filled                    -> pick another random spot
    // 6 - If not filled                -> give value of 'n' to that spot -> n+1

    // Loop over all the rows
    for(int r=0; r<size; r++)
    {
        // Pick an integer 'n' ranging [1,size]
        int n=1;
        while(n<=size)
        {
            // Check if integer 'n' is already in the row
            bool inrow = false;
            for(int c=0;c<size;c++)
            {
                if(n==sudoku[r][c])
                {
                    inrow=true;
                    break;
                }
            }

            // If integer not in row random select spots until a spot is found where no value is located yet
            while(!inrow)
            {
                // Select random spot
                int spot= int(dsfmt_genrand()*size);

                // Check whether spot is filled or not
                if(sudoku[r][spot]==0)
                {
                    sudoku[r][spot]=n;
                    inrow=true;
                }
            }

            n++;
        }

    }
}

int main() {

    // Intializing seed
    size_t seed = time(NULL);
    dsfmt_seed(seed);

    // Calculate block size
    blocksize = int(sqrt(double(size)));

    // Read Sudoku
    read_sudoku();

    // Print initial Sudoku
    print_sudoku(sudoku);

    // Randomly fill Sudoku and calculate score
    first_fill();
    current_score = get_score();


    // Monte Carlo simulation
    int step = 0;
    int accepted = 0;
    int old_score = 0;
    int stuck = 0;
    while (step < max_steps)
    {
        int n;
        for (n = 0; n < mc_steps; n++){
            accepted += attempt_change();

            // If optimal solution is found, print and stop
            if (current_score == -(size*size*2))
            {
                std::cout << "\nSOLUTION FOUND\n";
                std::cout << "step = " << step << ". T = " << T << "\n";
                print_sudoku(sudoku);
                return 1;
            }
        }

        if (step%output_steps ==0)
        {
            // Output
            //std::cout << "After " << step << " steps, the current sudoku is:\n";
            //print_sudoku(sudoku);
            std::cout << "Current score is " << current_score << "\n";
            double acceptanceRatio = double(accepted) / (double(mc_steps) * double(output_steps));
            printf("Step %d. acceptance: %lf. T: %lf\n", step, acceptanceRatio,T);
            accepted = 0;
        }

        step++;
        T*=deltaMin;

        // If stuck in a local minimum T is increased
        if (current_score == old_score){stuck++;}
        else {stuck = 0;}
        if (stuck == stuck_value) {T*=deltaPlus; }
        old_score = current_score;
    }

    return 1;
}