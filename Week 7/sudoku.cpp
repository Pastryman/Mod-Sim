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

/* Initializing sudoku*/
const int size      = 9;
int blocksize;
int sudoku[size][size];
int sudoku_initial[size][size];
int current_score;
const char*  init_filename = "sudoku.dat";

/* Parameters */
double T = 0.5;
const int output_steps = 1000;
const int mc_steps = 1000;
const int max_steps = 400000;

void read_sudoku(){
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
            if (sudoku[r][c]!=0) {sudoku_initial[r][c]=1;}
        }
    }
    fclose(fp);
}

void print_sudoku(int sud[size][size]){
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
    // High score is bad
    // Low score is good
    int score = 0;
    int count;

    for (int c = 0; c<size; c++){
        for (int val = 1; val<size+1; val++){
            count = 0;
            for (int r = 0; r<size; r++){
                if (sudoku[r][c] == val){count++;}
            }
            if (count == 1) {score--;}
        }
    }

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
    /*
     * 1. Pak een random rij
     * 2. Pak 2 random getallen
     * 3. Verwissel de getallen
     * 4. Nieuwe score berekenen
     * 5. Scores vergelijken (accepted?)
     * 6. Boltzmann factor (accepted?)
     * 7. Verandering wel/niet toepassen
     */

    // Pick a random row
    int row = int(dsfmt_genrand()*9.);
    // Pick two random positions (columns in that row)
    int pos1=0;
    int pos2=0;
    bool accept = false;
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
    // Fill first instance of sudoku
    // Method
    // The sudoku is filled stating that every row will contain the values from 1-size (9), so we loop over
    // every row in order to get a 'first fill' for the sudoku
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
            while(inrow==false)
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

    size_t seed = time(NULL);
    dsfmt_seed(seed);
    blocksize = int(sqrt(double(size)));
    std::cout<<"\n blocksize = " << blocksize << "\n";

    read_sudoku();

    print_sudoku(sudoku);
    print_sudoku(sudoku_initial);

    first_fill();
    current_score = get_score();
    print_sudoku(sudoku);

    int step = 0;
    int accepted = 0;
    while (step < max_steps)
    {
        int n;
        for (n = 0; n < mc_steps; n++){
            accepted += attempt_change();
            if (current_score == -(size*size*2))
            {
                std::cout << "SOLUTION FOUND\n";
                print_sudoku(sudoku);
                return 1;
            }
        }

        if (step%output_steps ==0)
        {
            std::cout << "After " << step << "steps, the current sudoku is:\n";
            print_sudoku(sudoku);
            std::cout << "Current score is " << current_score << "\n";
            double acceptanceRatio = double(accepted) / (double(mc_steps) * double(output_steps));
            printf("Step %d. acceptance: %lf. T: %lf\n", step, acceptanceRatio,T);
            accepted = 0;
        }

        step++;
        T*=0.99999;
    }

    // Sudoku vullen met random getallen, wel zodat elke rij (keuze) kloppend is

    // Iteraties doen totdat we een oplossing hebben
    // In elke iteratie wordt T verlaagd
    // Iteratie bestaat uit:
    //      1. Pak 2 getallen in een random rij en verwissel ze
    //      2. Bereken oude en nieuwe score en test old<new
    //      3. Test de Boltzmann factor
    //      4. Sta wissel toe of wissel terug
    //      5. Trek een bak
    // Stap 1-5 zouden in 1 methode kunnen, is wel zo mooi
    //      6. Verlaag T
    //      7. Check of we de oplossing hebben

    // Bedenk wat we willen printen tussendoor en wanneer
    // Elke 1000 iteraties: score, sudoku en acceptance rate

    return 1;
}