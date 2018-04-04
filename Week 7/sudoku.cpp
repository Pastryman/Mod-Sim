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
int sudoku[size][size];
int sudoku_old[size][size];
int sudoku_initial[size][size];
int score;
const char*  init_filename = "sudoku.dat";

/* Parameters */
const double T = 0.5;


void read_sudoku(void){
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

int get_score(int sud[9][9]){
    // Calculates score of a sudoku
}

void attempt_change()
{

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

    read_sudoku();

    print_sudoku(sudoku);

    first_fill();

    print_sudoku(sudoku);

    //print_sudoku(sudoku_initial);

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