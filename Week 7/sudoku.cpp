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
}

int main() {

    read_sudoku();

    print_sudoku(sudoku);
    print_sudoku(sudoku_initial);

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