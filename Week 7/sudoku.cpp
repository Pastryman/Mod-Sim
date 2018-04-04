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
        }
    }
    fclose(fp);
}

void print_sudoku(void){
    for(int r = 0; r < size; ++r){
        if(r*r % size == 0){std::cout << "---------------------\n"; }
        for(int c = 0; c < size; ++c){
            if(c*c % size == 0){std::cout << "|";}
            if(c==size-1){
                std::cout << sudoku[r][c] << "|\n";;
            }
            else{
                std::cout << sudoku[r][c] << " ";
            }
        }
    }
}


int main() {

    read_sudoku();

    print_sudoku();

    return 1;
}