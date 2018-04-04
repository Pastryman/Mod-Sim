#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "mt19937.h"

/* Initializing sudoku*/
const int size      = 9;
int sudoku[size][size];
const char*  sodoku_filename = "sudoku.csv";

/* Parameters */
const double T = 0.5;


void read_sudoku(void){
    FILE* fp = fopen(sodoku_filename, "r");
    int r,c;
    for(r = 0; r < size; ++r){
        for(c = 0; c < size; ++c){
            fscanf(fp, "%d;", &sudoku[r][c]);
            if(c==size-1){
                fscanf(fp, "%d\n", &sudoku[r][c]);
            }
        }
    }
    fclose(fp);
}

int main() {
    read_sudoku();

    for(int r = 0; r < size; ++r){
        for(int c = 0; c < size; ++c){
            if(c % 3 == 0){std::cout << "|";}
            std::cout << sudoku[r][c] << "\t";
            if(c==size-1){
                std::cout << sudoku[r][c] << "\n";;
            }
        }
    }
}