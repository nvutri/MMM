/* define macros as a neat way for index*/
#define A(x, y) A[(x) * N + (y)]
#define B(x, y) B[(x) * N + (y)]
#define C(x, y) C[(x) * N + (y)]

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <time.h>

#include "matMult.h"
#include "partA.h"
#include "partB.h"

/*MMM function declarations */
void matmult_ikj_a(double* A, double* B, double* C) ;
void matmult_jik_a(double* A, double* B, double* C) ;

void matmult_ikj_b_1_4(double* A, double* B, double* C) ;
void matmult_jik_b_1_3(double* A, double* B, double* C) ;
void matmult_jik_b_1_6(double* A, double* B, double* C) ;
void matmult_jik_b_1_4(double* A, double* B, double* C) ;
void matmult_jik_b_1_5(double* A, double* B, double* C) ;
void matmult_jik_b_2_3(double* A, double* B, double* C) ;

/* papi and matrix function declarations*/
void init_papi();
int begin_papi(int Event);
long_long end_papi(int EventSet);
void printMatrix(double* X);
void initialize(double* const X, const double VAL);
void flushCache(double* matrix) ;
double* alloc(int SIZE);

/*Global variable*/
int N;
int TYPE;
const int EVENT = PAPI_FP_OPS;

int main(int argc, char** argv) {

    N = atoi(argv[1]);
    TYPE = atoi(argv[2]);
    /*Look for number of floating point computations*/


    /* Run the measurement*/
    double* A = alloc(N);
    double* B = alloc(N);
    double* C = alloc(N);

    /*Initiallize arrays */
    initialize(A, 1);
    initialize(B, 2);
    initialize(C, 0);

    init_papi();
    /*Start clocking*/
    int eventSet = begin_papi(EVENT);
    long long ret;

    switch (TYPE) {
        case 1:  //part a
            matmult_ikj_a(A, B, C);
            break;
        case 2:  //part a
            matmult_jik_a(A, B, C);
            break;

        case 3:  //part b
            matmult_ikj_b_1_4(A, B, C);
            break;

        case 4:  //part b
            matmult_jik_b_2_3(A, B, C);
            break;
        case 5:  //part b
            matmult_jik_b_1_4(A, B, C);
            break;

        case 6:  //part b
            matmult_jik_b_1_5(A, B, C);
            break;

        case 7:  //part b
            matmult_jik_b_1_6(A, B, C);
            break;

        case 8:  //part b
            matmult_jik_b_1_3(A, B, C);
            break;
    }

    /* Stop  clocking    */
    ret = end_papi(eventSet);
    int t = clock();
    std::cout << TYPE << "\t"
              <<    N << "\t"
              << ret  << "\t\t "
              << ((float)t)/CLOCKS_PER_SEC
              << std::endl;

//    printMatrix(C);

    flushCache(A);
    flushCache(B);
    flushCache(C);

    return 0;
}
