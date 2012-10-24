/**
 * Project: Matrix Multiplication
 * Student: Tri Nguyen, Zang Ze
 * Notes: We took helps and advices from Professor Robert Van De Geijn.
 */

/**
 *  Define macros as a neat way for indexing
 *  Accessing elements in column order
 */
#define A(x, y) A[ (y) * N + (x)]
#define B(x, y) B[ (y) * N + (x)]
#define C(x, y) C[ (y) * N + (x)]

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <time.h>

#include "matMult.h"
#include "partA.h"
#include "partB.h"
#include "partC.h"
#include "partD.h"
#include "competition.h"

/*MMM function declarations */
//part A
void matmult_ikj_a(double* A, double* B, double* C, unsigned N);
void matmult_jik_a(double* A, double* B, double* C, unsigned N);
//part B
void matmult_ikj_b_1_4(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_4(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_4_1(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_8(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_8_1(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_4_4(double* A, double* B, double* C, unsigned N);
//part C
void matmult_jik_c(double* A, double* B, double* C, unsigned N, unsigned NB);
//part D
void matmult_jik_d(double* A, double* B, double* C, unsigned N);

/* papi and matrix function declarations*/
void init_papi();
int begin_papi(int Event);
long_long end_papi(int EventSet);
void printMatrix(const double* X, unsigned N);
void initialize(double* const X, const double VAL);
void flushCache(double* matrix);
double* alloc(int SIZE);


/* Enable to separate partB and partC */
void partA(double* A, double* B, double* C, unsigned N, unsigned TYPE);
void partB(double* A, double* B, double* C, unsigned N, unsigned TYPE);
void partC(double* A, double* B, double* c, unsigned N, unsigned NB);
void partD(double* A, double* B, double* c, unsigned N, unsigned NB);

// Final Version
//void matmult(double* A, double* B, double* C, unsigned N);
void matmult(const double* A, const double* B, double* C, unsigned N);

/*Global variable*/
const int EVENT = PAPI_FP_OPS;

int main(int argc, char** argv) {
    const unsigned N = atoi(argv[1]);
    const unsigned T = atoi(argv[2]);
    const char*  part = argv[3];

    double* A = alloc(N);
    double* B = alloc(N);
    double* C = alloc(N);

    /*Initiallize arrays */
    initialize(A, 1, N);
    initialize(B, 2, N);
    initialize(C, 0, N);

    init_papi();
    /*Start clocking*/
    int eventSet = begin_papi(EVENT);
    long long ret;
    int start = clock();

    if (strcmp(part, "A") == 0)
        partA(A, B, C, N, T);

    if (strcmp(part, "B") == 0)
        partB(A, B, C, N, T);

    if (strcmp(part, "C") == 0)
        partC(A, B, C, N, T);

    if (strcmp(part, "D") == 0)
        partD(A, B, C, N, T);

    if (strcmp(part, "F") == 0)
        matmult(A, B, C, N);

    if (argc > 4){
        const char*  debug = argv[4];
        if (strcmp(debug, "-v") == 0)
            printMatrix(C, N);
    }

    /* Stop  clocking    */
    int stop = clock();
    ret = end_papi(eventSet);
    int elapsed_time = stop - start;
    std::cout << "\t" << N << "\t" << ret << "\t\t "
            << ((float) elapsed_time) / CLOCKS_PER_SEC << std::endl;

    flushCache(A);
    flushCache(B);
    flushCache(C);

    return 0;
}

void partC(double* A, double* B, double* C, unsigned N, unsigned NB) {
    //  NB from 16 -> 23
    // Mini kernel jik_b_1_4
    matmult_jik_c(A, B, C, N, NB);
    std::cout << NB;
}

void partD(double* A, double* B, double* C, unsigned N, unsigned TYPE) {
    const unsigned NB = 64;
    matmult_jik_d(A, B, C, N, NB);
}

void partA(double* A, double* B, double* C, unsigned N, unsigned TYPE) {
    const std::string TYPE_STRING[5] = {
            " ",
            "ikj_a", "jik_a",
    };

    switch (TYPE) {
        case 1:  //part a
            matmult_ikj_a(A, B, C, N);
            break;
        case 2:  //part a
            matmult_jik_a(A, B, C, N);
            break;
    }
    std::cout << TYPE_STRING[ TYPE ] << std::endl;
}
void partB(double* A, double* B, double* C, unsigned N, unsigned TYPE) {
    const std::string TYPE_STRING[20] = {
            " ",
            "ikj_b_1_4", "jik_b_1_4",
            "jik_b_4_1", "jik_b_1_8",
            "jik_b_8_1", "jik_b_4_4",
    };

    switch (TYPE) {
        case 1:  //part b
            matmult_ikj_b_1_4(A, B, C, N);
            break;

        case 2:  //part b
            matmult_jik_b_1_4(A, B, C, N);
            break;

        case 3:  //part b
            matmult_jik_b_4_1(A, B, C, N);
            break;

        case 4:  //part b
            matmult_jik_b_1_8(A, B, C, N);
            break;

        case 5:  //part b
            matmult_jik_b_8_1(A, B, C, N);
            break;

        case 6:  //part b
            matmult_jik_b_4_4(A, B, C, N);
            break;
    }
    std::cout << TYPE_STRING[ TYPE ] << std::endl;
}
