#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "matMult.h"

/* define macros as a neat way for index*/
#define A(x, y) A[(x) * N + (y)]
#define B(x, y) B[(x) * N + (y)]
#define C(x, y) C[(x) * N + (y)]


/*MMM function declarations */
void matmult_ikj_a(double* A, double* B, double* C) ;
void matmult_jik_a(double* A, double* B, double* C) ;

void matmult_ikj_b(double* A, double* B, double* C) ;
void matmult_jik_b(double* A, double* B, double* C) ;
void matmult_jik_b_2(double* A, double* B, double* C);

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
            matmult_ikj_b(A, B, C);
            break;
        case 4:  //part b
            matmult_jik_b(A, B, C);
            break;
        case 5:  //part b
            matmult_jik_b_2(A, B, C);
            break;

    }

    /* Stop  clocking    */
    ret = end_papi(eventSet);
    std::cout << TYPE << "\t" << N << "\t"
              << ret  << std::endl;
    printMatrix(C);

    flushCache(A);
    flushCache(B);
    flushCache(C);

    return 0;
}

/**
 * Simple Matrix Multiplication.
 * For Part A
 */
void matmult_ikj_a(double* A, double* B, double* C) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned k = 0; k < N; ++k) {
            for (unsigned j = 0; j < N; ++j) {
                C(i,j) += A(i,k) * B(k,j);
            }
        }
    }
}

void matmult_jik_a(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned k = 0; k < N; ++k) {
                C(i,j) += A(i,k) * B(k,j);
            }
        }
    }
}

/**
 * Mini Kernel MMM for part B
 * Unrolling loop j, i
 * Regiser Blocking
 */
void matmult_ikj_b(double* A, double* B, double* C) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned k = 0; k < N; ++k) {
            register double a0 = A(i,k);
            for (unsigned j = 0; j < N; j+=4) {
                double* c = &C(i,j);
                double* b = &B(k,j);
                c[0] += a0 * b[0];
                c[1] += a0 * b[1];
                c[2] += a0 * b[2];
                c[3] += a0 * b[3];
            }
        }
    }
}

/**
 *  NU = 1; MU = 4;
 */
void matmult_jik_b(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i+=5) {

            register double c0, c1, c2, c3 ,c4;
            c0 = c1 = c2 = c3 = c4 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                b0  = B(k,   j);
                c0 += A(i,   k) * b0;
                c1 += A(i+1, k) * b0;
                c2 += A(i+2, k) * b0;
                c3 += A(i+3, k) * b0;
                c4 += A(i+4, k) * b0;
            }
            C(i,   j) += c0;
            C(i+1, j) += c1;
            C(i+2, j) += c2;
            C(i+3, j) += c3;
            C(i+4, j) += c4;
        }
    }
}

/**
 * NU = 2; MU = 4;
*/
void matmult_jik_b_2(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; j+=2) {
        for (unsigned i = 0; i < N; i+=4) {

            register double c00, c01, c02, c03,
                            c10, c11, c12, c13;
            c00 = c01 = c02 = c03 = 0;
            c10 = c11 = c12 = c13 = 0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0, b1;
                b0  = B(k,   j);
                b1  = B(k, j+1);
                c00 += A(i,   k) * b0;
                c01 += A(i+1, k) * b0;
                c02 += A(i+2, k) * b0;
                c03 += A(i+3, k) * b0;

                c10 += A(i,   k) * b1;
                c11 += A(i+1, k) * b1;
                c12 += A(i+2, k) * b1;
                c13 += A(i+3, k) * b1;
            }
            C(i,   j) += c00;
            C(i+1, j) += c01;
            C(i+2, j) += c02;
            C(i+3, j) += c03;

            C(i,   j+1) += c10;
            C(i+1, j+1) += c11;
            C(i+2, j+1) += c12;
            C(i+3, j+1) += c13;
        }
    }
}
