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
#include "partC.h"

/*MMM function declarations */
//part A
void matmult_ikj_a(double* A, double* B, double* C, unsigned N);
void matmult_jik_a(double* A, double* B, double* C, unsigned N);
//part B
void matmult_ikj_b_1_4(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_3(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_6(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_4(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_5(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_2_3(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_4_1(double* A, double* B, double* C, unsigned N);
void matmult_jik_b_1_8(double* A, double* B, double* C, unsigned N);
//part C
void matmult_jik_c(double* A, double* B, double* C, unsigned N, unsigned NB);

/* papi and matrix function declarations*/
void init_papi();
int begin_papi(int Event);
long_long end_papi(int EventSet);
void printMatrix(double* X, unsigned N);
void initialize(double* const X, const double VAL);
void flushCache(double* matrix);
double* alloc(int SIZE);


/* Enable to separate partB and partC */
void partB(double* A, double* B, double* C, unsigned N, unsigned TYPE);
void partC(double* A, double* B, double* c, unsigned N, unsigned NB);


/*Global variable*/
const int EVENT = PAPI_FP_OPS;

int main(int argc, char** argv) {
    const unsigned N = atoi(argv[1]);
    const unsigned T = atoi(argv[2]);
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

    partB(A, B, C, N, T);
//    partC(A, B, C, N, T);

    /* Stop  clocking    */
    ret = end_papi(eventSet);
    int t = clock();
    std::cout << "\t" << N << "\t" << ret << "\t\t "
            << ((float) t) / CLOCKS_PER_SEC << std::endl;

    flushCache(A);
    flushCache(B);
    flushCache(C);

    return 0;
}

void partC(double* A, double* B, double* C, unsigned N, unsigned NB) {
    //  NB from 16 -> 23
    matmult_jik_c(A, B, C, N, NB);
    std::cout << NB;
}

void partB(double* A, double* B, double* C, unsigned N, unsigned TYPE) {
    const std::string TYPE_STRING[11] = {
            " ", "ikj_a", "jik_a", "ikj_b_1_4",
            "jik_b_1_3", "jik_b_1_6", "jik_b_1_4",
            "jik_b_1_5", "jik_b_2_3",
            "jik_b_4_1", "jik_b_1_8" };

    switch (TYPE) {
        case 1:  //part a
            matmult_ikj_a(A, B, C, N);
            break;
        case 2:  //part a
            matmult_jik_a(A, B, C, N);
            break;

        case 3:  //part b
            matmult_ikj_b_1_4(A, B, C, N);
            break;

        case 4:  //part b
            matmult_jik_b_1_3(A, B, C, N);
            break;

        case 5:  //part b
            matmult_jik_b_1_6(A, B, C, N);
            break;

        case 6:  //part b
            matmult_jik_b_1_4(A, B, C, N);
            break;

        case 7:  //part b
            matmult_jik_b_1_5(A, B, C, N);
            break;

        case 8:  //part b
            matmult_jik_b_2_3(A, B, C, N);
            break;

        case 9:  //part b
            matmult_jik_b_4_1(A, B, C, N);
            break;

        case 10:  //part b
            matmult_jik_b_1_8(A, B, C, N);
            break;
    }
    std::cout << TYPE_STRING[ TYPE ];
}
