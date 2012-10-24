#include <algorithm>
#include <emmintrin.h>

void printMatrix(const double* X, unsigned N);
void matmult_jik_2_1(const double* A, const double* B, double* C, unsigned N,
                     unsigned JB, unsigned KB, unsigned IB);
void matmult_jik_1_2(const double* A, const double* B, double* C, unsigned N,
                     unsigned JB, unsigned KB, unsigned IB);

//void matmult(double* A, double* B, double* C, unsigned N) {
//    unsigned NB = 64;
//    for (unsigned j = 0; j < N; j += NB) {
//        int JB = std::min(N - j, NB);
//        for (unsigned i = 0; i < N; i += NB) {
//            int IB = std::min(N - i, NB);
//            for (unsigned k = 0; k < N; k += NB) {
//                int KB = std::min(N - k, NB);
//
//                /*sub matrix a, b, c*/
//                const double* a = (double*) &A(i, k);
//                const double* b = (double*) &B(k, j);
//                double* c = (double*) &C(i, j);
//                matmult_jik_2_1(a, b, c, N, JB, KB, IB);
//            }
//        }
//    }
//}
void matmult_jik_1_2(const double* A, const double* B, double* C, unsigned N,
                     unsigned JB, unsigned KB, unsigned IB) {
    for (unsigned j = 0; j < JB; ++j) {
        for (unsigned i = 0; i < IB; i += 2) {
            // micro kernel
            __m128d rA, rB, rC, rT;
            rC = _mm_load_sd(&C(i, j));

            for (unsigned k = 0; k < KB; ++k) {
                rA = _mm_load_sd(&A(i, k));
                rB = _mm_set_sd(B(k, j));
                rT = _mm_mul_sd(rA, rB);
                rC = _mm_add_sd(rT, rC);

            }
            _mm_store_sd(&C(i, j), rC);
        }
    }
}

/**
 * NU = 4; MU = 1;
 */
void matmult_jik_2_1(const double* A, const double* B, double* C, unsigned N,
                     unsigned JB, unsigned KB, unsigned IB) {
    for (unsigned j = 0; j < JB; j += 2) {
        for (unsigned i = 0; i < IB; ++i) {
            __m128d rA, rB, rC, rT;
            rC = _mm_load_sd(&C(i, j));
            for (unsigned k = 0; k < KB; ++k) {
                rA = _mm_load_sd(&A(i, k));
                rB = _mm_load_sd(&B(k, j));
                rT = _mm_mul_sd(rA, rB);
                rC = _mm_add_sd(rT, rC);
            }
            _mm_store_sd(&C(i, j), rC);
//            printMatrix(C, N);
        }
    }
}

void matmult(const double* A, const double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; j += 1) {
        for (unsigned i = 0; i < N; ++i) {
            __m128d rA, rB, rC, rT;
            rC = _mm_load_sd(&C(i, j));
            for (unsigned k = 0; k < N; ++k) {
                rA = _mm_load_sd(&A(i, k));
                rB = _mm_load_sd(&B(k, j));
                rT = _mm_mul_sd(rA, rB);
                rC = _mm_add_sd(rT, rC);
            }
            _mm_store_sd(&C(i, j), rC);
        }
    }
}

//void matmult_jik_1_4(const double* A, const double* B, double* C, unsigned N,
//                       unsigned JB, unsigned KB, unsigned IB)  {
//    for (unsigned j = 0; j < JB; ++j) {
//        for (unsigned i = 0; i < IB; i += 4) {
//            // micro kernel
////            __m128d rA, rB, rC;
//            register double c0, c1, c2, c3;
//            c0 = c1 = c2 = c3 = 0.0;
////            rA = __m_load_sd
//            for (unsigned k = 0; k < KB; ++k) {
//                register double b0;
//                const double* a = &A(i, k);
//                b0 = B(k, j);
//                c0 += a[0] * b0;
//                c1 += a[1] * b0;
//                c2 += a[2] * b0;
//                c3 += a[3] * b0;
//            }
//            double* c = &C(i, j);
//            c[0] += c0;
//            c[1] += c1;
//            c[2] += c2;
//            c[3] += c3;
//        }
//    }
//}

