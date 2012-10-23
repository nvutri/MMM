#include <algorithm>
#include <xmmintrin.h>


void matmult_jik_4_1(const float* A, const float* B, float* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB);

void matmult(double* A, double* B, double* C, unsigned N) {
    unsigned NB = 44;
    for (unsigned j = 0; j < N; j += NB) {
        int JB = std::min(N-j, NB);
        for (unsigned i = 0; i < N; i += NB) {
            int IB = std::min(N-i, NB);
            for (unsigned k = 0; k < N; k += NB) {
                int KB = std::min(N-k, NB);

                /*sub matrix a, b, c*/
                double* a = &A(i, k);
                double* b = &B(k, j);
                double* c = &C(i, j);
                matmult_jik_d_1_8(a, b, c, N, JB, KB, IB);
            }
        }
    }
}

/**
 * NU = 4; MU = 1;
 */
void matmult_jik_4_1(const float* A, const float* B, float* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB)  {
    for (unsigned j = 0; j < JB; j += 8) {
        for (unsigned i = 0; i < IB; ++i) {
//            register double c0, c1, c2, c3;
//            c0 = 0.0;
//            c1 = 0.0;
//            c2 = 0.0;
//            c3 = 0.0;
//            double* c = &C(i, j);
            __m128 rA, rB, rC;

            for (unsigned k = 0; k < KB; ++k) {
                register double a0;
                a0 = A(i, k);
//                double* b = &B(k, j);
                rA = _mm_load_ps(&A(i,k));
                rB = _mm_load_ps(&B(k,j));
                rC = _mm_mul_ps(rA, rB);
                _mm_store_ps(&C(i,j), rC);
            }
//            c[0] += c0;
//            c[1] += c1;
//            c[2] += c2;
//            c[3] += c3;
        }
    }
}

