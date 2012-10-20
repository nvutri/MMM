#include <algorithm>
void matmult_jik_c_1_4(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB);

void matmult_jik_c(double* A, double* B, double* C, unsigned N, unsigned NB) {
    for (unsigned k = 0; k < N; k += NB) {
        int KB = std::min(N-k, NB);
        for (unsigned i = 0; i < N; i += NB) {
            int IB = std::min(N-i, NB);
            for (unsigned j = 0; j < N; j += NB) {
                int JB = std::min(N-j, NB);
                /*sub matrix a, b, c*/
                double* a = &A(i, k);
                double* b = &B(k, j);
                double* c = &C(i, j);
                matmult_jik_c_1_4(a, b, c, N, JB, KB, IB);
            }
        }
    }

}
void matmult_jik_c_1_4(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB) {
    for (unsigned j = 0; j < JB; ++j) {
        for (unsigned i = 0; i < IB; i += 4) {
            // micro kernel
            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < KB; ++k) {
                register double b0;
                b0 = B(k, j);
                c0 += A(i, k) * b0;
                c1 += A(i + 1, k) * b0;
                c2 += A(i + 2, k) * b0;
                c3 += A(i + 3, k) * b0;
            }
            C(i, j) += c0;
            C(i + 1, j) += c1;
            C(i + 2, j) += c2;
            C(i + 3, j) += c3;
        }
    }
}
