void matmult_jik_d_4_1(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB);

void matmult_jik_d(double* A, double* B, double* C, unsigned N, unsigned NB) {
    for (unsigned k = 0; k < N; k += NB) {
        int KB = std::min(N-k, NB);
        for (unsigned i = 0; i < N; i += NB) {
            int IB = std::min(N-i, NB);
            // Copy block A(i,k) size NB here

            for (unsigned j = 0; j < N; j += NB) {
                int JB = std::min(N-j, NB);
                // Copy row B(k,j) NB here
                // Copy row C(k,j) NB here
                /*sub matrix a, b, c*/
                double* a = &A(i, k);
                double* b = &B(k, j);
                double* c = &C(i, j);
                matmult_jik_d_4_1(a, b, c, N, JB, KB, IB);
            }
        }
    }
}

/**
 * NU = 4; MU = 1;
 * mini-kernel MMM
 */
void matmult_jik_d_4_1(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB)  {
    for (unsigned j = 0; j < JB; j += 4) {
        for (unsigned i = 0; i < IB; ++i) {

            // micro-kernel MMM
            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < KB; ++k) {
                register double a0;
                a0 = A(i, k);
                c0 += a0 * B(k, j);
                c1 += a0 * B(k, j + 1);
                c2 += a0 * B(k, j + 2);
                c3 += a0 * B(k, j + 3);
            }
            C(i, j) += c0;
            C(i, j + 1) += c1;
            C(i, j + 2) += c2;
            C(i, j + 3) += c3;
        }
    }
}
