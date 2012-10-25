// define sub-block A_
#define _A(x, y) _A[ (y) * N + (x)]
//#define _B(x, y) _B[ (y) * JB + (x)]

void matmult_jik_d_1_8(double* A, double* B, double* C, unsigned N, unsigned JB,
                       unsigned KB, unsigned IB);

void matmult_jik_d(double* A, double* B, double* C, unsigned N, unsigned NB) {

    // Local storage
    double* _A = alloc(N);
    // Copy Full A
    memcpy(_A, A, N * N * sizeof(double));

    for (unsigned j = 0; j < N; j += NB) {
        int JB = std::min(N - j, NB);

        for (unsigned i = 0; i < N; i += NB) {
            int IB = std::min(N - i, NB);
            for (unsigned k = 0; k < N; k += NB) {
                int KB = std::min(N - k, NB);
                /*sub matrix a, b, c*/
                double* a = &_A(i, k);
                double* b = &B(k, j);
                double* c = &C(i, j);
                matmult_jik_d_1_8(a, b, c, N, JB, KB, IB);
            }
        }
    }
    //Deallocate A and B
    free(_A);
//    free(_B);
}

/**
 * NU = 1; MU = 8;
 * mini-kernel MMM
 */
void matmult_jik_d_1_8(double* _A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB) {
    double* _B = (double*) malloc(sizeof(double) * KB);
    for (unsigned j = 0; j < JB; ++j) {
        // Copy a panel of B
        for (unsigned k = 0; k < KB; ++k) {
            _B[k] = B(k, j);
        }
        for (unsigned i = 0; i < IB; i += 8) {

            register double c0, c1, c2, c3, c4, c5, c6, c7;
            c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;

            for (unsigned k = 0; k < KB; ++k) {
                register double b0;
                double* a = &_A(i, k);

                b0 = _B[k];

                c0 += a[0] * b0;
                c1 += a[1] * b0;
                c2 += a[2] * b0;
                c3 += a[3] * b0;
                c4 += a[4] * b0;
                c5 += a[5] * b0;
                c6 += a[6] * b0;
                c7 += a[7] * b0;
            }
            double* c = &C(i, j);
            c[0] += c0;
            c[1] += c1;
            c[2] += c2;
            c[3] += c3;
            c[4] += c4;
            c[5] += c5;
            c[6] += c6;
            c[7] += c7;
        }
    }
    free(_B);
}
#undef _A
#undef _B

