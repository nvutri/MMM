// define sub-block A_
#define _A(x, y) _A[ (y) * KB + (x)]

void matmult_jik_d_1_8(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB);

void matmult_jik_d_4_4(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB);

void matmult_jik_d(double* A, double* B, double* C, unsigned N, unsigned NB) {
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
 * NU = 1; MU = 8;
 * mini-kernel MMM
 */
void matmult_jik_d_1_8(double* A, double* B, double* C, unsigned N,
                       unsigned JB, unsigned KB, unsigned IB)  {
    // Local storage
    double* _A = alloc(IB);
    double* _B = (double*) malloc(sizeof(double) * KB);

    // Copy Full A
    for (unsigned k = 0; k< KB; ++k)
        for (unsigned i = 0; i< IB; ++i)
            _A(i, k) = A(i, k);

    for (unsigned j = 0; j < JB; ++j) {
        // Copy a panel of B
        for (unsigned k = 0; k < KB; ++k){
            _B[k] = B(k, j);
        }
        // Do computation, no copy a tile for C
        for (unsigned i = 0; i < IB; i += 8) {

            register double c0, c1, c2, c3, c4, c5, c6 , c7;
            c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;

            for (unsigned k = 0; k < KB; ++k) {
                register double b0;
                double* a = &_A(i, k);

                b0  = _B[k];

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

    //Deallocate A and B
    free(_A);
    free(_B);
}
//void matmult_jik_d_4_4(double* A, double* B, double* C, unsigned N,
//                       unsigned JB, unsigned KB, unsigned IB)  {
//    for (unsigned j = 0; j < JB; j += 4) {
//        for (unsigned i = 0; i < IB; ++i) {
//            register double c0, c1, c2, c3;
//            c0 = 0.0;
//            c1 = 0.0;
//            c2 = 0.0;
//            c3 = 0.0;
//            double* b0 = &B(0, j);
//            double* b1 = &B(0, j + 1);
//            double* b2 = &B(0, j + 2);
//            double* b3 = &B(0, j + 3);
//
//            for (unsigned k = 0; k < KB; k+=4) {
//                double a0 = A(i, k);
//
//                c0 += a0 * b0[0];
//                c1 += a0 * b1[0];
//                c2 += a0 * b2[0];
//                c3 += a0 * b3[0];
//
//                double a1 = A(i, k + 1);
//
//                c0 += a1 * b0[1];
//                c1 += a1 * b1[1];
//                c2 += a1 * b2[1];
//                c3 += a1 * b3[1];
//
//                double a2 = A(i, k + 2);
//
//                c0 += a2 * b0[2];
//                c1 += a2 * b1[2];
//                c2 += a2 * b2[2];
//                c3 += a2 * b3[2];
//
//                double a3 = A(i, k + 3);
//
//                c0 += a3 * b0[3];
//                c1 += a3 * b1[3];
//                c2 += a3 * b2[3];
//                c3 += a3 * b3[3];
//
//                ++b0;
//                ++b1;
//                ++b2;
//                ++b3;
//
//            }
//            C(i, j    ) += c0;
//            C(i, j + 1) += c1;
//            C(i, j + 2) += c2;
//            C(i, j + 3) += c3;
//        }
//    }
//}
//
