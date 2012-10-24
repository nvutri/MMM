//#include <algorithm>
void matmult_4_4(double* A, double* B, double* C, unsigned N, unsigned NB);
void matmult_1_1(double* A, double* B, double* C, unsigned N);

/**
 * Helps and Ideas from Professor Robert Van De Geijn.
 * NB = 16.
 * NU = KU = 4
 */
void matmult(double* A, double* B, double* C, unsigned N) {
    unsigned NB = 16; //TODO: Change back to 64
    int LB = N % NB;  // Left over tile
    int NN = N - LB;  // The nice size
    for (unsigned j = 0; j < NN; j += NB) {
        for (unsigned i = 0; i < NN; i += NB) {
            for (unsigned k = 0; k < NN; k += NB) {
                /*sub matrix a, b, c*/
                double* a = &A(i, k);
                double* b = &B(k, j);
                double* c = &C(i, j);
                matmult_4_4(a, b, c, N, NB);
            }
        }
    }
    // Resolving the leftover tile
    if (LB > 0) {
        // NBxNB of C
        for (unsigned j = 0; j < NN; ++j) {
            for (unsigned i = 0; i < NN; ++i) {
                for (unsigned k = NN; k < N; ++k) {
                    C(i,j) += A(i,k) * B(k,j);
                }
            }
        }
        // Bottom panel of C
        for (unsigned j = 0; j < N; ++j) {
            for (unsigned i = NN; i < N; ++i) {
                for (unsigned k = 0; k < N; ++k) {
                    C(i,j) += A(i,k) * B(k,j);
                }
            }
        }
        // Right tile of C.
        // Watch out for overlapping LBxLB at bottom right corner
        for (unsigned j = NN; j < N; ++j) {
            for (unsigned i = 0; i < N - LB; ++i) {
                for (unsigned k = 0; k < N; ++k) {
                    C(i,j) += A(i,k) * B(k,j);
                }
            }
        }

    }
}

void matmult_4_4(double* A, double* B, double* C, unsigned N, unsigned NB) {
    for (unsigned j = 0; j < NB; j += 4) {
        for (unsigned i = 0; i < NB; ++i) {
            register double c0, c1, c2, c3, a0;
            c0 = c1 = c2 = c3 = 0.0;
            double* a = &A(0, i);
            double* b0 = &B(0, j);
            double* b1 = &B(0, j + 1);
            double* b2 = &B(0, j + 2);
            double* b3 = &B(0, j + 3);

            for (unsigned k = 0; k < NB; k += 4) {
                a0 = *a++;

                c0 += a0 * b0[0];
                c1 += a0 * b1[0];
                c2 += a0 * b2[0];
                c3 += a0 * b3[0];

                a0 = *a++;

                c0 += a0 * b0[1];
                c1 += a0 * b1[1];
                c2 += a0 * b2[1];
                c3 += a0 * b3[1];

                a0 = *a++;

                c0 += a0 * b0[2];
                c1 += a0 * b1[2];
                c2 += a0 * b2[2];
                c3 += a0 * b3[2];

                a0 = *a++;

                c0 += a0 * b0[3];
                c1 += a0 * b1[3];
                c2 += a0 * b2[3];
                c3 += a0 * b3[3];

                b0 += 4;
                b1 += 4;
                b2 += 4;
                b3 += 4;

            }
            C(i, j) += c0;
            C(i, j + 1) += c1;
            C(i, j + 2) += c2;
            C(i, j + 3) += c3;
        }
    }
}
