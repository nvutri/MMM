/**
 * Mini Kernel MMM for part B
 * Unrolling loop j, i
 * Regiser Blocking
 */
// NU = 4
void matmult_ikj_b_1_4(double* A, double* B, double* C, unsigned N) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned k = 0; k < N; ++k) {
            register double a0 = A(i, k);
            for (unsigned j = 0; j < N; j += 4) {
                C(i, j)     += a0 * B(k, j);
                C(i, j + 1) += a0 * B(k, j + 1);
                C(i, j + 2) += a0 * B(k, j + 2);
                C(i, j + 3) += a0 * B(k, j + 3);
            }
        }
    }
}
/**
 * part b NU = 1 MU = 4
 */
void matmult_jik_b_1_4(double* A, double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 4) {
            // micro kernel

            register double c0, c1, c2, c3;
            c0 = c1 = c2 = c3 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                double* a = &A(i, k);
                b0 = B(k, j);
                c0 += a[0] * b0;
                c1 += a[1] * b0;
                c2 += a[2] * b0;
                c3 += a[3] * b0;
            }
            double* c = &C(i, j);
            c[0] += c0;
            c[1] += c1;
            c[2] += c2;
            c[3] += c3;
        }
    }
}


/**
 * NU = 4; MU = 1;
 */
void matmult_jik_b_4_1(double* A, double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; j += 4) {
        for (unsigned i = 0; i < N; ++i) {
            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double a0;
                double* b = &B(k, j);
                a0 = A(i, k);
                c0 += a0 * B(k, j);
                c1 += a0 * B(k, j + 1);
                c2 += a0 * B(k, j + 2);
                c3 += a0 * B(k, j + 3);
            }
            C(i, j    ) += c0;
            C(i, j + 1) += c1;
            C(i, j + 2) += c2;
            C(i, j + 3) += c3;
        }
    }
}

/* NU = 4 MU = 1 KU=4 */
void matmult_jik_b_4_4(double* A, double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; j += 4) {
        for (unsigned i = 0; i < N; ++i) {
            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < N; k+=4) {
                double a0 = A(i, k);
                double* b0 = &B(k, j);
                double* b1 = &B(k, j + 1);
                double* b2 = &B(k, j + 2);
                double* b3 = &B(k, j + 3);

                c0 += a0 * b0[0];
                c1 += a0 * b1[0];
                c2 += a0 * b2[0];
                c3 += a0 * b3[0];

                double a1 = A(i, k + 1);

                c0 += a1 * b0[1];
                c1 += a1 * b1[1];
                c2 += a1 * b2[1];
                c3 += a1 * b3[1];

                double a2 = A(i, k + 2);

                c0 += a2 * b0[2];
                c1 += a2 * b1[2];
                c2 += a2 * b2[2];
                c3 += a2 * b3[2];

                double a3 = A(i, k + 3);

                c0 += a3 * b0[3];
                c1 += a3 * b1[3];
                c2 += a3 * b2[3];
                c3 += a3 * b3[3];

            }
            C(i, j    ) += c0;
            C(i, j + 1) += c1;
            C(i, j + 2) += c2;
            C(i, j + 3) += c3;
        }
    }
}

/**
 * NU = 8; MU = 1;
 */
void matmult_jik_b_8_1(double* A, double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; j += 8) {
        for (unsigned i = 0; i < N; ++i) {
            register double c0, c1, c2, c3, c4, c5, c6, c7;
            c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double a0;
                a0 = A(i, k);
                c0 += a0 * B(k, j);
                c1 += a0 * B(k, j + 1);
                c2 += a0 * B(k, j + 2);
                c3 += a0 * B(k, j + 3);
                c4 += a0 * B(k, j + 4);
                c5 += a0 * B(k, j + 5);
                c6 += a0 * B(k, j + 6);
                c7 += a0 * B(k, j + 7);
            }
            C(i, j    ) += c0;
            C(i, j + 1) += c1;
            C(i, j + 2) += c2;
            C(i, j + 3) += c3;
            C(i, j + 4) += c4;
            C(i, j + 5) += c5;
            C(i, j + 6) += c6;
            C(i, j + 7) += c7;
        }
    }
}
/**
 * part b NU = 1 MU = 1 KU=8
 */
void matmult_jik_b_1_8(double* A, double* B, double* C, unsigned N) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 8) {

            register double c0, c1, c2, c3, c4, c5, c6 , c7;
            c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                double* a = &A(i, k);
                b0 = B(k, j);
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
}
