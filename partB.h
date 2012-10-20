/**
 * Mini Kernel MMM for part B
 * Unrolling loop j, i
 * Regiser Blocking
 */
void matmult_ikj_b_1_4(double* A, double* B, double* C) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned k = 0; k < N; ++k) {
            register double a0 = A(i, k);
            for (unsigned j = 0; j < N; j += 4) {
                double* c = &C(i, j);
                double* b = &B(k, j);
                c[0] += a0 * b[0];
                c[1] += a0 * b[1];
                c[2] += a0 * b[2];
                c[3] += a0 * b[3];
            }
        }
    }
}
/**
 * part b 1_6
 */
void matmult_jik_b_1_6(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 6) {

            register double c0, c1, c2, c3, c4, c5;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            c5 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                b0 = B(k, j);
                c0 += A(i, k) * b0;
                c1 += A(i + 1, k) * b0;
                c2 += A(i + 2, k) * b0;
                c3 += A(i + 3, k) * b0;
                c4 += A(i + 4, k) * b0;
                c5 += A(i + 5, k) * b0;
            }
            C(i, j) += c0;
            C(i + 1, j) += c1;
            C(i + 2, j) += c2;
            C(i + 3, j) += c3;
            C(i + 4, j) += c4;
            C(i + 5, j) += c5;
        }
    }
}
/**
 * part b 1 _ 4
 */
void matmult_jik_b_1_4(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 4) {

            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
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

/**
 *  NU = 1; MU = 5;
 */
void matmult_jik_b_1_5(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 5) {

            register double c0, c1, c2, c3, c4;
            c0 = c1 = c2 = c3 = c4 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                b0 = B(k, j);
                c0 += A(i, k) * b0;
                c1 += A(i + 1, k) * b0;
                c2 += A(i + 2, k) * b0;
                c3 += A(i + 3, k) * b0;
                c4 += A(i + 4, k) * b0;
            }
            C(i, j) += c0;
            C(i + 1, j) += c1;
            C(i + 2, j) += c2;
            C(i + 3, j) += c3;
            C(i + 4, j) += c4;
        }
    }
}

/**
 * NU = 2; MU = 3;
 */
void matmult_jik_b_2_3(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; j += 2) {
        for (unsigned i = 0; i < N; i += 3) {

            register double c00, c01, c02, c10, c11, c12;
            c00 = c01 = c02 = 0;
            c10 = c11 = c12 = 0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0, b1;
                b0 = B(k, j);
                b1 = B(k, j + 1);
                c00 += A(i, k) * b0;
                c01 += A(i + 1, k) * b0;
                c02 += A(i + 2, k) * b0;

                c10 += A(i, k) * b1;
                c11 += A(i + 1, k) * b1;
                c12 += A(i + 2, k) * b1;
            }
            C(i, j) += c00;
            C(i + 1, j) += c01;
            C(i + 2, j) += c02;

            C(i, j + 1) += c10;
            C(i + 1, j + 1) += c11;
            C(i + 2, j + 1) += c12;
        }
    }
}

/**
 * part b 1_3
 */
void matmult_jik_b_1_3(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; i += 3) {

            register double c0, c1, c2;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double b0;
                b0 = B(k, j);
                c0 += A(i, k) * b0;
                c1 += A(i + 1, k) * b0;
                c2 += A(i + 2, k) * b0;
            }
            C(i, j) += c0;
            C(i + 1, j) += c1;
            C(i + 2, j) += c2;
        }
    }
}

/**
 * NU = 4; MU = 1;
 */
void matmult_jik_b_4_1(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; j += 4) {
        for (unsigned i = 0; i < N; ++i) {
            register double c0, c1, c2, c3;
            c0 = 0.0;
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            for (unsigned k = 0; k < N; ++k) {
                register double a0;
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
