/**
 * Simple Matrix Multiplication.
 * For Part A
 */
void matmult_ikj_a(double* A, double* B, double* C) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned k = 0; k < N; ++k) {
            for (unsigned j = 0; j < N; ++j) {
                C(i,j) += A(i,k) * B(k,j);
            }
        }
    }
}

void matmult_jik_a(double* A, double* B, double* C) {
    for (unsigned j = 0; j < N; ++j) {
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned k = 0; k < N; ++k) {
                C(i,j) += A(i,k) * B(k,j);
            }
        }
    }
}
