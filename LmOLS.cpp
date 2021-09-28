// LmOLS.cpp
// Fits the log likelihood
// Kyle Dewsnap
// 16SEP21

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "LmOLS.h"

LmOLS::LmOLS(const gsl_matrix *X, const gsl_vector *y) {
    // Fits the model with OLS

    _betaHat = gsl_vector_alloc(X->size2); // k * 1
    _vCov = gsl_matrix_calloc(X->size2, X->size2); // k * k (Needs to be zeroed to prepare for triangular matrix)

    // Get the inverse of X'X and store it in _vCov (later, we'll scale _vCov to e'e/n-k to get the true vCov matrix)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0.0, _vCov); // _vCov = (X' X) (Only lower dig)
    // If X is invertible (e.g. X is full rank), then the matrix X'X is positive-definite (Can use cholesky)
    gsl_linalg_cholesky_decomp1(_vCov);
    gsl_linalg_cholesky_invert(_vCov); // XtX = (X'X)^-1

    // Estimate beta
    gsl_vector *Xty = gsl_vector_alloc(X->size2);
    gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, Xty); // Xty = (X'y)
    gsl_blas_dgemv(CblasNoTrans, 1.0, _vCov, Xty, 0.0, _betaHat); // _betaHat = (X'X)^-1 * (X'y)
    gsl_vector_free(Xty);

    // Find SSR = sum(e**2) = ||e||^2
    gsl_vector *e = gsl_vector_alloc(y->size); // n * 1
    gsl_vector_memcpy(e, y); // e now holds the actual Y values
    gsl_blas_dgemv(CblasNoTrans, -1.0, X, _betaHat, 1.0, e); // e = y_act - y_est = -(X*B) + y_act
    _ssr = pow(gsl_blas_dnrm2(e), 2); // = sqrt(sum(e^2))^2 = sum(e^2)
    gsl_vector_free(e);

    //VcoV matrix calculation
    gsl_matrix_scale(_vCov, (_ssr / (double)(X->size1 - X->size2))); // _vCov = SSR / (n - k) * (X'X)^-1

    for (int k = 0; k < _betaHat->size; k++) {
        printf("[%f]\n",
               gsl_vector_get(_betaHat, k)
        );
    }

    for (int i = 0; i < _betaHat->size; i++) {
        printf("[%f, %f, %f]\n",
               gsl_matrix_get(_vCov, i, 0), gsl_matrix_get(_vCov, i, 1), gsl_matrix_get(_vCov, i, 2)
        );
    }
    //TODO: Test constructor with "handwritten" case.
}
