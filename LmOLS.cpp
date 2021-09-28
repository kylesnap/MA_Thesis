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
    _vCov = gsl_matrix_alloc(X->size2, X->size2); // k * k

    // Need the inverse of X'X, so lets get that quickly...
    gsl_matrix *XtX = gsl_matrix_calloc(X->size2, X->size2); // k * k (Calloc because it's a triangular matrix)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0.0, XtX); // XtX = (X' X) (Only lower dig)
    // If X is invertible (e.g. X is full rank), then the matrix XtX = X'X is positive-definite (Can use cholesky)
    gsl_linalg_cholesky_decomp1(XtX);
    gsl_linalg_cholesky_invert(XtX); // XtX = (X'X)^-1

    // Estimate beta
    gsl_vector *Xty = gsl_vector_alloc(X->size2);
    gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, Xty); // Xty = (X'y)

    gsl_blas_dgemv(CblasNoTrans, 1.0, XtX, Xty, 0.0, _betaHat); // beta_h = (X'X)^-1 * (X'y)

    // Find yHat
    gsl_vector *yHat = gsl_vector_alloc(y->size);
    gsl_vector_memcpy(yHat, y); // yHat now holds the actual Y values
    gsl_blas_dgemv(CblasNoTrans, -1.0, X, _betaHat, 1.0, yHat); // = (-(X*Beta) + y_true)
    _ssr = pow(gsl_blas_dnrm2(yHat), 2); // = sqrt(sum(res^2))^2 = sum(res^2)

    //TODO: VcoV matrix calculation
    printf("%f\n", _ssr);

    for (int k = 0; k < _betaHat->size; k++) {
        printf("[%f]\n",
               gsl_vector_get(_betaHat, k)
        );
    }

    //TODO: Test constructor with "handwritten" case.
}
