// LmOLS.cpp
// Fits the log likelihood
// Kyle Dewsnap
// 16SEP21

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>

#include "LmOLS.h"

LmOLS::LmOLS(const gsl_matrix *X, const gsl_vector *y) {
    // Fits the model with OLS

    _betaHat = gsl_vector_alloc(X->size2); // k * 1
    _vCov = gsl_matrix_calloc(X->size2, X->size2); // k * k (Needs to be zeroed to prepare for triangular matrix)

    // Get the inverse of X'X and store it in _vCov (later, we'll scale _vCov to e'e/n-k to get the true vCov matrix)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0.0, _vCov); // _vCov = (X' X) (lower dig.)
    // If X is invertible (i.e. X is full rank), then the matrix X'X is positive-definite (Can use cholesky)
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

    // Find RSQ = 1 - (SSR/SST)
    double sst = gsl_stats_tss(y->data, y->stride, y->size);
    _rsq = 1 - (_ssr/sst);

    // VcoV matrix calculation _vCov = (SSR / (n - k))^2 * (X'X)^-1
    gsl_matrix_scale(_vCov, pow((_ssr / (double)(X->size1 - X->size2)), 2));
    //TODO: Test constructor with "handwritten" case.
}

void LmOLS::getBetaHat(std::vector<double> &v) {
    // Pushes the betahat vector onto std::vector v
    for (int i = 0; i < _betaHat->size; i++) {
        v.push_back(gsl_vector_get(_betaHat, i));
    }
}

void LmOLS::getBetaHat(std::string &s) {
    // Appends a comma separated string of beta hats (0, 1, 2) to string s
    if (_betaHat->size != 3) {
        std::cerr << "String method improper length: Beta";
        return;
    }
    char b[64];
    sprintf(b, "%f,%f,%f",
            gsl_vector_get(_betaHat, 0),
            gsl_vector_get(_betaHat, 1),
            gsl_vector_get(_betaHat, 2)
            );
    s.append(b);
}

void LmOLS::getBetaSE(std::vector<double> &v) {
    // Pushes the std errs of the beta coef (i.e. the squared diag. of vCov) to v.
    for (int i = 0; i < _vCov->size1; i++) {
        v.push_back(sqrt(gsl_matrix_get(_vCov, i, i)));
    }
}

void LmOLS::getBetaSE(std::string &s) {
    // Appends a comma separated string of beta ses (0, 1, 2) to string s
    if (_vCov->size1 != 3) {
        std::cerr << "String method improper length: Beta SE";
        return;
    }
    char b[64];
    sprintf(b, "%f,%f,%f",
            sqrt(gsl_matrix_get(_vCov, 0, 0)),
            sqrt(gsl_matrix_get(_vCov, 1, 1)),
            sqrt(gsl_matrix_get(_vCov, 2, 2))
    );
    s.append(b);
}

double LmOLS::getRSQ() {
    return _rsq;
}