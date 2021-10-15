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
    _betaSE = gsl_vector_alloc(X->size2); // k * 1

    gsl_matrix *vCov = gsl_matrix_calloc(X->size2, X->size2); // k * k (Needs to be zeroed to prepare for triangular matrix)

    // Get the inverse of X'X and store it in vCov (later, we'll scale vCov to e'e/n-k to get the true vCov matrix)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0.0, vCov); // _vCov = (X' X) (lower dig.)
    // If X is invertible (i.e. X is full rank), then the matrix X'X is positive-definite (Can use cholesky)
    gsl_linalg_cholesky_decomp1(vCov);
    gsl_linalg_cholesky_invert(vCov); // XtX = (X'X)^-1

    // Estimate beta
    gsl_vector *Xty = gsl_vector_alloc(X->size2);
    gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, Xty); // Xty = (X'y)
    gsl_blas_dgemv(CblasNoTrans, 1.0, vCov, Xty, 0.0, _betaHat); // _betaHat = (X'X)^-1 * (X'y)
    gsl_vector_free(Xty);

    // Find SSR = sum(e**2) = ||e||^2 // VALIDATED
    gsl_vector *e = gsl_vector_alloc(y->size); // n * 1
    gsl_vector_memcpy(e, y); // e now holds the actual Y values
    gsl_blas_dgemv(CblasNoTrans, -1.0, X, _betaHat, 1.0, e); // e = y_act - y_est = -(X*B) + y_act
    _ssr = pow(gsl_blas_dnrm2(e), 2); // = sqrt(sum(e^2))^2 = sum(e^2)
    gsl_vector_free(e);

    // SE of Coef. =  sqrt(var_of_resids. * diag((X'X)**-1) // VALIDATED.
    double resvar = (1 / (double)(X->size1 - X->size2)) * _ssr;
    for (int i = 0; i < vCov->size1; i++) {
        double se = sqrt(gsl_matrix_get(vCov, i, i) * resvar);
        gsl_vector_set(_betaSE, i, se);
    }

    // Find RSQ = 1 - (SSR/SST) // VALIDATED
    _rsq = 1 - (_ssr/gsl_stats_tss(y->data, y->stride, y->size));
}

void LmOLS::getBetaHat(std::vector<double> &v) {
    // Pushes the betahat vector onto std::vector v
    for (int i = 0; i < _betaHat->size; i++) {
        v.push_back(gsl_vector_get(_betaHat, i));
    }
}

void LmOLS::getBetaHat(std::string &s) {
    // Pushes the betahat vector to a comma delim. string.
    for (int i = 0; i < _betaHat->size; i++) {
        s.append(std::to_string(gsl_vector_get(_betaHat, i)));
        s.append(",");
    }
}

void LmOLS::getBetaSE(std::vector<double> &v) {
    // Pushes the betahat vector onto std::vector v
    for (int i = 0; i < _betaSE->size; i++) {
        v.push_back(gsl_vector_get(_betaSE, i));
    }
}

void LmOLS::getBetaSE(std::string &s) {
    // Pushes the betahat vector to a comma delim. string.
    for (int i = 0; i < _betaSE->size; i++) {
        s.append(std::to_string(gsl_vector_get(_betaSE, i)));
        s.append(",");
    }
}

double LmOLS::getRSQ() {
    return _rsq;
}