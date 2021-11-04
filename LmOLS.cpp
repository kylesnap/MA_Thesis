// LmOLS.cpp
// Fits the log likelihood
// Kyle Dewsnap
// 16SEP21

#include "LmOLS.h"

LmOLS::LmOLS(const gsl_matrix *X, const gsl_vector *y) {
    // Fits the model with OLS

    _betaHat = gsl_vector_alloc(X->size2); // k * 1
    _betaSE = gsl_vector_alloc(X->size2); // k * 1

    gsl_multifit_linear_workspace *wrk = gsl_multifit_linear_alloc(X->size1, X->size2);
    gsl_matrix *vCov = gsl_matrix_alloc(X->size2, X->size2); // k * k

    // Fit model with stock least-squares method
    gsl_multifit_linear(X, y, _betaHat, vCov, &_rsq, wrk);
    _rsq = 1 - (_rsq / gsl_stats_tss(y->data, y->stride, y->size));

    for (int i = 0; i < vCov->size1; i++) {
        double se = sqrt(gsl_matrix_get(vCov, i, i));
        gsl_vector_set(_betaSE, i, se);
    }


        /*
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
        _rsq = 1 - (_ssr/ */
}

void LmOLS::getBetaHat(std::vector<float> &v) {
    // Pushes the betahat vector onto std::vector v
    for (int i = 0; i < _betaHat->size; i++) {
        v.push_back((float) gsl_vector_get(_betaHat, i));
    }
}

void LmOLS::getBetaSE(std::vector<float> &v) {
    // Pushes the betahat vector onto std::vector v
    for (int i = 0; i < _betaSE->size; i++) {
        v.push_back((float) gsl_vector_get(_betaSE, i));
    }
}

float LmOLS::getRSQ() {
    return (float) _rsq;
}