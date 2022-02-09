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

    // Fit model with built-in least-squares method
    gsl_multifit_linear(X, y, _betaHat, vCov, &_rsq, wrk);
    _rsq = 1 - (_rsq / gsl_stats_tss(y->data, y->stride, y->size));

    for (int i = 0; i < vCov->size1; i++) {
        double se = sqrt(gsl_matrix_get(vCov, i, i));
        gsl_vector_set(_betaSE, i, se);
    }

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

float LmOLS::getRSQ() const {
    return (float) _rsq;
}