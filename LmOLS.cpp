// LmOLS.cpp
// Fits the log likelihood
// Kyle Dewsnap
// 16SEP21

#include "LmOLS.h"

LmOLS::LmOLS(const gsl_matrix *X, const gsl_vector *y) {
    // Fits the model with OLS

    _betaHat = gsl_vector_alloc(X->size2); // k * 1
    _betaSE = gsl_vector_alloc(X->size2); // k * 1
    _betaHatN = gsl_vector_alloc(X->size2); // k * 1
    _betaSEN = gsl_vector_alloc(X->size2); // k * 1

    gsl_matrix * XN = gsl_matrix_alloc(X->size1, X->size2);
    gsl_vector * yN = gsl_vector_alloc(X->size1);

    // Normalize XN
    gsl_matrix_set_all(XN, 1);
    for (int i = 1; i < XN->size2; i++) { // Start at the second column, because the first is for c. intercept.
        gsl_vector_const_view temp = gsl_matrix_const_column(X, i);
        double mean = gsl_stats_mean(temp.vector.data, temp.vector.stride, temp.vector.size);
        double sd = gsl_stats_sd(temp.vector.data, temp.vector.stride, temp.vector.size);
        gsl_vector_view norm = gsl_matrix_column(XN, i);
        gsl_vector_set_all(&norm.vector, mean);
        gsl_vector_axpby(1, &temp.vector, -1, &norm.vector); // x[i] = x[i] + (mean(x[i]) * -1)
        gsl_vector_scale(&norm.vector, (1 / sd)); // x[i] = x[i] * (1 / sd(x[i]))
    }

    // Normalize yN
    gsl_vector_set_all(yN, gsl_stats_mean(y->data, y->stride, y->size));
    gsl_vector_axpby(1, y, -1, yN); // y[i] = y[i] + (-1 * mean(y[i]))
    gsl_vector_scale(yN, (1 / gsl_stats_sd(y->data, y->stride, y->size)));

    // Fit unstandardized model
    gsl_multifit_linear_workspace * wrk = gsl_multifit_linear_alloc(X->size1, X->size2);
    gsl_matrix * vCov = gsl_matrix_alloc(X->size2, X->size2); // k * k

    gsl_multifit_linear(X, y, _betaHat, vCov, &_rsq, wrk);
    _rsq = 1 - (_rsq / gsl_stats_tss(y->data, y->stride, y->size));

    for (int i = 0; i < vCov->size1; i++) {
        double se = sqrt(gsl_matrix_get(vCov, i, i));
        gsl_vector_set(_betaSE, i, se);
    }

    // Fit standardized model
    gsl_matrix_set_zero(vCov);

    gsl_multifit_linear(XN, yN, _betaHatN, vCov, &_rsqN, wrk);
    _rsqN = 1 - (_rsqN / gsl_stats_tss(yN->data, yN->stride, yN->size));

    for (int i = 0; i < vCov->size1; i++) {
        double se = sqrt(gsl_matrix_get(vCov, i, i));
        gsl_vector_set(_betaSEN, i, se);
    }

}

void LmOLS::getBetaHat(std::vector<float> &v, bool norm) {
    // Pushes the betahat vector onto std::vector v
    if (norm) {
        for (int i = 0; i < _betaHatN->size; i++) {
            v.push_back((float) gsl_vector_get(_betaHatN, i));
        }
        return;
    }
    for (int i = 0; i < _betaHat->size; i++) {
        v.push_back((float) gsl_vector_get(_betaHat, i));
    }
}

void LmOLS::getBetaSE(std::vector<float> &v, bool norm) {
    // Pushes the betahat vector onto std::vector v
    if (norm) {
        for (int i = 0; i < _betaSEN->size; i++) {
            v.push_back((float) gsl_vector_get(_betaSEN, i));
        }
        return;
    }
    for (int i = 0; i < _betaSE->size; i++) {
        v.push_back((float) gsl_vector_get(_betaSE, i));
    }
}

float LmOLS::getRSQ(bool norm) const {
    if (norm) return (float) _rsqN;
    return (float) _rsq;
}