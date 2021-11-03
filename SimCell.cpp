// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include "SimCell.h"

SimCell::SimCell(int n, float rSq, std::tuple<float, float> betaP, std::tuple<float, float> betaQ,
                 std::tuple<float, float> propGrps, std::tuple<float, float, float, float> paramsTrue, gsl_rng *r) {

    // Error checks for RSQ
    if ((rSq - 1) * rSq >= 0) throw std::out_of_range("RSQ is a value outside of [0, 1]");
    _rSq = rSq;
    _r = r;

    // Handles odd sample sizes silently.
    n = n%2 == 1 ? n - 1 : n;

    // Construct objects -- Errors will be thrown by class methods.
    auto *bP = new BetaGen(std::get<0>(betaP), std::get<1>(betaP), r);
    auto *bQ = new BetaGen(std::get<0>(betaQ), std::get<1>(betaQ), r);
    _xMat = new DesignMat(n, bP, bQ, r, std::get<0>(propGrps), std::get<1>(propGrps));

    gsl_vector_set(_pTrue, 0, std::get<0>(paramsTrue));
    gsl_vector_set(_pTrue, 1, std::get<1>(paramsTrue));
    gsl_vector_set(_pTrue, 2, std::get<2>(paramsTrue));
    gsl_vector_set(_pTrue, 3, std::get<3>(paramsTrue));
}

void SimCell::toVec(std::vector<float> &v, bool print) {
    // Print SimCell params (All Params of Design Mat, true params, RSQ).
    _xMat->getDesignMat(v);
    for (int i = 0; i < 4; i++) {
        v.push_back(gsl_vector_get(_pTrue, i));
    }
    v.push_back(_rSq);
    if (print) {
        for (float i : v) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
}

void SimCell::run() {
    // Build a vector of 'true' Y responses.
    // RSQ = 1 - (var[Residuals] / var[Total])
    // var(e) = (1- RSQ) * sum(B[j]**2) / RSQ
    float err_var = pow(gsl_vector_get(_pTrue, 1), 2) + pow(gsl_vector_get(_pTrue, 2), 2) +
            pow(gsl_vector_get(_pTrue, 3), 2);
    err_var = (err_var * (1 - _rSq)) / _rSq;
    if (err_var <= 0) throw std::bad_function_call();

    // Add error terms to tY.
    gsl_vector *tY = gsl_vector_calloc(_xMat->getTX()->size1);
    for (int i = 0; i < tY->size; i++) {
        gsl_vector_set(tY, i, gsl_ran_gaussian(_r, err_var));
    }

    // Build tY = X * Params + e
    gsl_blas_dgemv(CblasNoTrans, 1.0, _xMat->getTX(), _pTrue, 1.0, tY);

    // TODO: Check RSQ of true data
    // TODO: Check model properties of response data and build logger
    // TODO: Build a monte carlo loop
    // TODO: Final testing!
}