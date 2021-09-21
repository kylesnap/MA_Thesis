//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>

#include "LogLik.h"

dlib::matrix<double> LogLik::logLik(const colvec &theta) {
    // The log likelihood of theta given the data for the logistic regression model

    // For one data point:
    // Pr(Y = y|X = x) = f(theta, x)**y * [1 - f(theta, x)]**(1 - y)
    // f(theta; x) = 1 / (1 + e**(-g(theta; x))
    // g(theta; x) = sum(i = 0 -> k) { x[i]theta[i] } = x * theta' (where x[1, k] and theta[1, k])

    //Find g(theta; x)
    //
    // TO DO::: and f(theta; x)
    dlib::matrix<double> g = _x * dlib::trans(theta);
    return g;

    // l(theta) = sum(i = 0 -> N) { y[i]ln(f(theta;x[i])) * [1 - f(theta; x[i])]**(1 - y[i])
    // f(theta;x) = 1 / (1 + e**-g(theta; x)); g(theta; x) = theta'x

    // Finding f(theta; x) for all i in x would make things speedier...
}

colvec LogLik::logLik_d(const colvec &theta) {
    return colvec();
}
