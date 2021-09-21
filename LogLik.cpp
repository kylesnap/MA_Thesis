//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>
#include "LogLik.h"

double logf(const double &z) {
    // Returns the logistic function at z; 1 / (1 + e**(-z))
    return 1 / (1 + exp(-z));
}

double LogLik::logLik(const colvec &theta) {
    // The log likelihood of theta given the data for the logistic regression model

    dlib::matrix<double> g(_x.nr(), 1);
    double loglik = 0;

    // FOR ONE DATA POINT
    // Pr(Y = y|X = x) = f(theta, x)**y * [1 - f(theta, x)]**(1 - y)
    // f(theta; x) = 1 / (1 + e**(-g(theta; x))
    // g(theta; x) = sum(i = 0 -> k) { x[i]theta[i] } = x * theta' (where x[1, k] and theta[1, k])

    // l(theta) = (i = 0 -> N) { y[i]ln(f(theta;x[i])) * [1 - f(theta; x[i])]*ln(1 - f(theta; x[i]))
    g = _x * dlib::trans(theta); // Linear combination of x and theta
    for (int i = 0; i <= g.nc(); i++) {
        double temp = logf(g(i)); // Applying logfunc to i-th element of the combination of x and theta.
        int y = int(_y(i)); // Getting i-th y var.
        loglik += (y * log(temp)) + ((1 - y) * log(1 - temp));
    }

    return loglik;
}

colvec LogLik::logLik_d(const colvec &theta) {
    return colvec();
}
