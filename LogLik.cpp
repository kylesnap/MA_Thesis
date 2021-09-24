//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>
#include "LogLik.h"

LogLik::LogLik(const dlib::matrix<double>& x, const dlib::matrix<double>& y) {
    // Resize matrices after checking they're legal (n of x == n of y)
    if (x.nc() != y.nc()) throw std::bad_function_call();

    _x = x; // k * n
    _y = y; // 1 * n
}

double p(const colvec &x_i, const colvec &theta) {
    // The logistic function for one observation
    // p(x; theta) = 1 / (1 + -exp(theta[0] + x * theta[1:k])
    double t0 = theta(0);
    double lc = dlib::dot(x_i, dlib::subm(theta, dlib::range(1, theta.nr()), dlib::range(0, 0)));
    return 1 / (1 + exp(-1 * (t0 + lc)));
}

double LogLik::logLik(const colvec &theta) {
    // Theta length has to be equal to the number of variables + 1 (the intercept)
    if (theta.nr() != (_x.nr() + 1)) { throw std::bad_function_call(); }

    // LogLik = sum[i -> n] { y[i] ln(p(x[i], theta)) + (1 - y)ln(1-p(x[i], theta))
    double ll = 0;
    for (int i = 0; i < _x.nc(); i++) {
        colvec x_i = dlib::colm(_x, i);
        double p_i = p(x_i, theta);
        ll += (_y(i) * log(p_i)) + ((1 - _y(i)) * (log(1 - p_i)));
    }

    return ll;
}

colvec LogLik::logLik_d(const colvec &theta) {
    // Theta length has to be equal to the number of variables + 1 (the intercept)
    if (theta.nr() != (_x.nr() + 1)) { throw std::bad_function_call(); }
    colvec r(theta.nr(), 1);
    std::vector<double> dt;

    // d(l) / d(theta_j) = sum[i -> n] { (y_i - p(x[i]))x[i,j]
    for (int j = 0; j < theta.nr(); j++) {
        double d = 0;
        for (int i = 0; i < _x.nc(); i++) {
            colvec x_i = dlib::colm(_x, i);
            double p_i = p(x_i, theta);
            if (j == 0) {
                d += (_y(i) - p_i);
            } else {
                d += (_y(i) - p_i) * x_i(j - 1);
            }
        }
        dt.push_back(d);
    }

    r = dlib::mat(dt);
    return(r);
}