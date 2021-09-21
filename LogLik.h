// LogLik.h
// Defines col_vec and LogLik
// Kyle Dewsnap
// 20SEP21

#include <dlib/matrix.h>

#include <utility>

#pragma once

typedef dlib::matrix<double, 1, 2> colvec;

class LogLik {
private:
    dlib::matrix<double> _x; // n * k matrix (first col are all ones for intercept term)
    dlib::matrix<double> _y; // n * 1 binary [0 or 1] responses
public:
    LogLik(dlib::matrix<double> x, dlib::matrix<double> y) : _x(std::move(x)), _y(std::move(y)) {};

    double logLik(const colvec &theta);
    colvec logLik_d(const colvec &theta);
};
