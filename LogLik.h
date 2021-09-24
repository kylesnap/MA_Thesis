// LogLik.h
// Defines col_vec and LogLik
// Kyle Dewsnap
// 20SEP21

#include <dlib/matrix.h>

#include <utility>

#pragma once

typedef dlib::matrix<double, 0, 1> colvec;

double p(const colvec &x_i, const colvec &theta);

class LogLik {
private:
    dlib::matrix<double> _x; // k * n matrix
    dlib::matrix<double> _y; // 1 * n binary [0 or 1] responses
public:
    LogLik(const dlib::matrix<double>& x, const dlib::matrix<double>& y);

    double logLik(const colvec &theta);
    colvec logLik_d(const colvec &theta);
};