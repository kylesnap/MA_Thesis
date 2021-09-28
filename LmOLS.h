//
// Created by Kyle Dewsnap on 2021-09-27.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#pragma once

class LmOLS {
private:
    gsl_vector * _betaHat;
    gsl_matrix * _vCov;
    double _ssr;
public:
    LmOLS(const gsl_matrix *X, const gsl_vector *Y);
};