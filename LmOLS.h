//
// Created by Kyle Dewsnap on 2021-09-27.
//

#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>

#pragma once

class LmOLS {
private:
    gsl_vector * _betaHat;
    gsl_vector * _betaSE;
    double _rsq;
public:
    LmOLS(const gsl_matrix *X, const gsl_vector *Y);

    void getBetaHat(std::vector<float> &v);
    void getBetaSE(std::vector<float> &v);

    float getRSQ() const;
};