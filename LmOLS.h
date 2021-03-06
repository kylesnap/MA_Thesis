//
// Created by Kyle Dewsnap on 2021-09-27.
// Defines the linear model object
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

    gsl_vector * _betaHatN;
    gsl_vector * _betaSEN;
    double _rsqN;
public:
    LmOLS(const gsl_matrix *X, const gsl_vector *Y);

    void getBetaHat(std::vector<float> &v, bool norm = false);
    void getBetaSE(std::vector<float> &v, bool norm = false);

    float getRSQ(bool norm = false) const;
};