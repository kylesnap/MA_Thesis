// SimCell.h
// Defines CellParam and Cell class.
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "DesignMat.h"
#include "LmOLS.h"

#pragma once

class SimCell {
private:
    int _REPS = 1000;
    
    float _varErr;
    gsl_rng *_r;
    DesignMat *_xMat;
    gsl_vector *_pTrue = gsl_vector_alloc(4);
public:
    SimCell(int n, float varErr, std::tuple<float, float> betaM, std::tuple<float, float> betaF, std::tuple<float, float> propGrps,
            std::tuple<float, float, float, float> paramsTrue, gsl_rng *r, bool repeat = true);

    void toVec(std::vector<float> &v, bool print = false);
    std::string run();
};