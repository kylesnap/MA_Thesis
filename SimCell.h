// SimCell.h
// Defines CellParam and Cell class.
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>


#include "DesignMat.h"

#pragma once

class SimCell {
private:
    int const REPS = 10;
    
    float _rSq;
    DesignMat *_xMat;
    gsl_vector *_pTrue = gsl_vector_alloc(3);
public:
    SimCell(int n, float rSq, std::tuple<float, float> betaP, std::tuple<float, float> betaQ, std::tuple<float, float> propGrps,
            std::tuple<float, float, float> paramsTrue, gsl_rng *r);

    void getSimCell(std::vector<float> &v, bool print = false);
    // void run()
};