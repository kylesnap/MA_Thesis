// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include "SimCell.h"

SimCell::SimCell(int n, float rSq, std::tuple<float, float> betaP, std::tuple<float, float> betaQ,
                 std::tuple<float, float> propGrps, std::tuple<float, float, float> paramsTrue, gsl_rng *r) {

    // Error checks for RSQ
    if ((rSq - 1) * rSq >= 0) throw std::out_of_range("RSQ is a value outside of [0, 1]");
    _rSq = rSq;

    // Handles odd sample sizes silently.
    n = n%2 == 1 ? n - 1 : n;

    // Construct objects -- Errors will be thrown by class methods.
    auto *bP = new BetaGen(std::get<0>(betaP), std::get<1>(betaP), r);
    auto *bQ = new BetaGen(std::get<0>(betaQ), std::get<1>(betaQ), r);
    _xMat = new DesignMat(n, bP, bQ, r, std::get<0>(propGrps), std::get<1>(propGrps));

    gsl_vector_set(_pTrue, 0, std::get<0>(paramsTrue));
    gsl_vector_set(_pTrue, 1, std::get<1>(paramsTrue));
    gsl_vector_set(_pTrue, 2, std::get<2>(paramsTrue));
}

void SimCell::getSimCell(std::vector<float> &v, bool print) {
    // Print SimCell params (All Params of Design Mat, true params, RSQ).
    _xMat->getDesignMat(v);
    for (int i = 0; i < 3; i++) {
        v.push_back(gsl_vector_get(_pTrue, i));
    }
    v.push_back(_rSq);
    if (print) {
        for (float i : v) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
}

// TODO: Build run() starting with true Y creation, then loop.