//
// Created by Kyle Dewsnap on 2021-10-15.
// Defines the design matrix object, and the beta-bern val gen.
//

# include <gsl/gsl_matrix.h>
# include <gsl/gsl_statistics_double.h>
# include <vector>
# include <tuple>
# include <iostream>
# include <algorithm>
# include <random>
# include <cmath>
# include <map>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_rng.h>

# pragma once

class BetaGen {
private:
    float _mode;
    float _conc;
    gsl_rng *_r;

    int _bh = 0;
    float _a = 0;
    float _b = 0;
public:
    BetaGen(float mode, float conc, gsl_rng *r, bool print = false);

    double betaR();
    int betaBernR();
    void getBetaDist(std::vector<float> &v) const;
};

class DesignMat {
private:
    int TK = 4; // Number of params in true model (B0, B1, BQ, BX)
    int MK = 3; // Number of params in estimated model (B0, B1, BQ)

    gsl_matrix *_tX;
    std::vector<char>_grps;
    gsl_rng *_r;
    int _n;
    BetaGen *_bM;
    BetaGen *_bF;
    float _pM;
    float _pF;
public:
    DesignMat(int n, BetaGen *bM, BetaGen *bF, gsl_rng *r, float pM, float pF);

    std::map<char, int> tallyGrps();
    std::vector<double> summary(bool head = false);
    void getDesignMat(std::vector<float> &v);

    const gsl_matrix * getTX();

    int genResponses(gsl_matrix *eX);
};