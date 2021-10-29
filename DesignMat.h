//
// Created by Kyle Dewsnap on 2021-10-15.
//

# include <gsl/gsl_matrix.h>
# include <gsl/gsl_statistics_double.h>
# include <vector>
# include <tuple>
# include <iostream>
# include <algorithm>
# include <random>
# include <map>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_rng.h>

# pragma once

class DesignMat {
private:
    int _TK = 4; // Number of params in true model (B0, B1, BQ, BX)
    int _MK = 3; // Number of params in estimated model (B0, B1, BQ)

    gsl_matrix *_tX;
    std::vector<char>_grps;
    gsl_rng *_r;
    int _n;
    float _pP;
    float _pQ;
public:
    DesignMat(int n, gsl_rng *r, float pP, float pQ);

    std::map<char, int> tally_grps();
    std::vector<double> summary(bool head = false);
};