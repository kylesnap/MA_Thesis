//
// Created by Kyle Dewsnap on 2021-10-15.
//

# include <gsl/gsl_matrix.h>
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
    gsl_matrix *_trueX;
    std::vector<char>_grps;
    gsl_rng *_r;
    int _n;
    float _pP;
    float _pQ;
public:
    DesignMat(int n, int k, gsl_rng *r, float pP, float pQ);

    std::map<char, int> tally_grps();

};