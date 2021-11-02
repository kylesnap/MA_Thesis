// SimCell.h
// Defines CellParam and Cell class.
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <gsl/gsl_rng.h>

#pragma once

struct CellParam {
    int n = std::numeric_limits<int>::max();
    float testFloat = std::numeric_limits<float>::max();
};

class SimCell {
private:
    const int REP = 1;
    const int K = 2;

    gsl_rng *_r;
    std::string _fileName;
    int _n;
    float _testFloat;
public:
    SimCell(CellParam p, gsl_rng *r, std::string fileName);
    void run() const;
};