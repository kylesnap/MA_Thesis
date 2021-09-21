// SimCell.h
// Defines CellParam and Cell class.
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <dlib/matrix.h>

#include "LogLik.h"

#pragma once

struct CellParam {
    int n = std::numeric_limits<int>::max();
    float b1 = std::numeric_limits<float>::max();
};

class SimCell {
private:
    const int _REP = 10;
    std::string _fileName;
    int _n;
    float _b1;
public:
    SimCell(CellParam p, std::string fileName);
    double logLik(const colvec &theta);
    colvec logLik_d(const colvec &theta);
    void run();
};