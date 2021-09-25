// SimCell.h
// Defines CellParam and Cell class.
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>

#pragma once

struct CellParam {
    int n = std::numeric_limits<int>::max();
    float testFloat = std::numeric_limits<float>::max();
};

class SimCell {
private:
    // TODO: Add random number generator object.
    const int _REP = 10;
    const int _K = 2;
    std::string _fileName;
    int _n;
    float _testFloat;
public:
    SimCell(CellParam p, std::string fileName);
    void run() const;
};