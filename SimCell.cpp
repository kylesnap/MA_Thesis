// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <utility>
#include <fstream>
#include <dlib/matrix.h>

#include "SimCell.h"

SimCell::SimCell(CellParam p, std::string fileName) {
    _fileName = fileName;
    _n = p.n;
    _b1 = p.b1;
    if (_n == std::numeric_limits<int>::max() || _b1 == std::numeric_limits<float>::max()) {
        throw std::invalid_argument("Missing cell parameters");
    }
}

void SimCell::run() {
    // Runs the simulation by opening an output stream, then loops the cell REPS number of times.
    std::streambuf *buf;
    std::ofstream of;

    // Make design matrix
    dlib::matrix<double> x(_n, 2);
    dlib::set_colm(x, 0) = 1;
    dlib::set_colm(x, 1) = dlib::randm(_n, 1);
    std::cout << x << std::endl;

    // Make design matrix
    for (int i = 0; i < _REP; i++) {
        ;
    }
/*
    // Opens a buffer!
    if(_fileName[0] != 0) {
        of.open(_fileName, std::ios::out | std::ios::app);
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream out(buf);



    if (of.is_open()) of.close();*/
}