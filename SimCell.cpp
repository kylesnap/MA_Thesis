// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <utility>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "SimCell.h"

SimCell::SimCell(CellParam p, std::string fileName) {
    _fileName = fileName;
    _n = p.n;
    _testFloat = p.testFloat;
    if (_n == std::numeric_limits<int>::max() || _testFloat == std::numeric_limits<float>::max()) {
        throw std::invalid_argument("Missing cell parameters");
    }
}

void set_x(gsl_matrix *x, gsl_rng *r) {
    // Sets design matrix (deterministic)
    // Curently, both X1 and X2 ~ N(0, 1)
    for (int i = 0; i < x->size1; i++) {
        for (int j = 0; j < x->size2; j++) {
            double f = gsl_ran_ugaussian(r);
            gsl_matrix_set(x, i, j, f);
        }
    }
}

void SimCell::run() const {
    // Runs the simulation by opening an output stream, then loops the cell REPS number of times.
    std::streambuf *buf;
    std::ofstream of;

    // Declare matrices
    gsl_matrix *X = gsl_matrix_alloc(_n, _K); // n * k
    gsl_matrix *Y = gsl_matrix_alloc(_n, 1); // n * 1

    gsl_matrix *beta = gsl_matrix_alloc(_K + 1, 1);

    // Declare random generator
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    // Fill X variable with deterministic values
    set_x(X, r);

    // DEBUG: Print
    for (int i = 0; i < _n; i++)
        for (int j = 0; j < _K; j++)
            printf ("m(%d,%d) = %g\n", i, j,
                    gsl_matrix_get (X, i, j));

    for (int i = 0; i < _REP; i++) {
        // TODO: Implement this pseudocode
            // (1) Write Y variable using X and BETA combination (needs a method)
            // (2) Fit LM (use old python code to inform method)
            // (3) Write string methods for cell and OLS results
            // (4) Uncomment file IO.
    }

    // CLEAN IT UP
    gsl_matrix_free(X);
    gsl_matrix_free(Y);
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