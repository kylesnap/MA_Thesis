// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include <iostream>
#include <string>
#include <utility>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include "SimCell.h"
#include "LmOLS.h"

SimCell::SimCell(CellParam p, gsl_rng *r, std::string fileName) {
    _r = r;
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
    gsl_matrix_set_all(x, 1); // ALWAYS set first col. to ones.
    for (int i = 0; i < x->size1; i++) {
        for (int j = 1; j < x->size2; j++) { // Starts at j = 1 because j = 0 represent the constant term.
            double f = gsl_ran_ugaussian(r);
            gsl_matrix_set(x, i, j, f);
        }
    }
}

void set_e(gsl_vector *e, gsl_rng *r) {
    // Sets the error terms
    // Currently, just returns a vector of zeros
    gsl_vector_set_all(e, 0);
}

void SimCell::run() const {
    // Runs the simulation by opening an output stream, then loops the cell REPS number of times.
    std::streambuf *buf;
    std::ofstream of;

    // Prepare the design matrices for the true and estimated models


    // Declare matrices, vectors, and intercepts
    gsl_matrix *X = gsl_matrix_alloc(_n, _K + 1); // n * k + 1 (Left most column are constant terms)

    gsl_vector *y = gsl_vector_alloc(_n); // n * 1
    gsl_vector *beta = gsl_vector_alloc(_K + 1); // k * 1

    // Fill X and beta with initial values
    set_x(X, _r);
    gsl_vector_set_all(beta, 1);

    LmOLS *mod;
    std::string ln;
    for (int i = 0; i < _REP; i++) {
        set_e(y, _r); // Error terms are first held in y
        gsl_blas_dgemv(CblasNoTrans, 1.0, X, beta, 1.0, y); // y = (X * Beta) + y
        mod = new LmOLS(X, y);
    }

    // TODO: Add string methods to write csv lines to output

    /* for (int i = 0; i < X->size1; i++) {
        printf("%f\n",
               gsl_matrix_get(X, i, i)
        );
    }
    gsl_vector_fprintf(stdout, y, "%f"); */

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