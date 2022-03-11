// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include "SimCell.h"

SimCell::SimCell(int n, float varErr, std::tuple<float, float> betaM, std::tuple<float, float> betaF,
                 std::tuple<float, float> propGrps, std::tuple<float, float, float, float> paramsTrue, gsl_rng *r,
                 bool repeat) {
    // Builds the simulation cell.
    _REPS = repeat ? _REPS : 1;

    // Error checks for error variance.
    if (varErr <= 0) throw std::out_of_range("Variance of the Errors is smaller than zero.");
    _varErr = varErr;
    _r = r;

    // Construct objects -- Errors will be thrown by class methods.
    auto *bM = new BetaGen(std::get<0>(betaM), std::get<1>(betaM), r);
    auto *bF = new BetaGen(std::get<0>(betaF), std::get<1>(betaF), r);
    _xMat = new DesignMat(n, bM, bF, r, std::get<0>(propGrps), std::get<1>(propGrps));

    gsl_vector_set(_pTrue, 0, std::get<0>(paramsTrue));
    gsl_vector_set(_pTrue, 1, std::get<1>(paramsTrue));
    gsl_vector_set(_pTrue, 2, std::get<2>(paramsTrue));
    gsl_vector_set(_pTrue, 3, std::get<3>(paramsTrue));
}

void SimCell::toVec(std::vector<float> &v, bool print) {
    // Print SimCell params (All Params of Design Mat, true params, sigma).
    _xMat->getDesignMat(v);
    for (int i = 0; i < 4; i++) {
        v.push_back((float) gsl_vector_get(_pTrue, i));
    }
    v.push_back(_varErr);
    if (print) {
        for (float i : v) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
}

std::string SimCell::run() {
    // Runs the simulation.
    int n = (int) _xMat->getTX()->size1;

    // Make observed response scores.
    gsl_vector *tY = gsl_vector_calloc(n);
    for (int i = 0; i < tY->size; i++) {
        gsl_vector_set(tY, i, gsl_ran_gaussian(_r, _varErr));
    }
    gsl_blas_dgemv(CblasNoTrans, 1, _xMat->getTX(), _pTrue, 1, tY);

    // Begin simulation loop
    std::vector<float> v;
    std::string out;
    gsl_matrix *resp = gsl_matrix_alloc(n, 3);
    LmOLS *mod;
    for (int i = 0; i < _REPS; i++) {
        // Change responses based on satisficing and fit new model.
        int ns = _xMat->genResponses(resp);
        mod = new LmOLS(resp, tY);

        // Assemble line to print.
        v.clear();
        v.push_back((float) i);
        toVec(v);
        mod->getBetaHat(v);
        mod->getBetaSE(v);
        v.push_back(mod->getRSQ());
        mod->getBetaHat(v, true);
        mod->getBetaSE(v, true);
        v.push_back(mod->getRSQ(true));
        v.push_back((float) ns);

        gsl_vector_view newCount = gsl_matrix_column(resp, 2);
        v.push_back((float) gsl_vector_sum(&newCount.vector));

        for (float j : v) {
            out.append(std::to_string(j));
            out.append(",");
        }
        out.append("\n");
    }

    return out;
}