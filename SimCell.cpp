// SimCell.cpp
// Implements
// Kyle Dewsnap
// 16SEP21

#include "SimCell.h"

SimCell::SimCell(int n, float varErr, std::tuple<float, float> betaP, std::tuple<float, float> betaQ,
                 std::tuple<float, float> propGrps, std::tuple<float, float, float, float> paramsTrue, gsl_rng *r) {

    // Error checks for RSQ
    if (varErr <= 0) throw std::out_of_range("Variance of the Errors is smaller than zero.");
    _varErr = varErr;
    _r = r;

    // Handles odd sample sizes silently.
    n = n%2 == 1 ? n - 1 : n;

    // Construct objects -- Errors will be thrown by class methods.
    auto *bP = new BetaGen(std::get<0>(betaP), std::get<1>(betaP), r);
    auto *bQ = new BetaGen(std::get<0>(betaQ), std::get<1>(betaQ), r);
    _xMat = new DesignMat(n, bP, bQ, r, std::get<0>(propGrps), std::get<1>(propGrps));

    gsl_vector_set(_pTrue, 0, std::get<0>(paramsTrue));
    gsl_vector_set(_pTrue, 1, std::get<1>(paramsTrue));
    gsl_vector_set(_pTrue, 2, std::get<2>(paramsTrue));
    gsl_vector_set(_pTrue, 3, std::get<3>(paramsTrue));
}

void SimCell::toVec(std::vector<float> &v, bool print) {
    // Print SimCell params (All Params of Design Mat, true params, sigma).
    _xMat->getDesignMat(v);
    for (int i = 0; i < 4; i++) {
        v.push_back(gsl_vector_get(_pTrue, i));
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
    int n = (int) _xMat->getTX()->size1;
    // Build a vector of 'true' Y responses.

    // Add error terms to tY.
    gsl_vector *tY = gsl_vector_calloc(n);
    for (int i = 0; i < tY->size; i++) {
        gsl_vector_set(tY, i, gsl_ran_gaussian(_r, _varErr));
    }

    /*
    // Build tY = X * Params + e
    gsl_blas_dgemv(CblasNoTrans, 1.0, _xMat->getTX(), _pTrue, 1.0, tY);
    LmOLS chk = LmOLS(_xMat->getTX(), tY);

    // Testing
    chk.getBetaHat(v);
    chk.getBetaSE(v);

    for (float i : v) {
        std::cout << i << " ";
    }
    std::cout << chk.getRSQ() << std::endl;
    std::cout << "=" << std::endl;

     */
    // Begin simulation loop
    std::vector<float> v;
    std::string out;
    gsl_matrix *resp = gsl_matrix_alloc(n, 3);
    LmOLS *mod;
    std::cout << "I,N,NP,NQ,NX,BP_A,BP_B,BQ_A,BQ_B,"
                 "BTRUE_0,BTRUE_1,BTRUE_Q,BTRUE_X,ERR_VAR,"
                 "BHAT_0,BHAT_1,BHAT_Q,BSE_0,BSE_1,BSE_Q,RSQ" << std::endl;
    for (int i = 0; i < REPS; i++) {
        // Change responses based on satisficing and fit new model.
        _xMat->genResponses(resp);
        mod = new LmOLS(resp, tY);

        // Assemble line to print.
        v.clear();
        v.push_back((float) i);
        toVec(v);
        mod->getBetaHat(v);
        mod->getBetaSE(v);
        v.push_back(mod->getRSQ());

        for (float j : v) {
            out.append(std::to_string(j));
            out.append(",");
        }
        out.append("\n");
    }

    return out;

    // TODO: Final testing!
}