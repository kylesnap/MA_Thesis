//
// Created by Kyle Dewsnap on 2021-10-15.
//

#include "DesignMat.h"

DesignMat::DesignMat(int n, BetaGen *bP, BetaGen *bQ, gsl_rng *r, float pP, float pQ) {
    // Assigns groups to respondents and builds true design mat.

    // Basic error checking.
    if (n <= std::min(TK, MK)) throw std::length_error("Too few observations for number of variables.");
    if ((pP - 1) * pP > 0 || (pQ - 1) * pQ > 0) {
        throw std::out_of_range("pP or pQ are outside a possible range [0, 1].");
    } else if (pP + pQ > 1) {
        throw std::out_of_range("pP and pQ sum to more than one.");
    }
    _n = n;
    _pP = pP;
    _pQ = pQ;

    _bP = bP;
    _bQ = bQ;
    _r = r;

    // Build list of group memberships. // VERIFIED
    auto it = _grps.end();
    _grps.insert(it, ceil((float) _n * _pP), 'P');
    it = _grps.end();
    _grps.insert(it, ceil((float) _n * _pQ), 'Q');
    if (_grps.size() != _n) {
        it = _grps.end();
        _grps.insert(it, (_n - _grps.size()), 'X');
    }
    std::shuffle(_grps.begin(), _grps.end(), std::default_random_engine(gsl_rng_get(_r)));

    tallyGrps(); // Currently, this will just check for errors.

    // Build true X // VERIFIED
    _tX = gsl_matrix_calloc(_n, TK);

    gsl_vector_view row_i;
    char group_i;
    for (int i = 0; i < _n; i++) {
        row_i = gsl_matrix_row(_tX, i);
        gsl_vector_set(&row_i.vector, 0, 1);
        gsl_vector_set(&row_i.vector, 1, gsl_ran_ugaussian(r));
        group_i = _grps[i];
        switch (group_i) {
            case 'P':
                // Do nothing, as the P will be the reference group.
                break;
            case 'Q':
                gsl_vector_set(&row_i.vector, 2, 1);
                break;
            case 'X':
                gsl_vector_set(&row_i.vector, 3, 1);
                break;
            default: // Should be impossible, but worth a check.
                throw std::length_error("There are members with an invalid group.");
        }
    }

}

std::map<char, int> DesignMat::tallyGrps() { // Verified
    // Makes a tally of the number of group members
    int nP = 0, nQ = 0, nX = 0;
    for (char i : _grps) {
        switch (i) {
            case 'P':
                ++nP;
                break;
            case 'Q':
                ++nQ;
                break;
            case 'X':
                ++nX;
                break;
            default: // Should be impossible, but worth a check.
                throw std::length_error("There are members with an invalid group.");
        }
    }
    if (nP + nQ + nX != _n) throw std::length_error("Not every individual has an assigned group.");
    return { {'P', nP}, {'Q', nQ}, {'X', nX} };
}

std::vector<double> DesignMat::summary(bool head) { // TESTED!
    // Provides a summary of the mean and var of the X1, and the counts of Q and X. If head, print first six lines.
    gsl_vector_view col;
    col = gsl_matrix_column(_tX, 1);
    double mean = gsl_stats_mean(col.vector.data, col.vector.stride, col.vector.size);
    double var = gsl_stats_variance(col.vector.data, col.vector.stride, col.vector.size);

    col = gsl_matrix_column(_tX, 2);
    double q = gsl_vector_sum(&col.vector);

    col = gsl_matrix_column(_tX, 3);
    double x = gsl_vector_sum(&col.vector);

    if (head) {
        for (int i = 0; i < 6; i++) {
            col = gsl_matrix_row(_tX, i);
            std::cout << i << ": ";
            for (int j = 0; j < 4; j++) {
                std::cout << gsl_vector_get(&col.vector, j) << " ";
            }
            std::cout << std::endl;
        }
    }
    return {mean, var, q, x};
}

const gsl_matrix* DesignMat::getTX() {
    return _tX;
}

int DesignMat::genResponses(gsl_matrix *eX) { // Validated!
   // Takes in an empty design matrix, then copies tX data into it; however, column = 2 is binary responses to survey q.
   // Returns number of satisficiers.
   int ns = 0;
   if ((eX->size1 != _tX->size1) || (eX->size2 != MK)) throw std::length_error("Dimension mismatch.");

   // Copy first two cols
   gsl_vector_view subeX = gsl_matrix_column(eX, 0);
   gsl_vector_view subtX = gsl_matrix_column(_tX, 0);
   gsl_vector_memcpy(&subeX.vector, &subtX.vector);

   subeX = gsl_matrix_column(eX, 1);
   subtX = gsl_matrix_column(_tX, 1);
   gsl_vector_memcpy(&subeX.vector, &subtX.vector);

   // Fill third col (which represents Qs) with responses
   subeX = gsl_matrix_column(eX, 2); // Column of item responses
   char grp_i;
   for (int i = 0; i < _tX->size1; i++) {
       subtX = gsl_matrix_subrow(_tX, i, 2, 2); // Subrow of i'th group membership
       grp_i = gsl_vector_isnull(&subtX.vector) ? 'P' : (gsl_vector_max_index(&subtX.vector) == 0) ? 'Q' : 'X';

       switch (grp_i) {
           case 'P':
               if (_bP->betaBernR()) {
                   gsl_vector_set(&subeX.vector, i, gsl_ran_bernoulli(_r, 0.5));
                   ns++;
               } else {
                   gsl_vector_set(&subeX.vector, i, 0);
               }
               break;
           case 'Q':
               if (_bQ->betaBernR()) {
                   gsl_vector_set(&subeX.vector, i, gsl_ran_bernoulli(_r, 0.5));
                   ns++;
               } else {
                   gsl_vector_set(&subeX.vector, i, 1);
               }
               break;
           case 'X':
               // Always satisfice.
               gsl_vector_set(&subeX.vector, i, gsl_ran_bernoulli(_r, 0.5));
               ++ns;
               break;
           default: // Should be impossible, but worth a check.
               throw std::length_error("There are members with an invalid group.");
       }
   }

   return ns;
}

void DesignMat::getDesignMat(std::vector<float> &v) {
    // Appends parameters of design matrix (N; Abs. Count of P, Q, and X; Beta Dists for P and Q.)
    v.push_back((float) _n);
    std::map<char, int> grps = tallyGrps();
    v.push_back((float) grps['P']);
    v.push_back((float) grps['Q']);
    v.push_back((float) grps['X']);
    _bP->getBetaDist(v);
    _bQ->getBetaDist(v);
}

BetaGen::BetaGen(float mode, float conc, gsl_rng *r, bool print) { // VALIDATED
    if (conc < 2) throw std::out_of_range("Concentration below 2.");
    _r = r;
    _mode = mode;
    _conc = conc;

    if (mode < 0) {
        if (print) std::cout << "Mode below zero: All variables generated will have p = 0." << std::endl;
        _bh = -1;
        return;
    } else if (mode > 1) {
        if (print) std::cout << "Mode above one: All variables generated will have p = 1" << std::endl;
        _bh = 1;
        return;
    }
    _a = _mode * (_conc - 2) + 1;
    _b = _conc - _a;
}

double BetaGen::betaR() { // VALIDATED
    // Return a single double from the beta distribution.
    return (_bh == -1) ? 0 : (_bh == 1) ? 1 : gsl_ran_beta(_r, _a, _b);
}

int BetaGen::betaBernR() { // VALIDATED
    return (_bh == -1) ? 0 : (_bh == 1) ? 1 : gsl_ran_bernoulli(_r, gsl_ran_beta(_r, _a, _b)); // NOLINT(cppcoreguidelines-narrowing-conversions)
}

void BetaGen::getBetaDist(std::vector<float> &v) { //
    // Prints A and B for distribution (A and B are both negative if beh == -1, both 0 if beh == 1
    (_bh == -1) ? v.push_back(-1) : (_bh == 1) ? v.push_back(0) : v.push_back(_a);
    (_bh == -1) ? v.push_back(-1) : (_bh == 1) ? v.push_back(0) : v.push_back(_b);
}