//
// Created by Kyle Dewsnap on 2021-10-15.
//

#include "DesignMat.h"

DesignMat::DesignMat(int n, int k, gsl_rng *r, float pP, float pQ) {
    // Assigns groups to respondents and builds true design mat.

    // Basic error checking.
    if (n <= k) throw std::length_error("Too few observations for number of variables.");
    if ((pP - 1) * pP > 0 || (pQ - 1) * pQ > 0) {
        throw std::out_of_range("pP or pQ are outside a possible range [0, 1].");
    } else if (pP + pQ > 1) {
        throw std::out_of_range("pP and pQ sum to more than one.");
    }
    _n = n;
    _pP = pP;
    _pQ = pQ;

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

    tally_grps(); // Currently, this will just check for errors.
}

std::map<char, int> DesignMat::tally_grps() { // Verified
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

