//
// Created by Kyle Dewsnap on 2021-10-18.
//

#include "catch.hpp"
#include "../DesignMat.h"
#include <gsl/gsl_rng.h>

const gsl_rng_type *T = gsl_rng_default;
gsl_rng *r = gsl_rng_alloc(T);

BetaP dbP = {0.5, 50};
BetaP dbQ = {0.5, 50};

TEST_CASE("Exceptions") {

    REQUIRE_THROWS_AS(DesignMat(2,  dbP,  dbQ, r, 0.5, 0.5), std::length_error); // N == K
    REQUIRE_THROWS_AS(DesignMat(2,  dbP,  dbQ, r, 0.5, 0.5), std::length_error); // N < K

    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, -0.5, 0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 0.5, -0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 1.0, 0.5), std::out_of_range); // pA > 1
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 0.6, 0.5), std::out_of_range); // pA + pB > 1

    BetaP invalid_bP = {2, 50};
    BetaP invalid_bQ = {0.5, 1};
    REQUIRE_THROWS_AS(DesignMat(100,  invalid_bP,  dbQ, r, 0.5, 0.5), std::out_of_range);
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  invalid_bQ, r, 0.5, 0.5), std::out_of_range);

    REQUIRE_NOTHROW(DesignMat(100,  dbP,  dbQ, r, 0.5, 0.5));
    REQUIRE_NOTHROW(DesignMat(100,  dbP,  dbQ, r, 0.4, 0.3));
    REQUIRE_NOTHROW(DesignMat(100,  dbP,  dbQ, r, 1, 0));
}

TEST_CASE("Member Generation") {

    SECTION("Single Member Groups") {

        DesignMat all_p = DesignMat(100,  dbP,  dbQ, r, 1, 0);
        DesignMat all_q = DesignMat(50,  dbP,  dbQ, r, 0, 1);
        DesignMat all_x = DesignMat(75,  dbP,  dbQ, r, 0, 0);

        REQUIRE(all_p.tallyGrps()['P'] == 100);
        REQUIRE(all_q.tallyGrps()['Q'] == 50);
        REQUIRE(all_x.tallyGrps()['X'] == 75);

    }

    SECTION("Complicated Group Splits") {

        DesignMat tst = DesignMat(100,  dbP,  dbQ, r, 0.75, 0.25); // Even N, equal p and q, no x
        REQUIRE(tst.tallyGrps()['P'] == 75);
        REQUIRE(tst.tallyGrps()['Q'] == 25);
        REQUIRE(tst.tallyGrps()['X'] == 0);

        tst = DesignMat(97,  dbP,  dbQ, r, 0.33, 0.33); // Odd N, equal p and q, add x
        REQUIRE(tst.tallyGrps()['P'] == 33);
        REQUIRE(tst.tallyGrps()['Q'] == 33);
        REQUIRE(tst.tallyGrps()['X'] == 31);

        tst = DesignMat(79,  dbP,  dbQ, r, 0.69, 0.0); // Odd N, no q, add x
        REQUIRE(tst.tallyGrps()['P'] == 55);
        REQUIRE(tst.tallyGrps()['Q'] == 0);
        REQUIRE(tst.tallyGrps()['X'] == 24);
    }

}

TEST_CASE("True Design Matrix Generation") {
    std::vector<double> res;

    DesignMat tst = DesignMat(100000,  dbP,  dbQ, r, 0.75, 0.25);
    res = {0, 1, (double) tst.tallyGrps()['Q'], (double) tst.tallyGrps()['X']};
    REQUIRE_THAT(tst.summary(), Catch::Matchers::Approx(res).margin(0.01));

    tst = DesignMat(97000,  dbP,  dbQ, r, 0.33, 0.33);
    res = {0, 1, (double) tst.tallyGrps()['Q'], (double) tst.tallyGrps()['X']};
    REQUIRE_THAT(tst.summary(), Catch::Matchers::Approx(res).margin(0.01));

    tst = DesignMat(79000,  dbP,  dbQ, r, 0.69, 0.0); // Odd N, no q, add x
    res = {0, 1, (double) tst.tallyGrps()['Q'], (double) tst.tallyGrps()['X']};
    REQUIRE_THAT(tst.summary(), Catch::Matchers::Approx(res).margin(0.01));
}

TEST_CASE("Response Matrix Generation") { // Verified

    SECTION("Check Satisficing Rates") {
        // Honest P and Q
        BetaP t1bP = {0, 500};
        BetaP t1bQ = {0, 500};

        DesignMat tst = DesignMat(10000, t1bP, t1bQ, r, 0.5, 0.5);
        gsl_matrix *err_mat = gsl_matrix_alloc(1000, 3);
        REQUIRE_THROWS(tst.fillResponses(err_mat));

        gsl_matrix *to_fill = gsl_matrix_alloc(10000, 3);
        REQUIRE(tst.fillResponses(to_fill) == 0); // No one should satisfice.

        // Honest P, no Qs, but Xs will satisfice.
        tst = DesignMat(10000, t1bP, t1bQ, r, 0.5, 0.0);
        REQUIRE(tst.fillResponses(to_fill) == 10000 / 2);

        // Everyone lies!
        BetaP t2bP = {0.99999, 500000};
        BetaP t2bQ = {0.99999, 500000};
        tst = DesignMat(10000, t2bP, t2bQ, r, 0.5, 0.5);
        REQUIRE(tst.fillResponses(to_fill) == 10000);

        // Only P and X lie; Q will be true.
        tst = DesignMat(10000, t2bP, t1bQ, r, 1.0 / 3.0, 1.0 / 3.0);
        REQUIRE(tst.fillResponses(to_fill) == Approx(2 * (10000.0 / 3.0)).margin(1));
    }

}
