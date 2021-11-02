//
// Created by Kyle Dewsnap on 2021-10-18.
//

#include "catch.hpp"
#include "../DesignMat.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>

const gsl_rng_type *T = gsl_rng_default;
gsl_rng *r = gsl_rng_alloc(T);

BetaGen *dbP = new BetaGen(0.5, 50, r);
BetaGen *dbQ = new BetaGen(0.5, 50, r);

TEST_CASE("Exceptions") {
    REQUIRE_THROWS_AS(BetaGen(0.5, 1, r), std::out_of_range);

    REQUIRE_THROWS_AS(DesignMat(2,  dbP,  dbQ, r, 0.5, 0.5), std::length_error); // N == K
    REQUIRE_THROWS_AS(DesignMat(2,  dbP,  dbQ, r, 0.5, 0.5), std::length_error); // N < K

    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, -0.5, 0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 0.5, -0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 1.0, 0.5), std::out_of_range); // pA > 1
    REQUIRE_THROWS_AS(DesignMat(100,  dbP,  dbQ, r, 0.6, 0.5), std::out_of_range); // pA + pB > 1

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

TEST_CASE("Response Matrix Generation") {

    SECTION("Beta Generation Check") {
        gsl_vector *vecP = gsl_vector_calloc(10000);
        gsl_vector *vecR = gsl_vector_calloc(10000);

        // All ones!
        BetaGen all1 = BetaGen(2, 50, r);
        for (int i = 0; i < vecP->size; i++) {
            gsl_vector_set(vecP, i, all1.betaR());
            gsl_vector_set(vecR, i, all1.betaBernR());
        }
        REQUIRE(gsl_vector_sum(vecP) == vecP->size);
        REQUIRE(gsl_vector_equal(vecP, vecR));

        // All zeros!
        BetaGen all0 = BetaGen(-1, 50, r);
        for (int i = 0; i < vecP->size; i++) {
            gsl_vector_set(vecP, i, all0.betaR());
            gsl_vector_set(vecR, i, all0.betaBernR());
        }
        REQUIRE(gsl_vector_sum(vecP) == 0);
        REQUIRE(gsl_vector_equal(vecP, vecR));

        // Three Trials
        BetaGen tst2 = BetaGen(0.25, 50, r);
        BetaGen tst3 = BetaGen(0.5, 500, r);

        for (int i = 0; i < vecP->size; i++) {
            gsl_vector_set(vecP, i, dbP->betaR());
        }
        CHECK(gsl_stats_mean(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.5).margin(0.01));
        CHECK(gsl_stats_variance(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.0025).margin(0.01));

        for (int i = 0; i < vecP->size; i++) {
            gsl_vector_set(vecP, i, tst2.betaR());
        }
        CHECK(gsl_stats_mean(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.26).margin(0.01));
        CHECK(gsl_stats_variance(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.0038).margin(0.01));


        for (int i = 0; i < vecP->size; i++) {
            gsl_vector_set(vecP, i, tst3.betaR());
        }
        CHECK(gsl_stats_mean(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.5).margin(0.01));
        CHECK(gsl_stats_variance(vecP->data, vecP->stride, vecP->size) == Catch::Detail::Approx(0.0005).margin(0.01));
    }

    SECTION("First Rows Copy + Dimensionality Check") {
        DesignMat tst = DesignMat(100,  dbP,  dbQ, r, 0.75, 0.25);
        gsl_matrix *bad_x = gsl_matrix_calloc(101, 3);
        REQUIRE_THROWS(tst.fillResponses(bad_x));
        gsl_matrix *bad_x2 = gsl_matrix_calloc(100, 4);
        REQUIRE_THROWS(tst.fillResponses(bad_x2));
        gsl_matrix_free(bad_x);
        gsl_matrix_free(bad_x2);

        /*
        gsl_matrix *ex = gsl_matrix_calloc(100, 3);
        tst.fillResponses(ex);

        gsl_matrix_view tx2 = gsl_matrix_submatrix(tst.getTrueX(), 0, 0, 100, 2);
        gsl_matrix_view ex2 = gsl_matrix_submatrix(ex, 0, 0, 100, 2);

        REQUIRE(gsl_matrix_equal(&tx2.matrix, &ex2.matrix)); */
    }
}
