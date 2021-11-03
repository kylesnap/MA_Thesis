//
// Created by Kyle Dewsnap on 2021-11-03.
//

#include "catch.hpp"
#include "../SimCell.h"

const gsl_rng_type *T2 = gsl_rng_default;
gsl_rng *r2 = gsl_rng_alloc(T2);

TEST_CASE("Error Throws") {
    std::tuple<float, float> goodBeta = {0.5, 50};
    std::tuple<float, float> badBeta = {0.5, 1};

    std::tuple<float, float> goodP = {0.5, 0.4};
    std::tuple<float, float> badP = {0.5, 0.6};

    std::tuple<float, float, float, float> goodTheta = {0, 0, 0, 0};

    // Bad N or Rsq
    REQUIRE_THROWS_AS(SimCell(2, 0.8, goodBeta, goodBeta, goodP, goodTheta, r2), std::length_error);
    REQUIRE_THROWS_AS(SimCell(100, 1.2, goodBeta, goodBeta, goodP, goodTheta, r2), std::out_of_range);

    // Bad Beta Params
    REQUIRE_THROWS_AS(SimCell(100, 0.8, badBeta, goodBeta, goodP, goodTheta, r2), std::out_of_range);
    REQUIRE_THROWS_AS(SimCell(100, 0.8, goodBeta, badBeta, goodP, goodTheta, r2), std::out_of_range);

    // Bad Group Props.
    REQUIRE_THROWS_AS(SimCell(100, 0.8, goodBeta, goodBeta, badP, goodTheta, r2), std::out_of_range);

    // All Good
    SimCell *tst;
    std::vector<float> v;
    REQUIRE_NOTHROW(tst = new SimCell(100, 0.8, goodBeta, goodBeta, goodP, goodTheta, r2));
}

TEST_CASE("Test Design Mat Construction") {
    std::tuple<float, float> b1 = {0.5, 50};
    std::tuple<float, float> b2 = {0.25, 40};

    std::tuple<float, float> p1 = {0.5, 0.5};
    std::tuple<float, float> p2 = {0.25, 0.4};

    std::tuple<float, float, float, float> t1 = {0, 0, 0, 0};
    std::tuple<float, float, float, float> t2 = {1, 0, 1, 0};

    std::vector<float> v;
    SimCell tst1 = SimCell(100, 0.8, b1, b1, p1, t1, r2);
    tst1.toVec(v);
    // N = 100, RSQ = 0.8, (A, B) = (25, 25), PROP = (50, 50, 0) THETA = 0, 0, 0
    CHECK_THAT(v, Catch::Matchers::Equals(std::vector<float>{100, 50, 50, 0, 25, 25, 25, 25, 0, 0, 0, 0, 0.8}));

    v.clear();
    SimCell tst2 = SimCell(200, 0.5, b2, b2, p2, t2, r2);
    tst2.toVec(v);
    // N = 200, RSQ = 0.5, (A, B) = (10.5, 29.5), PROP = (50, 80, 70), THETA = 1, 0, 1
    CHECK_THAT(v, Catch::Matchers::Equals(std::vector<float>{200, 50, 80, 70, 10.5, 29.5, 10.5, 29.5, 1, 0, 1, 0, 0.5}));

    v.clear();
    SimCell tst3 = SimCell(420, 0.69, b1, b2, p1, t2, r2);
    tst3.toVec(v);
    // N = 420, RSQ = 0.69, BQ/BP are different, PROP = (210, 210, 0), THETA = 1, 0, 1
    CHECK_THAT(v, Catch::Matchers::Equals(std::vector<float>{420, 210, 210, 0, 25, 25, 10.5, 29.5, 1, 0, 1, 0, 0.69}));
}