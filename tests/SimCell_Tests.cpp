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

    std::tuple<float, float, float> goodTheta = {0, 0, 0};

    // Bad N or Rsq
    REQUIRE_THROWS_AS(SimCell(2, 0.8, goodBeta, goodBeta, goodP, goodTheta, r2), std::length_error);
    REQUIRE_THROWS_AS(SimCell(100, 1.2, goodBeta, goodBeta, goodP, goodTheta, r2), std::out_of_range);

    // Bad Beta Params
    REQUIRE_THROWS_AS(SimCell(100, 0.8, badBeta, goodBeta, goodP, goodTheta, r2), std::out_of_range);
    REQUIRE_THROWS_AS(SimCell(100, 0.8, goodBeta, badBeta, goodP, goodTheta, r2), std::out_of_range);

    // Bad Group Props.
    REQUIRE_THROWS_AS(SimCell(100, 0.8, goodBeta, goodBeta, badP, goodTheta, r2), std::out_of_range);

    // All Good
    REQUIRE_NOTHROW(SimCell(100, 0.8, goodBeta, goodBeta, goodP, goodTheta, r2));
}