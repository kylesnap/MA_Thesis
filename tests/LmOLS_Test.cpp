//
// Created by Kyle Dewsnap on 2021-09-28.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "catch.hpp"
#include "../LmOLS.h"

double datX[] = { 1, 2.6200, 41.9847,
                  1, 2.8750, 38.2609,
                  1, 2.3200, 40.0862,
                  1, 3.2150, 34.2146,
                  1, 3.4400, 50.8721,
                  1, 3.4600, 30.3468,
                  1, 3.5700, 68.6275,
                  1, 3.1900, 19.4357,
                  1, 3.1500, 30.1587,
                  1, 3.4400, 35.7558,
                  1, 3.4400, 35.7558,
                  1, 4.0700, 44.2260,
                  1, 3.7300, 48.2574,
                  1, 3.7800, 47.6190,
                  1, 5.2500, 39.0476,
                  1, 5.4240, 39.6386,
                  1, 5.3450, 43.0309,
                  1, 2.2000, 30.0000,
                  1, 1.6150, 32.1981,
                  1, 1.8350, 35.4223,
                  1, 2.4650, 39.3509,
                  1, 3.5200, 42.6136,
                  1, 3.4350, 43.6681,
                  1, 3.8400, 63.8021,
                  1, 3.8450, 45.5137,
                  1, 1.9350, 34.1085,
                  1, 2.1400, 42.5234,
                  1, 1.5130, 74.6861,
                  1, 3.1700, 83.2808,
                  1, 2.7700, 63.1769,
                  1, 3.5700, 93.8375,
                  1, 2.7800, 39.2086 };

double datY[] = { 4.7619,
                  4.7619,
                  4.3860,
                  4.6729,
                  5.3476,
                  5.5249,
                  6.9930,
                  4.0984,
                  4.3860,
                  5.2083,
                  5.6180,
                  6.0976,
                  5.7803,
                  6.5789,
                  9.6154,
                  9.6154,
                  6.8027,
                  3.0864,
                  3.2895,
                  2.9499,
                  4.6512,
                  6.4516,
                  6.5789,
                  7.5188,
                  5.2083,
                  3.6630,
                  3.8462,
                  3.2895,
                  6.3291,
                  5.0761,
                  6.6667,
                  4.6729 };

TEST_CASE("Hand-Written Test") {
    double x[] = { 1,  0.7, 2.0,
                   1, -1.0, 0.4,
                   1 , 0.7 , 1.5,
                   1 , 0.6 , 0.3 };
    double y[] = { 0.4,
                   1.1,
                   -0.6,
                   1.0 };

    gsl_matrix_view X = gsl_matrix_view_array(x, 4, 3);
    gsl_vector_view Y = gsl_vector_view_array(y, 4);

    LmOLS mod = LmOLS(&X.matrix, &Y.vector);

    std::vector<float> betas;
    mod.getBetaHat(betas);
    REQUIRE_THAT(betas,
                 Catch::Matchers::Approx(std::vector<float>{1.071, -0.246, -0.510}).margin(0.001));

    std::vector<float> s;
    mod.getBetaSE(s);
    REQUIRE_THAT(s,
                 Catch::Matchers::Approx(std::vector<float>{0.861, 0.778, 0.778}).margin(0.001));

    REQUIRE(mod.getRSQ() == Catch::Detail::Approx(0.528).margin(0.001));
}

TEST_CASE("Test LM with Data") {

    gsl_matrix_view X = gsl_matrix_view_array(datX, 32, 3);
    gsl_vector_view y = gsl_vector_view_array(datY, 32);

    LmOLS mod = LmOLS(&X.matrix, &y.vector);

    std::vector<float> betas;
    mod.getBetaHat(betas);
    REQUIRE_THAT(betas,
                 Catch::Matchers::Approx(std::vector<float>{-0.401, 1.472, 0.02400}).margin(0.001));

    std::vector<float> s;
    mod.getBetaSE(s);
    REQUIRE_THAT(s,
                 Catch::Matchers::Approx(std::vector<float>{0.512, 0.122, 0.007}).margin(0.001));

    REQUIRE(mod.getRSQ() == Catch::Detail::Approx(0.848).margin(0.001));
}