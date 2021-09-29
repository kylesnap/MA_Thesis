//
// Created by Kyle Dewsnap on 2021-09-28.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "catch.hpp"
#include "../LmOLS.h"

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

    std::vector<double> betas;
    mod.getBetaHat(betas);
    REQUIRE_THAT(betas,
                 Catch::Matchers::Approx(std::vector<double>{1.071, -0.246, -0.510}).margin(0.001));

    std::vector<double> bses;
    mod.getBetaSE(bses);
    REQUIRE_THAT(bses,
                 Catch::Matchers::Approx(std::vector<double>{0.800, 0.723, 0.723}).margin(0.001));

    REQUIRE(mod.getRSQ() == Catch::Detail::Approx(0.528).margin(0.001));
}