//
// Created by Kyle Dewsnap on 2021-10-18.
//

#include "catch.hpp"
#include "../DesignMat.h"
#include <gsl/gsl_rng.h>

const gsl_rng_type *T = gsl_rng_default;
gsl_rng *r = gsl_rng_alloc(T);

TEST_CASE("Exceptions") {

    REQUIRE_THROWS_AS(DesignMat(2, 2, r, 0.5, 0.5), std::length_error); // N == K
    REQUIRE_THROWS_AS(DesignMat(2, 3, r, 0.5, 0.5), std::length_error); // N < K

    REQUIRE_THROWS_AS(DesignMat(100, 3, r, -0.5, 0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100, 3, r, 0.5, -0.5), std::out_of_range); // pA < 0
    REQUIRE_THROWS_AS(DesignMat(100, 3, r, 1.0, 0.5), std::out_of_range); // pA > 1
    REQUIRE_THROWS_AS(DesignMat(100, 3, r, 0.6, 0.5), std::out_of_range); // pA + pB > 1

    REQUIRE_NOTHROW(DesignMat(100, 3, r, 0.5, 0.5));
    REQUIRE_NOTHROW(DesignMat(100, 3, r, 0.4, 0.3));
    REQUIRE_NOTHROW(DesignMat(100, 3, r, 1, 0));
}
