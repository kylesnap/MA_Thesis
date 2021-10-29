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

TEST_CASE("Member Generation") {

    SECTION("Single Member Groups") {

        DesignMat all_p = DesignMat(100, 3, r, 1, 0);
        DesignMat all_q = DesignMat(50, 3, r, 0, 1);
        DesignMat all_x = DesignMat(75, 3, r, 0, 0);

        REQUIRE(all_p.tally_grps()['P'] == 100);
        REQUIRE(all_q.tally_grps()['Q'] == 50);
        REQUIRE(all_x.tally_grps()['X'] == 75);

    }

    SECTION("Complicated Group Splits") {

        DesignMat tst = DesignMat(100, 3, r, 0.75, 0.25); // Even N, equal p and q, no x
        REQUIRE(tst.tally_grps()['P'] == 75);
        REQUIRE(tst.tally_grps()['Q'] == 25);
        REQUIRE(tst.tally_grps()['X'] == 0);

        tst = DesignMat(97, 3, r, 0.33, 0.33); // Odd N, equal p and q, add x
        REQUIRE(tst.tally_grps()['P'] == 33);
        REQUIRE(tst.tally_grps()['Q'] == 33);
        REQUIRE(tst.tally_grps()['X'] == 31);

        tst = DesignMat(79, 3, r, 0.69, 0.0); // Odd N, no q, add x
        REQUIRE(tst.tally_grps()['P'] == 55);
        REQUIRE(tst.tally_grps()['Q'] == 0);
        REQUIRE(tst.tally_grps()['X'] == 24);
    }


}
