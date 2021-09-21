//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../LogLik.h"

TEST_CASE("Tests for LogLik Function") {
    SECTION("#1") {

        dlib::matrix<double, 2, 2> x;
        dlib::matrix<double, 2, 1> y;
        colvec t;

        x = 1, 1,
        0, 0;
        y = 1,
        0;
        t = 1,
        1;

        LogLik test = LogLik(x, y);
        REQUIRE(test.logLik(t) == -0.8200752);
    }
}