//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../LogLik.h"

TEST_CASE("One dimensional tests") {
    colvec theta;
    theta = 1, 2;

    dlib::matrix<double, 1, 2> x_1, x_2;
    x_1 = 0, 0;
    x_2 = 5, 6;

    dlib::matrix<double> y;

    LogLik t1 = LogLik(x_1, y);
    LogLik t2 = LogLik(x_2, y);

    REQUIRE(t1.logLik(theta) == 0);
    REQUIRE(t2.logLik(theta) == 17);
}

TEST_CASE("Two dimensional tests") {
    colvec theta;
    theta = 1, 1;

    dlib::matrix<double, 2, 2> x_1, x_2;
    x_1 = 0, 0,
    0, 0;
    x_2 = 1, 2,
    3, 4;

    dlib::matrix<double> y;

    dlib::matrix<double, 2, 1> r_1, r_2;
    r_1 = 0,
    0;
    r_2 = 3,
    7;

    LogLik t1 = LogLik(x_1, y);
    LogLik t2 = LogLik(x_2, y);

    CHECK(t1.logLik(theta) == r_1);
    CHECK(t2.logLik(theta) == r_2);
}