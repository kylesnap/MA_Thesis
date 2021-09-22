//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../LogLik.h"

TEST_CASE("Tests for LogLik Function") {
    dlib::matrix<double, 2, 2> x;
    dlib::matrix<double, 2, 1> y;
    colvec t;

    SECTION("#1") {

        x = 1,1,
        0,0;
        y = 1,
        0;
        t = 1, 1;

        LogLik test = LogLik(x, y);
        CHECK(test.logLik(t) == Catch::Detail::Approx(-0.8200752));
    }

    SECTION("#2") {

        x = 1.1,0.2,
        -6.9,4.20;
        y = 0,
        1;
        t = 0.25, 7.5;

        LogLik test = LogLik(x, y);
        REQUIRE(test.logLik(t) == Catch::Detail::Approx(-1.931562));
    }

    dlib::matrix<double, 3, 2> x2;
    dlib::matrix<double, 3, 1> y2;
    SECTION("#3") {

        x2 = 1,23,
        1,22,
        1,21;
        y2 = 1,
        0,
        1;
        t = 1, -1;

        LogLik test = LogLik(x2, y2);
        REQUIRE(test.logLik(t) == Catch::Detail::Approx(-42));
    }
}

TEST_CASE("Tests for LogLik Derivative Function") {
    dlib::matrix<double, 2, 2> x;
    dlib::matrix<double, 2, 1> y;
    colvec t;

    SECTION("#1") {

        x = 1,1,
        0,0;
        y = 1,
        0;
        t = 1, 1;

        LogLik test = LogLik(x, y);
        test.logLik_d(t);
    }

    SECTION("#2") {

        x = 1.1,0.2,
                -6.9,4.20;
        y = 0,
                1;
        t = 0.25, 7.5;

        LogLik test = LogLik(x, y);
        test.logLik_d(t);
    }

    dlib::matrix<double, 3, 2> x2;
    dlib::matrix<double, 3, 1> y2;
    SECTION("#3") {

        x2 = 1,23,
                1,22,
                1,21;
        y2 = 1,
                0,
                1;
        t = 1, -1;

        LogLik test = LogLik(x2, y2);
        test.logLik_d(t);
    }
}