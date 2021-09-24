//
// Created by Kyle Dewsnap on 2021-09-20.
//

#include <dlib/matrix.h>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../LogLik.h"

TEST_CASE("Test logistic function") {
    SECTION("#1") {
        colvec t = {1, 1};
        colvec x = {1};

        CHECK(p(x, t) == Catch::Detail::Approx(0.8808));
    }

    SECTION("#2") {
        colvec t = {0, -1};
        colvec x = {0.5};

        CHECK(p(x, t) == Catch::Detail::Approx(0.3775407));
    }

    SECTION("#3") {
        colvec t = {0, 1, 6.9};
        colvec x = {2, -2};

        CHECK(p(x, t) == Catch::Detail::Approx(7.50450e-6));
    }
}

TEST_CASE("Dimensionality errors") {
    SECTION("N of X is not equal to N of Y") {
        dlib::matrix<double> x(1, 10);
        dlib::matrix<double> y(1, 10);
        dlib::matrix<double> x2(2, 20);
        dlib::matrix<double> y2(1, 20);

        CHECK_NOTHROW(LogLik(x, y));
        CHECK_NOTHROW(LogLik(x2, y2));

        CHECK_THROWS(LogLik(x, y2));
        CHECK_THROWS(LogLik(x2, y));
    }

    SECTION("K of Theta is not equal to K of X + 1") {
        dlib::matrix<double> x(1, 20);
        dlib::matrix<double> y(1, 20);

        colvec t(2, 1);
        colvec t2(1, 1); // One too few
        colvec t3(3, 1); // One too many

        LogLik test(x, y);

        CHECK_NOTHROW(test.logLik(t));
        CHECK_THROWS(test.logLik(t2));
        CHECK_THROWS(test.logLik(t3));
    }
}

TEST_CASE("Test LogLik") {
    SECTION("#1") {
        colvec t = {1, 1};
        dlib::matrix<double, 1, 1> x;
        dlib::matrix<double> y(1, 1);

        x = 1;
        y = 1;

        LogLik test(x, y);
        CHECK(test.logLik(t) == Catch::Detail::Approx(-0.126928));
    }

    SECTION("#2") {
        colvec t = {0.1, 1, 0.6};
        dlib::matrix<double, 2, 3> x;
        dlib::matrix<double> y(1, 3);

        x = 1.5, 4.5, 6.9,
        3.6, 9.0, -2.5;
        y = 1, 0, 1;

        LogLik test(x, y);
        CHECK(test.logLik(t) == Catch::Detail::Approx(-10.0271));
    }
}

TEST_CASE("Test LogLik Derivative") {
    SECTION("#1") {
        colvec t = {1, 1};
        dlib::matrix<double, 1, 1> x;
        dlib::matrix<double> y(1, 1);

        x = 1;
        y = 1;

        LogLik test(x, y);
        CHECK(test.logLik_d(t)(0) == Catch::Detail::Approx(0.1192029));
        CHECK(test.logLik_d(t)(1) == Catch::Detail::Approx(0.1192029));
    }

    SECTION("#2") {
        colvec t = {0.1, 1, 0.6};
        dlib::matrix<double, 2, 3> x;
        dlib::matrix<double> y(1, 3);

        x = 1.5, 4.5, 6.9,
                3.6, 9.0, -2.5;
        y = 1, 0, 1;

        LogLik test(x, y);
        CHECK(test.logLik_d(t)(0) == Catch::Detail::Approx(-0.9731305));
        CHECK(test.logLik_d(t)(1) == Catch::Detail::Approx(-4.43758));
        CHECK(test.logLik_d(t)(2) == Catch::Detail::Approx(-8.92785));
    }
}
