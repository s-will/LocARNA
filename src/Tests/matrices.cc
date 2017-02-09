#include "catch.hpp"

#include <cassert>
#include <algorithm>
#include <../LocARNA/matrices.hh>

using namespace LocARNA;

/** @file some unit tests for the Matrix classes
*/

struct mul2 {
    size_t
    operator()(size_t x) {
        return 2 * x;
    }
};

TEST_CASE("Matrix can be resized, filled, and transformed") {
    size_t x = 3;
    size_t y = 4;

    Matrix<size_t> m;
    m.resize(x, y);

    m.fill(2);

    for (size_t i = 0; i < x; i++) {
        for (size_t j = 0; j < y; j++) {
            m(i, j) += i + 2 * j;
        }
    }

    m.transform(mul2());

    bool reread_ok = true;
    for (size_t i = 0; i < x; i++) {
        for (size_t j = 0; j < y; j++) {
            reread_ok &= (m(i, j) == 2 * (2 + i + 2 * j));
        }
    }
    REQUIRE(reread_ok);
}

TEST_CASE("OMatrix can be filled and read again") {
    size_t x = 3;
    size_t y = 4;
    size_t xo = 2;
    size_t yo = 1;

    SECTION("case xdim<ydim") {
        OMatrix<size_t> m;

        m.resize(x, y, xo, yo);

        for (size_t i = xo; i < xo + x; i++) {
            for (size_t j = yo; j < yo + y; j++) {
                m(i, j) = i * (j + 1);
            }
        }

        bool reread_ok = true;
        for (size_t i = xo; i < xo + x; i++) {
            for (size_t j = yo; j < yo + y; j++) {
                reread_ok &= (m(i, j) == i * (j + 1));
            }
        }
        REQUIRE(reread_ok);
    }
    SECTION("case xdim>ydim") {
        OMatrix<size_t> m;
        std::swap(x, y);

        m.resize(x, y, xo, yo);

        for (size_t i = xo; i < xo + x; i++) {
            for (size_t j = yo; j < yo + y; j++) {
                m(i, j) = i * (j + 1);
            }
        }

        bool reread_ok = true;
        for (size_t i = xo; i < xo + x; i++) {
            for (size_t j = yo; j < yo + y; j++) {
                reread_ok &= (m(i, j) == i * (j + 1));
            }
        }
        REQUIRE(reread_ok);
    }
}
