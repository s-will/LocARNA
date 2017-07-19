#include "catch.hpp"

#include <cassert>
#include <vector>
#include <../LocARNA/zip.hh>

using namespace LocARNA;

/** @file some unit tests for zip and enumerate
*/

TEST_CASE("zip and enumerate work with std::vector") {
    const std::vector<size_t> v0 = {1,2,3,4};
    const std::vector<size_t> v1 = {1,2,3,4,5};
    const std::vector<size_t> v2 = {2,4,6,8};

    SECTION("zip works with containers of same length") {
        int i=0;
        for ( const auto &x : zip(v0,v2) ) {
            REQUIRE(x.first * 2 == x.second);
            i++;
        }
        REQUIRE(i == 4);
    }

    SECTION("zip works with containers of unequal length (first is longer)") {
        int i=0;
        for ( const auto &x : zip(v1,v2) ) {
            REQUIRE(x.first * 2 == x.second);
            i++;
        }
        REQUIRE(i == 4);

    }

    SECTION("zip works with containers of unequal length (second is longer)") {
        int i=0;
        for ( const auto &x : zip(v2,v1) ) {
            REQUIRE(x.first == 2 * x.second);
            i++;
        }
        REQUIRE(i == 4);

    }

    SECTION("enumerate works") {
        for ( const auto &x : enumerate(v2) ) {
            REQUIRE((x.first+1) * 2 == x.second);
        }
    }
}
