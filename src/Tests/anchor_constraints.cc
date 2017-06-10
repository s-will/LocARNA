#include "catch.hpp"

#include <cassert>
#include <algorithm>
#include <../LocARNA/anchor_constraints.hh>

using namespace LocARNA;

TEST_CASE("Anchors (name length 1) can be constructed from strings with correct strict constraint semantics") {
    AnchorConstraints c(10,"A...B.C.EF",10,"..A.BDEF..",true);

    SECTION("allowed_match() works as expected") {
        REQUIRE(c.allowed_match(5,5));
        REQUIRE(!c.allowed_match(4,5));
        REQUIRE(!c.allowed_match(5,4));
        REQUIRE(!c.allowed_match(5,7));
        REQUIRE(!c.allowed_match(3,6));
        REQUIRE(!c.allowed_match(6,6)); // violates lex order on conraints 3 and 4
        REQUIRE(c.allowed_match(8,6));
        REQUIRE(!c.allowed_match(7,6));
    }

    SECTION("allowed_del() works as expected") {
        REQUIRE(c.allowed_del(6,5));
        REQUIRE(!c.allowed_del(7,6));
        REQUIRE(!c.allowed_del(5,5));
        REQUIRE(!c.allowed_del(6,4)); // 6 can be deleted 'right of' 4, but not *immediately right of* 4, because 5-5 is anchored!
        REQUIRE(!c.allowed_del(6,7));
        REQUIRE(!c.allowed_del(6,6));
    }


    SECTION("allowed_ins() works as expected") {
        REQUIRE(c.allowed_ins(7,6));
        REQUIRE(!c.allowed_ins(5,6));
    }

}
