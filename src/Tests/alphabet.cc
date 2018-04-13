#include "catch.hpp"

#include <array>
#include <string>
#include <vector>
#include <../LocARNA/alphabet.hh>

using namespace LocARNA;

/** @file some unit tests for Alphabet
 */

TEST_CASE("Construct char alphabet from array is working") {

    std::array<char,4> arr { {'A','C','G','U'} };
    Alphabet<char,4> a(arr);

    REQUIRE( a.in('A') );
    REQUIRE( a.in('C') );
    REQUIRE( a.in('G') );
    REQUIRE( a.in('U') );
    REQUIRE( ! a.in('T') );

    REQUIRE( a.idx('A') == 0 );
    REQUIRE( a.idx('C') == 1 );
    REQUIRE( a.idx('G') == 2 );
    REQUIRE( a.idx('U') == 3 );
}

TEST_CASE("Construct char alphabet from string is working") {

    Alphabet<char,4> a("ACGU");

    REQUIRE( a.in('A') );
    REQUIRE( a.in('C') );
    REQUIRE( a.in('G') );
    REQUIRE( a.in('U') );
    REQUIRE( ! a.in('T') );

    REQUIRE( a.idx('A') == 0 );
    REQUIRE( a.idx('C') == 1 );
    REQUIRE( a.idx('G') == 2 );
    REQUIRE( a.idx('U') == 3 );
}

TEST_CASE("Construct string alphabet from vector is working") {

    std::vector<std::string> v = { "A", "C", "G", "U" };

    Alphabet<std::string,4> a(v);

    REQUIRE( a.in("A") );
    REQUIRE( a.in("C") );
    REQUIRE( a.in("G") );
    REQUIRE( a.in("U") );
    REQUIRE( ! a.in("T") );

    REQUIRE( a.idx("A") == 0 );
    REQUIRE( a.idx("C") == 1 );
    REQUIRE( a.idx("G") == 2 );
    REQUIRE( a.idx("U") == 3 );
}
