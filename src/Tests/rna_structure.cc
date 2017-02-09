#include "catch.hpp"

#include <iostream>
#include <../LocARNA/rna_structure.hh>

using namespace LocARNA;


extern "C" {
#  include <ViennaRNA/energy_const.h> // import TURN
}


//                            12345678901234567890
std::string structure_string="..(((......)))......";
//                               123456789012345678901234567
std::string structure_pk_string=".((((.[[[.{{..)))).]]].}}.";


TEST_CASE("dot bracket string without pk can be parsed") {
    RnaStructure structure(structure_string);

    REQUIRE(structure.to_string() == structure_string);

    SECTION("dot bracket string with pk can be parsed") {
        RnaStructure structure_pk(structure_pk_string);

        REQUIRE(structure_pk.to_string() == structure_pk_string);

        SECTION("structures with pk can be copied and compared") {
            RnaStructure structure_pk_copy(structure_pk);
            REQUIRE(structure_pk == structure_pk_copy);
            REQUIRE(!( structure == structure_pk_copy ));
        }

        SECTION("nested and non-nested structures are distinguished") {
            REQUIRE(structure.nested());
            REQUIRE(!structure_pk.nested());
        }

        SECTION("crossing structures are recognized") {
            REQUIRE(structure.crossing());
            REQUIRE(structure_pk.crossing());
        }
    }

    SECTION("unlimited structures are recognized") {
        SECTION("structure is unlimited if left ends touch") {
            structure.insert(RnaStructure::bp_t(3,9));
            REQUIRE(!structure.crossing());
        }
        SECTION("structure is unlimited if right ends touch") {
            structure.insert(RnaStructure::bp_t(7,12));
            REQUIRE(!structure.crossing());
        }
        SECTION("structure is unlimited if left end touches right end") {
            structure.insert(RnaStructure::bp_t(14,19));
            REQUIRE(!structure.crossing());
        }
    }
}
