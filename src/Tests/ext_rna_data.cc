#include "catch.hpp"

#include <fstream>
#include <sstream>

#include <../LocARNA/pfold_params.hh>
#include <../LocARNA/ext_rna_data.hh>

using namespace LocARNA;

/** @file some unit tests for ExtRnaData
*/


TEST_CASE("ExtRnaData can fold alignments, write to file and read again") {
    PFoldParams pfparams(true,true,-1,2);

    std::ostringstream sizeinfo1;
    std::ostringstream sizeinfo2;

    std::string outfilename="ext-archaea.pp";

    SECTION("fold alignment") {
        ExtRnaData *rna_data=0L;
        REQUIRE_NOTHROW( rna_data =
                         new ExtRnaData("archaea.aln",0.01,0.0001,0.0001,5,10,10,pfparams) );

        rna_data->write_size_info(sizeinfo1);

        REQUIRE(sizeinfo1.str()
                == "arcs: 13  stackings: 10  arcs in loops: 18  unpaireds in loops: 64");

        SECTION("write to pp file") {
            std::ofstream out(outfilename.c_str());
            REQUIRE(out.good());
            rna_data->write_pp(out);
            out.close();

            SECTION("and read again") {
                if (rna_data) delete rna_data;
                REQUIRE_NOTHROW( rna_data =
                                 new ExtRnaData(outfilename,
                                                0.01,0.0001,0.0001,5,10,10,pfparams) );

                rna_data->write_size_info(sizeinfo2);

                REQUIRE(sizeinfo1.str() == sizeinfo2.str());
            }
            std::remove(outfilename.c_str());
        }
        delete rna_data;
    }
}


TEST_CASE("ExtRnaData can initialize from fixed structure, even in the context of restricted maxBPspan") {
    SECTION("Write clustal-like input with fixed structure") {

        std::string filename="test_ext_fixed_structure.pp";
        // write test example
        std::ofstream out(filename.c_str());

        out
            << "test CCCUCGGG"<<std::endl
            << "#FS  ((...).)"<<std::endl;

        out.close();

        SECTION("read without maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(false,false,-1,2);
            ExtRnaData rd(filename,
                          0.0, 0.0, 0.0,
                          -1, -1, -1, pfoldparams);

            REQUIRE( rd.arc_prob(2,6) > 0.99 );
            REQUIRE( rd.arc_prob(1,8) > 0.99 );
            REQUIRE( rd.arc_prob(2,7) < 0.01 );

            REQUIRE( rd.unpaired_in_loop_prob(3,2,6) > 0.99 );
            REQUIRE( rd.unpaired_in_loop_prob(3,1,8) < 0.01 );

            REQUIRE( rd.arc_in_loop_prob(2,6,1,8) > 0.99 );
            REQUIRE( rd.arc_in_loop_prob(2,7,1,8) < 0.01 );

            REQUIRE( rd.arc_in_loop_prob(2,6,0,9) < 0.01 );
            REQUIRE( rd.arc_in_loop_prob(2,7,0,9) < 0.01 );

        }

        SECTION("read with maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(false,false,6,2);
            ExtRnaData rd(filename,
                          0.0, 0.0, 0.0,
                          -1, -1, -1, pfoldparams);

            REQUIRE( rd.arc_prob(2,6) > 0.99 );
            REQUIRE( rd.arc_prob(1,8) < 0.01 );
            REQUIRE( rd.arc_prob(2,7) < 0.01 );

            REQUIRE( rd.unpaired_in_loop_prob(3,2,6) > 0.99 );
            REQUIRE( rd.unpaired_in_loop_prob(8,0,9) > 0.99 );

            REQUIRE( rd.arc_in_loop_prob(2,6,0,9) > 0.99 );
            REQUIRE( rd.arc_in_loop_prob(2,7,0,9) < 0.01 );
        }

        std::remove(filename.c_str());
    }
}
