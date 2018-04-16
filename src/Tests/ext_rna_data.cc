#include "catch.hpp"

#include <fstream>
#include <sstream>

#include <../LocARNA/pfold_params.hh>
#include <../LocARNA/ext_rna_data.hh>

using namespace LocARNA;

/** @file some unit tests for ExtRnaData
*/

TEST_CASE("ExtRnaData can fold alignments, write to file and read again") {
    PFoldParams pfoldparams(PFoldParams::args::noLP(true),
                        PFoldParams::args::stacking(true));

    std::ostringstream sizeinfo1;
    std::ostringstream sizeinfo2;

    std::string outfilename = "ext-archaea.pp";

    SECTION("fold alignment") {
        std::unique_ptr<ExtRnaData> rna_data;
        REQUIRE_NOTHROW(rna_data =
                            std::make_unique<ExtRnaData>("archaea.aln", 0.01,
                                                         0.0001, 0.0001, 5, 10,
                                                         10, pfoldparams));

        rna_data->write_size_info(sizeinfo1);

        REQUIRE(sizeinfo1.str() ==
                "arcs: 13  stackings: 10  arcs in loops: 18  unpaireds in "
                "loops: 64");

        SECTION("write to pp file") {
            std::ofstream out(outfilename.c_str());
            REQUIRE(out.good());
            rna_data->write_pp(out);
            out.close();

            SECTION("and read again") {
                REQUIRE_NOTHROW(rna_data = std::make_unique<ExtRnaData>(outfilename, 0.01,
                                                                        0.0001, 0.0001, 5, 10,
                                                                        10, pfoldparams));

                rna_data->write_size_info(sizeinfo2);

                REQUIRE(sizeinfo1.str() == sizeinfo2.str());
            }
            std::remove(outfilename.c_str());
        }
    }
}

TEST_CASE(
    "ExtRnaData can initialize from fixed structure, even in the context of "
    "restricted maxBPspan") {
    SECTION("Write clustal-like input with fixed structure") {
        std::string filename = "test_ext_fixed_structure.pp";
        // write test example
        std::ofstream out(filename.c_str());

        out << "test CCCUCGGG" << std::endl << "#FS  ((...).)" << std::endl;

        out.close();

        SECTION(
            "read without maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(PFoldParams::args::noLP(false),
                                    PFoldParams::args::stacking(false));
            ExtRnaData rd(filename, 0.0, 0.0, 0.0, -1, -1, -1, pfoldparams);

            REQUIRE(rd.arc_prob(2, 6) > 0.99);
            REQUIRE(rd.arc_prob(1, 8) > 0.99);
            REQUIRE(rd.arc_prob(2, 7) < 0.01);

            REQUIRE(rd.unpaired_in_loop_prob(3, 2, 6) > 0.99);
            REQUIRE(rd.unpaired_in_loop_prob(3, 1, 8) < 0.01);

            REQUIRE(rd.arc_in_loop_prob(2, 6, 1, 8) > 0.99);
            REQUIRE(rd.arc_in_loop_prob(2, 7, 1, 8) < 0.01);

            REQUIRE(rd.arc_in_loop_prob(2, 6, 0, 9) < 0.01);
            REQUIRE(rd.arc_in_loop_prob(2, 7, 0, 9) < 0.01);
        }

        SECTION("read with maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(PFoldParams::args::noLP(false),
                                    PFoldParams::args::stacking(false),
                                    PFoldParams::args::max_bp_span(6));

            ExtRnaData rd(filename, 0.0, 0.0, 0.0, -1, -1, -1, pfoldparams);

            REQUIRE(rd.arc_prob(2, 6) > 0.99);
            REQUIRE(rd.arc_prob(1, 8) < 0.01);
            REQUIRE(rd.arc_prob(2, 7) < 0.01);

            REQUIRE(rd.unpaired_in_loop_prob(3, 2, 6) > 0.99);
            REQUIRE(rd.unpaired_in_loop_prob(8, 0, 9) > 0.99);

            REQUIRE(rd.arc_in_loop_prob(2, 6, 0, 9) > 0.99);
            REQUIRE(rd.arc_in_loop_prob(2, 7, 0, 9) < 0.01);
        }

        std::remove(filename.c_str());
    }
}
