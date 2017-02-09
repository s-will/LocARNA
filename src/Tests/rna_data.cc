#include "catch.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include <../LocARNA/pfold_params.hh>
#include <../LocARNA/sequence.hh>
#include <../LocARNA/alignment.hh>
#include <../LocARNA/rna_ensemble.hh>
#include <../LocARNA/rna_data.hh>

using namespace LocARNA;

/** @file some unit tests for RnaData and ExtRnaData
*/

PFoldParams pfparams(true, true, -1, 2);

TEST_CASE("RnaData can fold alignments, write to file and reead again") {
    std::ostringstream sizeinfo1;
    std::ostringstream sizeinfo2;

    SECTION("fold alignment") {
        RnaData *rna_data = 0L;
        REQUIRE_NOTHROW(rna_data =
                            new RnaData("archaea.aln", 0.1, 2, pfparams));
        rna_data->write_size_info(sizeinfo1);

        SECTION("write to pp file") {
            std::ofstream out("archaea.pp");
            REQUIRE(out.is_open());

            REQUIRE_NOTHROW(rna_data->write_pp(out));

            SECTION("and read again") {
                if (rna_data)
                    delete rna_data;
                REQUIRE_NOTHROW(
                    rna_data = new RnaData("archaea.pp", 0.1, 2, pfparams));

                rna_data->write_size_info(sizeinfo2);

                REQUIRE(sizeinfo1.str() == sizeinfo2.str());
            }
            std::remove("archaea.pp");
        }
        delete rna_data;
    }
}

TEST_CASE("RnaData can construct pairwise (averaged) consensus dot plots") {
    SECTION("create pairwise alignment") {
        std::string alistrA =
            "CCUCG--AGGGGAACCCGA-------------AAGGGACC-CGA-GAGG";
        std::string alistrB =
            "CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGG-AGCAGGUGGCG";

        Sequence seqA("seqA", "CCUCGAGGGGAACCCGAAAGGGACCCGAGAGG");
        Sequence seqB("seqB",
                      "CGCCACCCUGCGAACCCAAUAUAAAAUAAUACAAGGGAGCAGGUGGCG");

        Alignment alignment(seqA, seqB,
                            Alignment::edges_t(Alignment::alistr_to_edge_ends(
                                                   alistrA),
                                               Alignment::alistr_to_edge_ends(
                                                   alistrB)));

        REQUIRE(alignment.seqA().is_proper());
        REQUIRE(alignment.seqB().is_proper());

        SECTION("and compute the averaged consensus dot plot") {
            RnaEnsemble ensA(seqA, pfparams, false, true);
            RnaData rna_dataA(ensA, 0.05, 3, pfparams);
            RnaEnsemble ensB(seqB, pfparams, false, true);
            RnaData rna_dataB(ensB, 0.05, 3, pfparams);

            RnaData consensus(rna_dataA, rna_dataB, alignment, 0.01, 0.01);

            // missing requirements for consensus

            SECTION("write to file and read again") {
                std::string filename = "test_avg_dp.pp";
                std::ofstream out(filename.c_str());

                consensus.write_pp(out);
                out.close();

                RnaData consensus2(filename, 0.1, 2, pfparams);

                std::remove(filename.c_str());

                std::ostringstream sizeinfo1;
                std::ostringstream sizeinfo2;
                consensus.write_size_info(sizeinfo1);
                consensus.write_size_info(sizeinfo2);

                REQUIRE(sizeinfo1.str() == sizeinfo2.str());
            }
        }

        SECTION(
            "RnaEnsemble can compute the alifold consensus dotplot of the "
            "alignment") {
            MultipleAlignment ma(alignment);

            RnaEnsemble ens_ma(ma, pfparams, false, true);

            RnaData ali_consensus(ens_ma, 0.1, 0, pfparams);

            // missing test of ensemble

            SECTION("this can be written to file and read again by RnaData") {
                std::string filename = "test_alicons_dp.pp";
                // write and read again for the ali_consensus object
                std::ofstream out(filename.c_str());
                // std::cout<<"Consensus: "<<std::endl;

                ali_consensus.write_pp(out);
                out.close();

                RnaData ali_consensus2(filename, 0.1, 2, pfparams);
                std::remove(filename.c_str());

                std::ostringstream sizeinfo1;
                std::ostringstream sizeinfo2;
                ali_consensus.write_size_info(sizeinfo1);
                ali_consensus2.write_size_info(sizeinfo2);

                REQUIRE(sizeinfo1.str() == sizeinfo2.str());
            }
        }
    }
}

TEST_CASE(
    "RnaData can initialize from fixed structure, even in the context of "
    "restricted maxBPspan") {
    SECTION("Write clustal-like input with fixed structure") {
        std::string filename = "test_fixed_structure.pp";
        // write test example
        std::ofstream out(filename.c_str());

        out << "test CCCUCGGG" << std::endl << "#FS  ((...).)" << std::endl;

        out.close();

        SECTION(
            "read without maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(false, false, -1, 2);
            RnaData rd(filename, 0.0, 0.0, pfoldparams);

            REQUIRE(rd.arc_prob(2, 6) > 0.99);
            REQUIRE(rd.arc_prob(1, 8) > 0.99);
            REQUIRE(rd.arc_prob(2, 7) < 0.01);
        }

        SECTION("read with maxBPspan restriction and check some base pairs") {
            PFoldParams pfoldparams(false, false, 6, 2);
            RnaData rd(filename, 0.0, 0.0, pfoldparams);

            REQUIRE(rd.arc_prob(2, 6) > 0.99);
            REQUIRE(rd.arc_prob(1, 8) < 0.01);
            REQUIRE(rd.arc_prob(2, 7) < 0.01);
        }

        std::remove(filename.c_str());
    }
}
