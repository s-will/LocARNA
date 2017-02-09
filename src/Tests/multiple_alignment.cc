#include "catch.hpp"

#include <iostream>
#include <sstream>
#include <../LocARNA/multiple_alignment.hh>
#include <../LocARNA/alignment.hh>
#include <../LocARNA/sequence.hh>

using namespace LocARNA;

/** @file some unit tests for MultipleAlignment and Sequence

    Tests writing from fasta and clustalw format, basic functionality
    and simple checks
*/

TEST_CASE("a multiple alignment can be constructed as empty") {
    MultipleAlignment ma;

    REQUIRE(ma.empty());
}

TEST_CASE("2-way multiple alignment can be created from strings") {
    MultipleAlignment ma("seqA", "seqB", "A-CGT-U", "CCCG-CU");

    REQUIRE(ma.is_proper());
}

TEST_CASE("a multiple alignment can be created from file (clustalw)") {
    MultipleAlignment *ma = 0L;

    ma = new MultipleAlignment("archaea.aln");

    SECTION("the alignment is properly constructed, as intended") {
        REQUIRE(ma->is_proper());
        REQUIRE(!ma->empty());
        REQUIRE(ma->num_of_rows() == 6);
        REQUIRE(ma->contains("fdhA"));
        REQUIRE(ma->contains("fruA"));
    }

    SECTION("a sequence object can be created from the multiple alignment") {
        Sequence seq = (*ma).as_sequence();
        delete ma;

        SECTION(
            "a further alignment string can be appended, after deleting the "
            "original multiple alignment") {
            std::string name_str = "hdrA";
            std::string seq_str =
                "GG--CACCACUCGAAGGCUA-------------AG-CCAAAGUGGUG--CU";
            seq.append(Sequence::SeqEntry(name_str, seq_str));

            REQUIRE(seq.is_proper());
            REQUIRE(7 == seq.num_of_rows());
            REQUIRE(seq.contains("fdhA"));
            REQUIRE(seq.contains(name_str));

            size_t index = seq.index(name_str);
            std::string seq_str_at_index = seq.seqentry(index).seq().str();
            REQUIRE(seq_str_at_index == seq_str);
        }
    }
}

TEST_CASE("creating multiple alignment from file with wrong format fails") {
    MultipleAlignment ma("archaea-wrong.fa",
                         MultipleAlignment::FormatType::FASTA);

    REQUIRE(!ma.is_proper());
}

TEST_CASE("multiple alignment can be constructed from file (fasta)") {
    MultipleAlignment *ma;

    ma = new MultipleAlignment("archaea-aln.fa",
                               MultipleAlignment::FormatType::FASTA);
    REQUIRE(ma->is_proper());
    REQUIRE(!ma->empty());

    SECTION("the multiple alignment can be converted to a sequence") {
        Sequence seq = (*ma).as_sequence();
        delete ma;

        SECTION(
            "the sequence is proper and has 6 entries, after deleting the "
            "multiple alignment") {
            //! test whether seq is proper
            REQUIRE(seq.is_proper());

            //! test whether seq has correct size
            REQUIRE(seq.num_of_rows() == 6);
        }

        SECTION("the sequence can be written and read again (clustalw)") {
            std::ostringstream out;
            out << "CLUSTAL W" << std::endl;
            seq.write(out);

            std::istringstream in(out.str());
            try {
                MultipleAlignment seq2(in,
                                       MultipleAlignment::FormatType::CLUSTAL);

                REQUIRE(seq2.is_proper());
                REQUIRE(seq.num_of_rows() == seq2.num_of_rows());
                REQUIRE(seq.length() == seq2.length());

            } catch (failure &f) {
                std::cerr << "Catched exception: " << f.what() << std::endl;
                std::cerr << "... while reading" << std::endl
                          << out.str() << std::endl;
                REQUIRE(false);
            }
        }

        SECTION("the sequence can be written and read again (stockholm)") {
            std::ostringstream out;

            out << "# STOCKHOLM 1.0" << std::endl;
            seq.write(out, MultipleAlignment::FormatType::STOCKHOLM);

            std::istringstream in(out.str());
            try {
                MultipleAlignment
                    seq2(in, MultipleAlignment::FormatType::STOCKHOLM);

                REQUIRE(seq2.is_proper());
                REQUIRE(seq.num_of_rows() == seq2.num_of_rows());
                REQUIRE(seq.length() == seq2.length());

            } catch (failure &f) {
                std::cerr << "Catched exception: " << f.what() << std::endl;
                std::cerr << "... while reading" << std::endl
                          << out.str() << std::endl;
                REQUIRE(false);
            }
        }
    }
}
