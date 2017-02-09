#include "catch.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include <../LocARNA/sequence.hh>
#include <../LocARNA/multiple_alignment.hh>
#include <../LocARNA/trace_controller.hh>

using namespace LocARNA;

TEST_CASE("TracController is correctly initialized with simple alignment") {
    const std::string example_ma =
        "CLUSTAL W --- LocARNA 1.5.4 - Local Alignment of RNA\n"
        "\n"
        "\n"
        "seqA A--CTTG\n"
        "seqB ACCT--G\n";

    //   0 1 2 3 4 5
    // 0 *
    // 1   * * *
    // 2         *
    // 3         *
    // 4         *
    // 5          *

    //   0 1 2 3 4 5
    // 0 *
    // 1   * * *
    // 2         *
    // 3         *
    // 4         *
    // 5          *

    std::istringstream example_ma_istream(example_ma.c_str());

    MultipleAlignment ma(example_ma_istream);

    Sequence seqA;
    seqA.append(Sequence::SeqEntry("seqA", "A-CTTG"));

    Sequence seqB;
    seqB.append(Sequence::SeqEntry("seqB", "ACCT-G"));

    TraceController tc(seqA, seqB, &ma, 0);

    std::string expected_debug =
        "min_col_vector:   0   1   1   4   4   4   6 \n"
        "max_col_vector:   0   3   3   5   5   5   6 \n";

    std::ostringstream observed_debug;
    tc.print_debug(observed_debug);
    REQUIRE(observed_debug.str() == expected_debug);
}
