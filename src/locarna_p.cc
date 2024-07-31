/**
* \file locarna_p.cc
*
* \brief Defines main function of locarna_p
*
* LocARNA-P: global and LOCal Alignment of RNA - Partitition function variant
*
* Copyright (C) Sebastian Will 2005-
*
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "LocARNA/sequence.hh"
#include "LocARNA/basepairs.hh"
#include "LocARNA/aligner_p.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/edge_probs.hh"
#include "LocARNA/ribosum.hh"
#include "LocARNA/ribofit.hh"
#include "LocARNA/anchor_constraints.hh"
#include "LocARNA/trace_controller.hh"
#include "LocARNA/multiple_alignment.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/global_stopwatch.hh"
#include "LocARNA/main_helper.icc"

using namespace std;

//! Version string (from configure.ac via autoconf system)
const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

using namespace LocARNA;

// ------------------------------------------------------------
// Parameter

struct command_line_parameters
    : public MainHelper::std_command_line_parameters {
    bool write_arcmatch_probs;  //!< write_arcmatch_probs
    bool write_basematch_probs; //!< write_basematch_probs

    // ------------------------------------------------------------
    // File arguments

    std::string arcmatch_probs_file;  //!< arcmatch_probs_file
    std::string basematch_probs_file; //!< basematch_probs_file

    bool include_am_in_bm; //!< include_am_in_bm

    std::string fragment_match_probs; //!< fragment_match_probs

    double min_am_prob; //!< minimal arc match probability

    double min_bm_prob; //!< minimal arc match probability

    /**
     * @brief Scale for partition function
     *
     * All partition functions are multiplied by this factor.  Note
     * that this is much less flexible compared to the use of pf_scale
     * in RNAfold, which adapts the scale to the sequence length.
     *
     * @note IMPORTANT: here use double (not pf_score_t), since option
     * parser cannot interpret pf_score_t, when using long double; later
     * this will be converted to pf_score_t.
     * @note renamed to avoid conflict with vrna
     */
    double pf_scale;

    int temperature_alipf; //!< temperature for alignment partition functions
};

//! \brief holds command line parameters of locarna
command_line_parameters clp;

//! defines command line parameters
option_def my_options[] =
    {{"help", 'h', &clp.help, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["help"]},
     {"version", 'V', &clp.version, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["version"]},
     {"verbose", 'v', &clp.verbose, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["verbose"]},
     {"quiet", 'q', &clp.quiet, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["quiet"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Scoring parameters"},

     {"indel", 'i', 0, O_ARG_INT, &clp.indel, "-150", "score",
      clp.help_text["indel"]},
     {"indel-opening", 0, 0, O_ARG_INT, &clp.indel_opening, "-750", "score",
      clp.help_text["indel_opening"]},
     {"ribosum-file", 0, 0, O_ARG_STRING, &clp.ribosum_file, "RIBOSUM85_60",
      "f", clp.help_text["ribosum_file"]},
     {"use-ribosum", 0, 0, O_ARG_BOOL, &clp.use_ribosum, "true", "bool",
      clp.help_text["use_ribosum"]},
     {"match", 'm', 0, O_ARG_INT, &clp.match, "50", "score",
      clp.help_text["match"]},
     {"mismatch", 'M', 0, O_ARG_INT, &clp.mismatch, "0", "score",
      clp.help_text["mismatch"]},
     {"struct-weight", 's', 0, O_ARG_INT, &clp.struct_weight, "200", "score",
      clp.help_text["struct_weight"]},
     {"exp-prob", 'e', &clp.exp_prob_given, O_ARG_DOUBLE, &clp.exp_prob,
      O_NODEFAULT, "prob", clp.help_text["exp_prob"]},
     {"tau", 't', 0, O_ARG_INT, &clp.tau, "50", "factor", clp.help_text["tau"]},

     {"temperature-alipf", 0, 0, O_ARG_INT, &clp.temperature_alipf, "300",
      "int", clp.help_text["temperature_alipf"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Partition function representation"},

     {"pf-scale", 0, 0, O_ARG_DOUBLE, &clp.pf_scale, "1.0", "scale",
      clp.help_text["pf_scale"]},
     {"extended-pf", 0, &clp.extended_pf, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["extended_pf"]},
     {"quad-pf", 0, &clp.quad_pf, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["quad_pf"]
#if !defined(_GLIBCXX_USE_FLOAT128) || defined(__clang__)
      +" Quad precision (128 bit, __float128) is not available for your binary. Falls back to extended-pf."
#endif
     },

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Output"},

     {"write-arcmatch-probs", 0, &clp.write_arcmatch_probs, O_ARG_STRING,
      &clp.arcmatch_probs_file, O_NODEFAULT, "file",
      "Write arcmatch probabilities"},
     {"write-basematch-probs", 0, &clp.write_basematch_probs, O_ARG_STRING,
      &clp.basematch_probs_file, O_NODEFAULT, "file",
      "Write basematch probabilities"},
     {"min-am-prob", 'a', 0, O_ARG_DOUBLE, &clp.min_am_prob, "0.001", "amprob",
      clp.help_text["min_am_prob"]},
     {"min-bm-prob", 'b', 0, O_ARG_DOUBLE, &clp.min_bm_prob, "0.001", "bmprob",
      clp.help_text["min_bm_prob"]},
     {"include-am-in-bm", 0, &clp.include_am_in_bm, O_NO_ARG, 0, O_NODEFAULT,
      "", "Include arc match cases in base match probabilities"},
     {"stopwatch", 0, &clp.stopwatch, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["stopwatch"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "",
      "Heuristics for speed accuracy trade off"},

     {"min-prob", 'p', 0, O_ARG_DOUBLE, &clp.min_prob, "0.001", "prob",
      clp.help_text["min_prob"]},
     {"max-bps-length-ratio", 0, 0, O_ARG_DOUBLE, &clp.max_bps_length_ratio,
      "0.0", "factor", clp.help_text["max_bps_length_ratio"]},
     {"max-diff-am", 'D', 0, O_ARG_INT, &clp.max_diff_am, "-1", "diff",
      clp.help_text["max_diff_am"]},
     {"max-diff", 'd', 0, O_ARG_INT, &clp.max_diff, "-1", "diff",
      clp.help_text["max_diff"]},
     {"max-diff-at-am", 0, 0, O_ARG_INT, &clp.max_diff_at_am, "-1", "diff",
      clp.help_text["max_diff_at_am"]},
     {"max-diff-aln", 0, 0, O_ARG_STRING, &clp.max_diff_alignment_file, "",
      "aln file", clp.help_text["max_diff_alignment_file"]},
     {"max-diff-pw-aln", 0, 0, O_ARG_STRING, &clp.max_diff_pw_alignment, "",
      "alignment", clp.help_text["max_diff_pw_alignment"]},
     {"max-diff-relax", 0, &clp.max_diff_relax, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["max_diff_relax"]},
     {"min-trace-probability", 0, 0, O_ARG_DOUBLE, &clp.min_trace_probability,
      "1e-5", "probability", clp.help_text["min_trace_probability"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Computed probabilities"},

     {"fragment-match-probs", 0, 0, O_ARG_STRING, &clp.fragment_match_probs, "",
      "\"i j k l\"",
      "Requests probabilities for the match of fragments [i..j] and [k..l]. "
      "Accepts a ';' separated list of ranges."},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Constraints"},
     {"maxBPspan", 0, 0, O_ARG_INT, &clp.max_bp_span, "-1", "span",
      clp.help_text["max_bp_span"]},
     {"relaxed-anchors", 0, &clp.relaxed_anchors, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["relaxed_anchors"]},

     {"", 0, 0, O_SECTION_HIDE, 0, O_NODEFAULT, "", "Hidden Options"},

     {"ribofit", 0, 0, O_ARG_BOOL, &clp.ribofit, "false", "bool",
      clp.help_text["ribofit"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Input files"},

     {"", 0, 0, O_ARG_STRING, &clp.fileA, O_NODEFAULT, "Input 1",
      clp.help_text["fileA"]},
     {"", 0, 0, O_ARG_STRING, &clp.fileB, O_NODEFAULT, "Input 2",
      clp.help_text["fileB"]},

     {"", 0, 0, O_TEXT, 0, O_NODEFAULT, "", clp.help_text["files"]},

     {"", 0, 0, 0, 0, O_NODEFAULT, "", ""}};

// ------------------------------------------------------------

// ------------------------------------------------------------
// MAIN


/**
 * \brief Helper of main() of executable locarna_p
 *
 * @return success
 */
template <typename pf_score_t>
int
run_and_report();

/**
 * \brief Main function of executable locarna_p
 *
 * @param argc argument counter
 * @param argv argument vector
 *
 * @return success
 */
int
main(int argc, char **argv) {
    stopwatch.start("total");

    clp.no_lonely_pairs =
        false; //! @todo currently not a command line option of locarna_p

    // ------------------------------------------------------------
    // Process options

    bool process_success = process_options(argc, argv, my_options);

    if (clp.help) {
        cout << "locarna_p - pairwise partition function of "
             << "RNA alignments." << std::endl;
        cout << std::endl
             << "Computes base and base pair match probabilities " << std::endl
             << "from alignment partitition functions." << std::endl
             << std::endl;

        // cout << VERSION_STRING<<std::endl;

        print_help(argv[0], my_options);

        return 0;
    }

    if (clp.quiet) {
        clp.verbose = false;
    } // quiet overrides verbose

    if (clp.version || clp.verbose) {
        cout << "locarna_p (" << VERSION_STRING << ")" << std::endl;
        if (clp.version)
            return 0;
        else
            cout << std::endl;
    }

    if (!process_success) {
        std::cerr << O_error_msg << std::endl;
        print_usage(argv[0], my_options);
        return -1;
    }

    if (clp.stopwatch) {
        stopwatch.set_print_on_exit(true);
    }

    if (clp.verbose)
        print_options(my_options);

    if (clp.ribofit) {
        clp.use_ribosum = false;
    }

    if (clp.quad_pf) {
        return
            run_and_report<quad_pf_score_t>();
    } else if (clp.extended_pf) {
        return
            run_and_report<extended_pf_score_t>();
    } else {
        return
            run_and_report<standard_pf_score_t>();
    }
}


template <typename pf_score_t>
int
run_and_report() {
    typedef std::vector<int>::size_type size_type;

    check_score_t<pf_score_t> check(clp);

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    PFoldParams pfoldparams(PFoldParams::args::noLP(clp.no_lonely_pairs),
                            PFoldParams::args::stacking(clp.stacking),
                            PFoldParams::args::max_bp_span(clp.max_bp_span));

    std::unique_ptr<RnaData> rna_dataA;
    try {
        rna_dataA =
            std::make_unique<RnaData>(clp.fileA, clp.min_prob,
                                      clp.max_bps_length_ratio, pfoldparams);
    } catch (failure &f) {
        std::cerr << "ERROR: failed to read from file " << clp.fileA
                  << std::endl
                  << "       " << f.what() << std::endl;
        return -1;
    }

    std::unique_ptr<RnaData> rna_dataB;
    try {
        rna_dataB =
            std::make_unique<RnaData>(clp.fileB, clp.min_prob,
                                      clp.max_bps_length_ratio, pfoldparams);
    } catch (failure &f) {
        std::cerr << "ERROR: failed to read from file " << clp.fileB
                  << std::endl
                  << "       " << f.what() << std::endl;
        return -1;
    }

    const Sequence &seqA = rna_dataA->sequence();
    const Sequence &seqB = rna_dataB->sequence();

    size_type lenA = seqA.length();
    size_type lenB = seqB.length();

    // --------------------
    // handle max_diff restriction

    // missing: proper error handling in case that lenA, lenB, and
    // max_diff_alignment are incompatible

    // do inconsistency checking for max_diff_pw_alignment and
    // max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment != "" && clp.max_diff_alignment_file != "") {
        std::cerr
            << "Cannot simultaneously use both options --max-diff-pw-alignemnt"
            << " and --max-diff-alignment-file." << std::endl;
        return -1;
    }

    // construct TraceController and check inconsistency for with multiplicity
    // of sequences
    //

    std::unique_ptr<MultipleAlignment> multiple_ref_alignment = NULL;

    if (clp.max_diff_alignment_file != "") {
        multiple_ref_alignment =
            std::make_unique<MultipleAlignment>(clp.max_diff_alignment_file);
    } else if (clp.max_diff_pw_alignment != "") {
        if (seqA.num_of_rows() != 1 || seqB.num_of_rows() != 1) {
            std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning of "
                         "alignments."
                      << std::endl;
            return -1;
        }

        std::vector<std::string> alistr;
        split_at_separator(clp.max_diff_pw_alignment, '&', alistr);

        if (alistr.size() != 2) {
            std::cerr
                << "Invalid argument to --max-diff-pw-alignemnt; "
                << "require exactly one '&' separating the alignment strings."
                << std::endl;
            return -1;
        }

        if (alistr[0].length() != alistr[1].length()) {
            std::cerr << "Invalid argument to --max-diff-pw-alignemnt; "
                      << "alignment strings have unequal lengths." << std::endl;
            return -1;
        }

        multiple_ref_alignment =
            std::make_unique<MultipleAlignment>(seqA.seqentry(0).name(),
                                                seqB.seqentry(0).name(),
                                                alistr[0], alistr[1]);
    }

    // ----------------------------------------
    // Ribosum matrix
    //
    std::unique_ptr<RibosumFreq> ribosum;
    std::unique_ptr<Ribofit> ribofit;
    MainHelper::init_ribo_matrix(clp, ribosum, ribofit);


    AnchorConstraints seq_constraints(lenA, "", lenB, "", !clp.relaxed_anchors);

    TraceController trace_controller(seqA, seqB, multiple_ref_alignment.get(),
                                     clp.max_diff, clp.max_diff_relax);

    trace_controller.restrict_by_anchors(seq_constraints);

    restrict_trace_by_probabilities(clp, rna_dataA.get(), rna_dataB.get(),
                                    ribosum.get(), ribofit.get(),
                                    &trace_controller,pf_score_t());

    multiple_ref_alignment.release();

    // ----------------------------------------
    // construct set of relevant arc matches
    //

    // initialize from RnaData
    std::unique_ptr<ArcMatches> arc_matches =
        std::make_unique<ArcMatches>(*rna_dataA.get(), *rna_dataB.get(),
                                     clp.min_prob, clp.max_diff_am != -1
                                         ? (size_type)clp.max_diff_am
                                         : std::max(lenA, lenB),
                                     clp.max_diff_at_am != -1
                                         ? (size_type)clp.max_diff_at_am
                                         : std::max(lenA, lenB),
                                     trace_controller, seq_constraints);

    // ----------------------------------------
    // report on input in verbose mode
    if (clp.verbose)
        MainHelper::report_input(seqA, seqB, *arc_matches);

    // ----------------------------------------
    // construct scoring

    double my_exp_probA = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenA);
    double my_exp_probB = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenB);

    ScoringParams scoring_params(
        ScoringParams::match(clp.match),
        ScoringParams::mismatch(clp.mismatch),
        ScoringParams::indel(clp.indel),
        ScoringParams::indel_opening(clp.indel_opening),
        ScoringParams::ribosum(ribosum.get()),
        ScoringParams::ribofit(ribofit.get()),
        ScoringParams::struct_weight(clp.struct_weight),
        ScoringParams::tau_factor(clp.tau),
        ScoringParams::exp_probA(my_exp_probA),
        ScoringParams::exp_probB(my_exp_probB),
        ScoringParams::temperature_alipf(clp.temperature_alipf));

    PFScoring<pf_score_t> scoring(seqA, seqB, *rna_dataA, *rna_dataB, *arc_matches, nullptr,
                                  scoring_params);

    // ------------------------------------------------------------
    // Computation of the alignment score
    //

    using apparams_t =  AlignerPParams<pf_score_t>;
    // initialize aligner-p object, which does the alignment computation
    AlignerP<pf_score_t> aligner(
        apparams_t(AlignerParams::seqA(&seqA), AlignerParams::seqB(&seqB),
                   AlignerParams::scoring(&scoring),
                   AlignerParams::max_diff_am(clp.max_diff_am),
                   AlignerParams::max_diff_at_am(clp.max_diff_at_am),
                   AlignerParams::trace_controller(&trace_controller),
                   AlignerParams::constraints(&seq_constraints),
                   typename apparams_t::min_am_prob(clp.min_am_prob),
                   typename apparams_t::min_bm_prob(clp.min_bm_prob),
                   typename apparams_t::pf_scale((pf_score_t)clp.pf_scale)));

    if (clp.verbose) {
        std::cout << "Run inside algorithm." << std::endl;
    }

    pf_score_t pf = aligner.align_inside();

    if (!clp.quiet) {
        std::cout << "Partition function: " << pf << std::endl;
    }

    if (clp.verbose) {
        std::cout << "Run outside algorithm." << std::endl;
    }

    aligner.align_outside();

    if (clp.verbose) {
        std::cout << "Compute probabilities." << std::endl;
    }

    aligner.compute_arcmatch_probabilities();

    if (clp.write_arcmatch_probs) {
        if (clp.verbose) {
            std::cout << "Write Arc-match probabilities to file "
                      << clp.arcmatch_probs_file << "." << std::endl;
        }
        ofstream out(clp.arcmatch_probs_file.c_str());
        if (out.good()) {
            aligner.write_arcmatch_probabilities(out);
        } else {
            cerr << "Cannot write to " << clp.arcmatch_probs_file << "! Exit."
                 << std::endl;
            return -1;
        }
    }

    aligner.compute_basematch_probabilities(clp.include_am_in_bm);

    if (clp.write_basematch_probs) {
        if (clp.verbose) {
            std::cout << "Write Base-match probabilities to file "
                      << clp.basematch_probs_file << "." << std::endl;
        }
        ofstream out(clp.basematch_probs_file.c_str());
        if (out.good()) {
            aligner.write_basematch_probabilities(out);
        } else {
            cerr << "Cannot write to " << clp.basematch_probs_file << "! Exit."
                 << std::endl;
            return -1;
        }
    }

    // ----------------------------------------
    // optionally, compute probabilities that certain ranges are matched
    //

    if (clp.fragment_match_probs != "") {
        // parse the string fragment_match_probs
        // the string must have the format
        // fragment_list ::= empty | pos pos pos pos {; pos pos pos pos}

        std::istringstream fragments(clp.fragment_match_probs);

        size_type i;
        size_type j;
        size_type k;
        size_type l;

        while (!fragments.eof()) {
            if (!(fragments >> i >> j >> k >> l)) {
                std::cerr << "WARNING: expected 4 positions or end"
                          << " when parsing argument fragment-match-probs"
                          << "\"" << clp.fragment_match_probs << "\""
                          << std::endl;
                break;
            }

            std::cout << "FRAGMENT_MATCH_PROB " << i << " " << j << " " << k
                      << " " << l << " : ";
            std::cout.flush();
            std::cout << aligner.compute_fragment_match_prob(i, j, k, l)
                      << std::endl;

            if (!fragments.eof()) {
                char sep;
                fragments >> sep;
                if (sep != ';') {
                    std::cerr << "WARNING: expected ';' or end"
                              << " when parsing argument fragment-match-probs"
                              << std::endl;
                    break;
                }
            }
        }
    }

    stopwatch.stop("total");

    // DONE
    return 0;
}
