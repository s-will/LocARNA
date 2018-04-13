/**
 * \file sparse.cc
 *
 * \brief Defines main function of SPARSE
 *
 * Copyright (C) Milad Miladi <miladim(@)informatik.uni-freiburg.de>
 */

#include <iostream>
#include <fstream>
#include <vector>

//#include <math.h>

#include "LocARNA/sequence.hh"
#include "LocARNA/basepairs.hh"
#include "LocARNA/alignment.hh"
#include "LocARNA/aligner_n.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/edge_probs.hh"
#include "LocARNA/ribosum.hh"
#include "LocARNA/ribofit.hh"
#include "LocARNA/anchor_constraints.hh"
#include "LocARNA/sequence_annotation.hh"
#include "LocARNA/trace_controller.hh"
#include "LocARNA/multiple_alignment.hh"
#include "LocARNA/sparsification_mapper.hh"
#include "LocARNA/global_stopwatch.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/main_helper.icc"
#include "LocARNA/aligner_params.hh"

using namespace std;
using namespace LocARNA;

//! Version string (from configure.ac via autoconf system)
const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

// ------------------------------------------------------------
// Parameter

// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

//! \brief Switch on/off trace back
//! @note never made it into command line
const bool DO_TRACE = true;

//! \brief Structure for command line parameters of locarna
//!
//! Encapsulating all command line parameters in a common structure
//! avoids name conflicts and makes downstream code more informative.
//!
struct command_line_parameters
    : public MainHelper::std_command_line_parameters,
      public MainHelper::mea_command_line_parameters {
    double max_uil_length_ratio;  // max unpaired in loop length ratio
    double max_bpil_length_ratio; // max base pairs in loop length ratio

    int indel_loop; //!< indel extension score

    int indel_opening_loop; //!< indel opening score for loops

    double prob_unpaired_in_loop_threshold; //!< threshold for
                                            //! prob_unpaired_in_loop
    double prob_basepair_in_loop_threshold; //!< threshold for
                                            //! prob_basepait_in_loop

    bool special_gap_symbols; //!< whether to use special gap
                              //! symbols in the alignment result
};

//! \brief holds command line parameters of locarna
command_line_parameters clp;

//! defines command line parameters
option_def my_options[] = {
    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "cmd_only"},

    {"help", 'h', &clp.help, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["help"]},
    {"galaxy-xml", 0, &clp.galaxy_xml, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["galaxy_xml"]},
    {"version", 'V', &clp.version, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["version"]},
    {"verbose", 'v', &clp.verbose, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["verbose"]},
    {"quiet", 'q', &clp.quiet, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["quiet"]},

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Scoring parameters"},

    {"indel", 'i', 0, O_ARG_INT, &clp.indel, "-150", "score",
     clp.help_text["indel"]},
    {"indel-loop", 'i', 0, O_ARG_INT, &clp.indel_loop, "-300", "score",
     "Score for insertions and deletions of loops per base"},
    {"indel-opening", 0, 0, O_ARG_INT, &clp.indel_opening, "-750", "score",
     clp.help_text["indel_opening"]},
    {"indel-opening-loop", 0, 0, O_ARG_INT, &clp.indel_opening_loop, "-900",
     "score", "Opening score for insertions and deletions of loops"},
    {"ribosum-file", 0, 0, O_ARG_STRING, &clp.ribosum_file, "RIBOSUM85_60", "f",
     clp.help_text["ribosum_file"]},
    {"use-ribosum", 0, 0, O_ARG_BOOL, &clp.use_ribosum, "true", "bool",
     clp.help_text["use_ribosum"]},
    {"match", 'm', 0, O_ARG_INT, &clp.match, "50", "score",
     clp.help_text["match"]},
    {"mismatch", 'M', 0, O_ARG_INT, &clp.mismatch, "0", "score",
     clp.help_text["mismatch"]},
    {"unpaired-penalty", 0, 0, O_ARG_INT, &clp.unpaired_penalty, "0", "score",
     clp.help_text["unpaired_penalty"]},
    {"struct-weight", 's', 0, O_ARG_INT, &clp.struct_weight, "200", "score",
     clp.help_text["struct_weight"]},
    {"exp-prob", 'e', &clp.exp_prob_given, O_ARG_DOUBLE, &clp.exp_prob,
     O_NODEFAULT, "prob", clp.help_text["exp_prob"]},
    {"tau", 't', 0, O_ARG_INT, &clp.tau, "100", "factor", clp.help_text["tau"]},
    {"exclusion", 'E', 0, O_ARG_INT, &clp.exclusion, "0", "score",
     clp.help_text["exclusion"]},
    {"stacking", 0, &clp.stacking, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["stacking"]},
    {"new-stacking", 0, &clp.new_stacking, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["new_stacking"]},

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Partition function representation (for sequence envelopes)"},

    {"extended-pf", 0, &clp.extended_pf, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["extended_pf_sequence_only"]+" [default]"},
    {"quad-pf", 0, &clp.quad_pf, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["quad_pf"]
#if !defined(_GLIBCXX_USE_FLOAT128) || defined(__clang__)
     +" Quad precision (128 bit, __float128) is not available for your binary. Falls back to extended-pf."
#endif
    },

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Controlling_output"},

    {"width", 'w', 0, O_ARG_INT, &clp.width, "120", "columns",
     clp.help_text["width"]},
    {"clustal", 0, &clp.clustal_given, O_ARG_STRING, &clp.clustal, O_NODEFAULT,
     "file", clp.help_text["clustal"]},
    {"stockholm", 0, &clp.stockholm_given, O_ARG_STRING, &clp.stockholm,
     O_NODEFAULT, "file", clp.help_text["stockholm"]},
    {"pp", 0, &clp.pp_given, O_ARG_STRING, &clp.pp, O_NODEFAULT, "file",
     clp.help_text["pp"]},
    {"alifold-consensus-dp", 0, &clp.alifold_consensus_dp, O_NO_ARG, 0,
     O_NODEFAULT, "", clp.help_text["alifold_consensus_dp"]},
    {"consensus-structure", 0, 0, O_ARG_STRING, &clp.cons_struct_type,
     "alifold", "type", clp.help_text["cons_struct_type"]},
    {"local-output", 'L', &clp.local_output, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["local_output"]},
    {"local-file-output", 0, &clp.local_file_output, O_NO_ARG, 0, O_NODEFAULT,
     "", clp.help_text["local_file_output"]},
    {"pos-output", 'P', &clp.pos_output, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["pos_output"]},
    {"write-structure", 0, &clp.write_structure, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["write_structure"]},
    {"special-gap-symbols", 0, &clp.special_gap_symbols, O_NO_ARG, 0,
     O_NODEFAULT, "",
     "Special distinct gap symbols for loop gaps or gaps caused by "
     "sparsification"},
    {"stopwatch", 0, &clp.stopwatch, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["stopwatch"]},

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "",
     "Heuristics for speed accuracy trade off"},

    {"min-prob", 'p', 0, O_ARG_DOUBLE, &clp.min_prob, "0.001", "prob",
     clp.help_text["min_prob"]},
    {"prob-unpaired-in-loop-threshold", 0, 0, O_ARG_DOUBLE,
     &clp.prob_unpaired_in_loop_threshold, "0.00005", "threshold",
     "Threshold for prob_unpaired_in_loop"},
    // todo: is the default threshold value reasonable?
    {"prob-basepair-in-loop-threshold", 0, 0, O_ARG_DOUBLE,
     &clp.prob_basepair_in_loop_threshold, "0.0001", "threshold",
     "Threshold for prob_basepair_in_loop"},
    {"max-bps-length-ratio", 0, 0, O_ARG_DOUBLE, &clp.max_bps_length_ratio,
     "0.0", "factor", clp.help_text["max_bps_length_ratio"]},
    {"max-uil-length-ratio", 0, 0, O_ARG_DOUBLE, &clp.max_uil_length_ratio,
     "0.0", "factor",
     "Maximal ratio of #unpaired bases in loops divided by sequence length "
     "(def: no effect)"},
    {"max-bpil-length-ratio", 0, 0, O_ARG_DOUBLE, &clp.max_bpil_length_ratio,
     "0.0", "factor",
     "Maximal ratio of #base pairs in loops divided by loop length (def: no "
     "effect)"},
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

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "MEA score"},

    {"mea-alignment", 0, &clp.mea_alignment, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["mea_alignment"]},
    {"match-prob-method", 0, 0, O_ARG_INT, &clp.match_prob_method, "0", "int",
     clp.help_text["match_prob_method"]},
    {"probcons-file", 0, &clp.probcons_file_given, O_ARG_STRING,
     &clp.probcons_file, O_NODEFAULT, "file", clp.help_text["probcons_file"]},
    {"temperature-alipf", 0, 0, O_ARG_INT, &clp.temperature_alipf, "300", "int",
     clp.help_text["temperature_alipf"]},
    {"pf-struct-weight", 0, 0, O_ARG_INT, &clp.pf_struct_weight, "200",
     "weight", clp.help_text["pf_struct_weight"]},
    {"mea-gapcost", 0, &clp.mea_gapcost, O_NO_ARG, 0, O_NODEFAULT, "",
     "Use gap cost in mea alignment"},
    {"mea-alpha", 0, 0, O_ARG_INT, &clp.mea_alpha, "0", "weight",
     clp.help_text["mea_alpha"]},
    {"mea-beta", 0, 0, O_ARG_INT, &clp.mea_beta, "200", "weight",
     clp.help_text["mea_beta"]},
    {"mea-gamma", 0, 0, O_ARG_INT, &clp.mea_gamma, "100", "weight",
     clp.help_text["mea_gamma"]},
    {"probability-scale", 0, 0, O_ARG_INT, &clp.probability_scale, "10000",
     "scale", clp.help_text["probability_scale"]},
    {"write-match-probs", 0, &clp.write_matchprobs, O_ARG_STRING,
     &clp.matchprobs_outfile, O_NODEFAULT, "file",
     clp.help_text["write_matchprobs"]},
    {"read-match-probs", 0, &clp.read_matchprobs, O_ARG_STRING,
     &clp.matchprobs_infile, O_NODEFAULT, "file",
     clp.help_text["read_matchprobs"]},
    {"write-arcmatch-scores", 0, &clp.write_arcmatch_scores, O_ARG_STRING,
     &clp.arcmatch_scores_outfile, O_NODEFAULT, "file",
     clp.help_text["write_arcmatch_scores"]},
    {"read-arcmatch-scores", 0, &clp.read_arcmatch_scores, O_ARG_STRING,
     &clp.arcmatch_scores_infile, O_NODEFAULT, "file",
     clp.help_text["read_arcmatch_scores"]},
    {"read-arcmatch-probs", 0, &clp.read_arcmatch_probs, O_ARG_STRING,
     &clp.arcmatch_scores_infile, O_NODEFAULT, "file",
     clp.help_text["read_arcmatch_probs"]},

    {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Constraints"},

    {"noLP", 0, &clp.no_lonely_pairs, O_NO_ARG, 0, O_NODEFAULT, "",
     clp.help_text["no_lonely_pairs"]},
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
 * \brief Main method of executable locarna
 *
 * @param argc argument counter
 * @param argv argument vector
 *
 * @return success
 */
int
main(int argc, char **argv) {
    stopwatch.start("total");

    // make sure that unsupported features are turned off; consider moving
    // to standard clp class
    clp.local_file_output = false;
    clp.pos_output = false;
    clp.local_output = false;
    clp.struct_local = false;
    clp.sequ_local = false;
    clp.free_endgaps = "";
    // clp.normalized=0;
    clp.stacking = false;
    clp.new_stacking = false;

    // ------------------------------------------------------------
    // Process options
    bool process_success = process_options(argc, argv, my_options);

    if (clp.help) {
        cout << "sparse - fast pairwise fast alignment of RNAs." << endl
             << endl;

        // cout << VERSION_STRING<<endl<<endl;

        print_help(argv[0], my_options);

        cout << "Report bugs to <miladim (at) informatik.uni-freiburg.de>."
             << endl
             << endl;
        return 0;
    }

    if (clp.quiet) {
        clp.verbose = false;
    } // quiet overrides verbose

    if (clp.galaxy_xml) {
        print_galaxy_xml((char *)"sparse", my_options);
        return 0;
    }

    if (clp.version || clp.verbose) {
        cout << "sparse (" << VERSION_STRING << ")" << endl;
        if (clp.version)
            return 0;
        else
            cout << endl;
    }

    if (!process_success) {
        std::cerr << "ERROR --- " << O_error_msg << std::endl;
        print_usage(argv[0], my_options);
        return -1;
    }

    if (clp.stopwatch) {
        stopwatch.set_print_on_exit(true);
    }

    if (clp.verbose) {
        print_options(my_options);
    }

    // --------------------
    // Forbid unsupported option of SPARSE
    if (clp.struct_local) {
        std::cerr << "Exclusions is not supported" << std::endl;
        return -1;
    }

    // noLP is not supported by sparse recursion but yet useful for calculating
    // probablities with RNAfold
    if (clp.no_lonely_pairs) {
        // std::cerr << "WARNING: No lonely pairs option is not supported by
        // sparse algortihm" << std::endl;
        //      return -1;
    }
    if (clp.sequ_local) {
        std::cerr << "Local sequence alignment is not supported" << std::endl;
        return -1;
    }
    if (clp.stacking || clp.new_stacking) {
        std::cerr << "Stacking is not supported" << std::endl;
        return -1;
    }
    /*  if(clp.free_endgaps.compare("----") != 0 ) {
        std::cerr << "Free end gaps is not supported" << std::endl;
        return -1;
        }
    */

    // ------------------------------------------------------------
    // parameter consistency
    if (clp.read_arcmatch_scores && clp.read_arcmatch_probs) {
        std::cerr << "You cannot specify arc match score and probabilities "
                     "file simultaneously."
                  << std::endl;
        return -1;
    }

    if (clp.probability_scale <= 0) {
        std::cerr << "Probability scale must be greater 0." << std::endl;
        return -1;
    }

    if (clp.struct_weight < 0) {
        std::cerr << "Structure weight must be greater equal 0." << std::endl;
        return -1;
    }

    // ----------------------------------------
    // temporarily turn off stacking unless background prob is set
    //
    if (clp.stacking && !clp.exp_prob_given) {
        std::cerr << "WARNING: stacking turned off. "
                  << "Stacking requires setting a background probability "
                  << "explicitely (option --exp-prob)." << std::endl;
        clp.stacking = false;
    }

    if (clp.ribofit) {
        clp.use_ribosum = false;
    }

    clp.extended_pf = 1; //default on

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

    // ----------------------------------------
    // Ribosum matrix
    //
    std::unique_ptr<RibosumFreq> ribosum;
    std::unique_ptr<Ribofit> ribofit;
    MainHelper::init_ribo_matrix(clp, ribosum, ribofit);

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    PFoldParams pfoldparams(PFoldParams::args::noLP(clp.no_lonely_pairs),
                            PFoldParams::args::stacking(clp.stacking || clp.new_stacking),
                            PFoldParams::args::max_bp_span(clp.max_bp_span));

    std::unique_ptr<ExtRnaData> rna_dataA;
    try {
        rna_dataA =
            std::make_unique<ExtRnaData>(clp.fileA, clp.min_prob,
                                         clp.prob_basepair_in_loop_threshold,
                                         clp.prob_unpaired_in_loop_threshold,
                                         clp.max_bps_length_ratio,
                                         clp.max_uil_length_ratio,
                                         clp.max_bpil_length_ratio, pfoldparams);
    } catch (failure &f) {
        std::cerr << "ERROR: failed to read from file " << clp.fileA
                  << std::endl
                  << "       " << f.what() << std::endl;
        return -1;
    }

    std::unique_ptr<ExtRnaData> rna_dataB = 0;
    try {
        rna_dataB =
            std::make_unique<ExtRnaData>(clp.fileB, clp.min_prob,
                                         clp.prob_basepair_in_loop_threshold,
                                         clp.prob_unpaired_in_loop_threshold,
                                         clp.max_bps_length_ratio,
                                         clp.max_uil_length_ratio,
                                         clp.max_bpil_length_ratio, pfoldparams);
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

    //---------------------------------------------------------------
    // Anchor constraint alignment is not supported by sparse yet
    if (!(seqA.annotation(MultipleAlignment::AnnoType::anchors)
              .single_string() == "") ||
        !(seqB.annotation(MultipleAlignment::AnnoType::anchors)
              .single_string() == "")) {
        std::cout << "WARNING sequence constraints found in the input but will "
                     "be ignored."
                  << std::endl;
    }
    LocARNA::SequenceAnnotation emptyAnnotation;
    rna_dataA->set_anchors(emptyAnnotation);
    rna_dataB->set_anchors(emptyAnnotation);
    //---------------------------------------------------------------

    // --------------------
    // handle max_diff restriction

    // missing: proper error handling in case that lenA, lenB, and
    // max_diff_pw_alignment/max_diff_alignment_file are incompatible

    // do inconsistency checking for max_diff_pw_alignment and
    // max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment != "" && clp.max_diff_alignment_file != "") {
        std::cerr
            << "Cannot simultaneously use both options --max-diff-pw-alignment"
            << " and --max-diff-alignment-file." << std::endl;
        return -1;
    }

    // construct TraceController and check inconsistency for with
    // multiplicity of sequences
    //

    std::unique_ptr<MultipleAlignment> multiple_ref_alignment;

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
                << "Invalid argument to --max-diff-pw-alignemnt;"
                << " require exactly one '&' separating the alignment strings."
                << std::endl;
            return -1;
        }

        if (alistr[0].length() != alistr[1].length()) {
            std::cerr << "Invalid argument to --max-diff-pw-alignemnt;"
                      << " alignment strings have unequal lengths."
                      << std::endl;
            return -1;
        }

        multiple_ref_alignment =
            std::make_unique<MultipleAlignment>(seqA.seqentry(0).name(),
                                                seqB.seqentry(0).name(),
                                                alistr[0], alistr[1]);
    }

    // ------------------------------------------------------------
    // Handle constraints (optionally)

    AnchorConstraints seq_constraints(
        lenA,
        seqA.annotation(MultipleAlignment::AnnoType::anchors).single_string(),
        lenB,
        seqB.annotation(MultipleAlignment::AnnoType::anchors).single_string(),
        !clp.relaxed_anchors);

    if (clp.verbose) {
        if (!seq_constraints.empty()) {
            std::cout << "Found sequence constraints." << std::endl;
        }
    }

    TraceController trace_controller(seqA, seqB, multiple_ref_alignment.get(),
                                     clp.max_diff, clp.max_diff_relax);

    trace_controller.restrict_by_anchors(seq_constraints);

    restrict_trace_by_probabilities(clp, rna_dataA.get(), rna_dataB.get(),
                                    ribosum.get(), ribofit.get(),
                                    &trace_controller,
                                    pf_score_t());

    // ----------------------------------------
    // construct set of relevant arc matches
    //
    std::unique_ptr<ArcMatches> arc_matches;

    // ------------------------------------------------------------
    // handle reading and writing of arcmatch_scores
    //
    // (needed for mea alignment with probabilistic consistency
    // transformation of arc match scores)
    //
    if (clp.read_arcmatch_scores || clp.read_arcmatch_probs) {
        if (clp.verbose) {
            std::cout << "Read arcmatch scores from file "
                      << clp.arcmatch_scores_infile << "." << std::endl;
        }
        arc_matches = std::make_unique<ArcMatches>(
            seqA, seqB, clp.arcmatch_scores_infile, clp.read_arcmatch_probs
                ? ((clp.mea_beta * clp.probability_scale) / 100)
                : -1,
            clp.max_diff_am != -1 ? (size_type)clp.max_diff_am
                                  : std::max(lenA, lenB),
            clp.max_diff_at_am != -1 ? (size_type)clp.max_diff_at_am
                                     : std::max(lenA, lenB),
            trace_controller, seq_constraints);
    } else {
        // initialize from RnaData
        arc_matches =
            std::make_unique<ArcMatches>(*rna_dataA.get(), *rna_dataB.get(),
                                         clp.min_prob, clp.max_diff_am != -1
                                             ? (size_type)clp.max_diff_am
                                             : std::max(lenA, lenB),
                                         clp.max_diff_at_am != -1
                                             ? (size_type)clp.max_diff_at_am
                                             : std::max(lenA, lenB),
                                         trace_controller, seq_constraints);
    }

    const BasePairs &bpsA = arc_matches->get_base_pairsA();
    const BasePairs &bpsB = arc_matches->get_base_pairsB();

    // ----------------------------------------
    // report on input in verbose mode
    if (clp.verbose)
        MainHelper::report_input(seqA, seqB, *arc_matches);

    // construct sparsification mapper for seqs A,B
    SparsificationMapper mapperA(bpsA, *rna_dataA.get(),
                                 clp.prob_unpaired_in_loop_threshold,
                                 clp.prob_basepair_in_loop_threshold, true);
    SparsificationMapper mapperB(bpsB, *rna_dataB.get(),
                                 clp.prob_unpaired_in_loop_threshold,
                                 clp.prob_basepair_in_loop_threshold, true);

    // ------------------------------------------------------------
    // Sequence match probabilities (for MEA-Alignment)
    //
    // perform parameter consistency checks
    if (clp.read_matchprobs && !clp.mea_alignment) {
        std::cerr
            << "Warning: clp.read_matchprobs ignored for non-mea alignment.\n";
    }
    if ((clp.write_matchprobs || clp.mea_alignment) && ribosum == NULL &&
        ribofit == NULL) {
        std::cerr << "ERROR: Ribosum/fit is required for mea_alignment"
                  << " and computing matchprobs." << std::endl;
        exit(-1);
    }
    //
    std::unique_ptr<MatchProbs> match_probs;
    if (clp.write_matchprobs || clp.mea_alignment) {
        match_probs =
            MainHelper::init_match_probs(clp, rna_dataA.get(), rna_dataB.get(),
                                         &trace_controller, ribosum.get(),
                                         ribofit.get(), pf_score_t());
    }
    if (clp.write_matchprobs) {
        MainHelper::write_match_probs(clp, match_probs.get());
        if (!clp.write_arcmatch_scores) {
            return 0;
        } // return from main()
    }
    //

    // ----------------------------------------
    // construct scoring

    // Scoring Parameter
    //
    double my_exp_probA = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenA);
    double my_exp_probB = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenB);
    //
    ScoringParams scoring_params(
        ScoringParams::match(clp.match),
        ScoringParams::mismatch(clp.mismatch),
        // In true mea alignment gaps are only scored
        // for computing base match probs.
        // Consequently, we set the indel and indel
        // opening cost to 0 for the case of mea
        // alignment!
        ScoringParams::indel(
            (clp.mea_alignment && !clp.mea_gapcost)
                ? 0
                : (clp.indel *
                   (clp.mea_gapcost ? clp.probability_scale / 100 : 1))),
        ScoringParams::indel_opening(
            (clp.mea_alignment && !clp.mea_gapcost)
                ? 0
                : (clp.indel_opening *
                   (clp.mea_gapcost ? clp.probability_scale / 100 : 1))),
        ScoringParams::ribosum(ribosum.get()),
        ScoringParams::ribofit(ribofit.get()),
        ScoringParams::unpaired_penalty(clp.unpaired_penalty),
        ScoringParams::struct_weight(clp.struct_weight),
        ScoringParams::tau_factor(clp.tau),
        ScoringParams::exclusion(clp.exclusion),
        ScoringParams::exp_probA(my_exp_probA),
        ScoringParams::exp_probB(my_exp_probB),
        ScoringParams::temperature_alipf(clp.temperature_alipf),
        ScoringParams::stacking(clp.stacking),
        ScoringParams::new_stacking(clp.new_stacking),
        ScoringParams::mea_scoring(clp.mea_alignment),
        ScoringParams::mea_alpha(clp.mea_alpha),
        ScoringParams::mea_beta(clp.mea_beta),
        ScoringParams::mea_gamma(clp.mea_gamma),
        ScoringParams::probability_scale(clp.probability_scale),

        // sparse specific scoring paramers:
        ScoringParams::indel_loop(
            (clp.mea_alignment && !clp.mea_gapcost)
                ? 0
                : (clp.indel_loop *
                   (clp.mea_gapcost ? clp.probability_scale / 100 : 1))),
        ScoringParams::indel_opening_loop(
            (clp.mea_alignment && !clp.mea_gapcost)
                ? 0
                : (clp.indel_opening_loop *
                   (clp.mea_gapcost ? clp.probability_scale / 100 : 1))));

    Scoring scoring(seqA, seqB, *rna_dataA.get(), *rna_dataB.get(),
                    *arc_matches.get(), match_probs.get(), scoring_params);

    if (clp.write_arcmatch_scores) {
        if (clp.verbose) {
            std::cout << "Write arcmatch scores to file "
                      << clp.arcmatch_scores_outfile << " and exit."
                      << std::endl;
        }
        arc_matches->write_arcmatch_scores(clp.arcmatch_scores_outfile,
                                           scoring);
        return 0;
    }

    // ------------------------------------------------------------
    // Computation of the alignment score
    //

    // initialize aligner object, which does the alignment computation
    auto aligner = std::make_unique<AlignerN>(
        AlignerNParams(AlignerParams::seqA(&seqA), AlignerParams::seqB(&seqB),
                       AlignerParams::scoring(&scoring),
                       AlignerParams::struct_local(clp.struct_local),
                       AlignerParams::sequ_local(clp.sequ_local),
                       AlignerParams::free_endgaps(FreeEndgaps(clp.free_endgaps)),
                       AlignerParams::max_diff_am(clp.max_diff_am),
                       AlignerParams::max_diff_at_am(clp.max_diff_at_am),
                       AlignerParams::trace_controller(&trace_controller),
                       AlignerParams::stacking(clp.stacking ||
                                               clp.new_stacking),
                       AlignerParams::constraints(&seq_constraints),
                       AlignerNParams::sparsification_mapperA(&mapperA),
                       AlignerNParams::sparsification_mapperB(&mapperB)));

    infty_score_t score;

    // otherwise compute the best alignment
    score = aligner->align();

    // ----------------------------------------
    // report score
    //
    if (!clp.quiet) {
        std::cout << "Score: " << score << std::endl << std::endl;
    }

    std::unique_ptr<Alignment> alignment;
    bool return_code = 0;

    // ------------------------------------------------------------
    // Traceback
    //
    if (DO_TRACE) {
        aligner->trace();

        alignment = std::make_unique<Alignment>(aligner->get_alignment());

        aligner.reset(); // delete the aligner with DP matrices
        arc_matches.reset(); // delete the arc_matches

        // ----------------------------------------
        // write alignment in different output formats
        //
        std::string consensus_structure = "";

        std::unique_ptr<RnaData> consensus =
            MainHelper::consensus(clp, pfoldparams, my_exp_probA, my_exp_probB,
                                  rna_dataA.get(), rna_dataB.get(), *alignment,
                                  consensus_structure);

        return_code =
            MainHelper::write_alignment(clp, score, consensus_structure,
                                        consensus.get(), *alignment,
                                        multiple_ref_alignment.get());

        // ----------------------------------------
        // write alignment to screen

        if (!clp.quiet) {
            MultipleAlignment ma(*alignment, clp.local_output,
                                 clp.special_gap_symbols);

            if (clp.write_structure) {
                // annotate multiple alignment with structures
                std::string structureA =
                    alignment->dot_bracket_structureA(clp.local_output);
                std::string structureB =
                    alignment->dot_bracket_structureB(clp.local_output);
                ma.prepend(MultipleAlignment::SeqEntry("", structureA));
                ma.append(MultipleAlignment::SeqEntry("", structureB));
            }

            if (consensus_structure != "") {
                ma.append(MultipleAlignment::SeqEntry(clp.cons_struct_type,
                                                      consensus_structure));
            }

            ma.write(std::cout, clp.width,
                     MultipleAlignment::FormatType::CLUSTAL);

            std::cout << endl;
        }
    }

    stopwatch.stop("total");

    // ----------------------------------------
    // DONE
    return return_code;
}
