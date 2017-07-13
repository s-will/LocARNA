/**
 * \file locarna.cc
 *
 * \brief Defines main function of locarna
 *
 * LocARNA: global and LOCal Alignment of RNA
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de>
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <memory>

//#include <math.h>

#include "LocARNA/sequence.hh"
#include "LocARNA/scoring.hh"
#include "LocARNA/basepairs.hh"
#include "LocARNA/alignment.hh"
#include "LocARNA/aligner.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/match_probs.hh"
#include "LocARNA/ribosum.hh"
#include "LocARNA/ribofit.hh"
#include "LocARNA/anchor_constraints.hh"
#include "LocARNA/trace_controller.hh"
#include "LocARNA/global_stopwatch.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/rna_ensemble.hh"
#include "LocARNA/main_helper.icc"

// using namespace std;
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
//! @note never made it into command line (since no one asked me kindly)
const bool DO_TRACE = true;

//! \brief Structure for command line parameters of locarna
//!
//! Encapsulating all command line parameters in a common structure
//! avoids name conflicts and makes downstream code more informative.
//!
struct command_line_parameters
    : public MainHelper::std_command_line_parameters,
      public MainHelper::mea_command_line_parameters {
    //! find suboptimal solution (either k-best or all solutions
    //! better than a threshold)
    bool subopt;

    int kbest_k;          //!< kbest_k
    int subopt_threshold; //!< subopt_threshold

    bool normalized; //!< whether to do normalized alignment

    bool penalized;       //!< whether to penalize alignment positions
    int position_penalty; //!< position-wise penalty

    int normalized_L; //!< normalized_L

    bool score_components; //!< whether to report score components

    command_line_parameters()
        : MainHelper::std_command_line_parameters(),
          MainHelper::mea_command_line_parameters(help_text) {
        help_text["normalized"] =
            "Perform normalized local alignment with parameter L. "
            "This causes locarna to compute the best local alignment according "
            "to "
            "'Score' / ( L + 'length' ), "
            "where length is the sum of the lengths of the two locally aligned "
            "subsequences. "
            "Thus, the larger L, the larger the local alignment; the size of "
            "value L is in the order "
            "of local alignment lengths. Verbose yields info on the iterative "
            "optimizations.";
        help_text["penalized"] = "Penalized local alignment with penalty PP";
        help_text["normalized_L"] =
            "Parameter L for normalized local alignment. "
            "Larger values produce larger alignments.";
        help_text["score_components"] =
            "Output components of the score (experimental).";
    }
};

//! \brief holds command line parameters of locarna
command_line_parameters clp;

//! defines command line parameters
option_def my_options[] =
    {{"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "cmd_only"},

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

     {"indel", 'i', 0, O_ARG_INT, &clp.indel, "-350", "score",
      clp.help_text["indel"]},
     {"indel-opening", 0, 0, O_ARG_INT, &clp.indel_opening, "-500", "score",
      clp.help_text["indel_opening"]},
     {"ribosum-file", 0, 0, O_ARG_STRING, &clp.ribosum_file, "RIBOSUM85_60",
      "f", clp.help_text["ribosum_file"]},
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
     {"tau", 't', 0, O_ARG_INT, &clp.tau, "0", "factor", clp.help_text["tau"]},
     {"exclusion", 'E', 0, O_ARG_INT, &clp.exclusion, "0", "score",
      clp.help_text["exclusion"]},
     {"stacking", 0, &clp.stacking, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["stacking"]},
     {"new-stacking", 0, &clp.new_stacking, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["new_stacking"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Locality"},

     {"struct-local", 0, &clp.struct_local_given, O_ARG_BOOL, &clp.struct_local,
      "false", "bool", clp.help_text["struct_local"]},
     {"sequ-local", 0, &clp.sequ_local_given, O_ARG_BOOL, &clp.sequ_local,
      "false", "bool", clp.help_text["sequ_local"]},
     {"free-endgaps", 0, 0, O_ARG_STRING, &clp.free_endgaps, "----", "spec",
      clp.help_text["free_endgaps"]},
     {"normalized", 0, &clp.normalized, O_ARG_INT, &clp.normalized_L, "0", "L",
      clp.help_text["normalized"]},

     {"penalized", 0, &clp.penalized, O_ARG_INT, &clp.position_penalty, "0",
      "PP", clp.help_text["penalized"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Output"},

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
     {"score-components", 0, &clp.score_components, O_NO_ARG, 0, O_NODEFAULT,
      "", clp.help_text["score_components"]},
     {"stopwatch", 0, &clp.stopwatch, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["stopwatch"]},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "",
      "Heuristics for speed accuracy trade off"},

     {"min-prob", 'p', 0, O_ARG_DOUBLE, &clp.min_prob, "0.0005", "prob",
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

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "Special sauce options"},
     {"kbest", 0, &clp.subopt, O_ARG_INT, &clp.kbest_k, "-1", "k",
      "Enumerate k-best alignments"},
     {"better", 0, &clp.subopt, O_ARG_INT, &clp.subopt_threshold, "-1000000",
      "t", "Enumerate alignments better threshold t"},

     {"", 0, 0, O_SECTION, 0, O_NODEFAULT, "", "MEA score"},

     {"mea-alignment", 0, &clp.mea_alignment, O_NO_ARG, 0, O_NODEFAULT, "",
      clp.help_text["mea_alignment"]},
     {"match-prob-method", 0, 0, O_ARG_INT, &clp.match_prob_method, "0", "int",
      clp.help_text["match_prob_method"]},
     {"probcons-file", 0, &clp.probcons_file_given, O_ARG_STRING,
      &clp.probcons_file, O_NODEFAULT, "file", clp.help_text["probcons_file"]},
     {"temperature-alipf", 0, 0, O_ARG_INT, &clp.temperature_alipf, "150",
      "int", clp.help_text["temperature_alipf"]},
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

    typedef std::vector<int>::size_type size_type;

    // ------------------------------------------------------------
    // Process options

    bool process_success = process_options(argc, argv, my_options);

    if (clp.help) {
        std::cout << "locarna - pairwise (global and local) alignment of RNA."
                  << std::endl
                  << std::endl;

        // std::cout << VERSION_STRING<<std::endl<<std::endl;

        print_help(argv[0], my_options);

        std::cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."
                  << std::endl
                  << std::endl;
        return 0;
    }

    if (clp.quiet) {
        clp.verbose = false;
    } // quiet overrides verbose

    if (clp.galaxy_xml) {
        print_galaxy_xml((char *)"locarna", my_options);
        return 0;
    }

    if (clp.version || clp.verbose) {
        std::cout << VERSION_STRING << std::endl;
        if (clp.version)
            return 0;
        else
            std::cout << std::endl;
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

    // ------------------------------------------------------------
    // parameter consistency
    if (clp.read_arcmatch_scores && clp.read_arcmatch_probs) {
        std::cerr << "One cannot specify files for arc match scores "
                  << "and probabilities simultaneously." << std::endl;
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

    if (clp.normalized && clp.penalized) {
        std::cerr << "One cannot specify penalized and normalized "
                  << "simultaneously." << std::endl;
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

    // ------------------------------------------------------------
    // if normalized or penalized alignment shall be computed,
    // automatically turn on sequ_local unless sequ_local mode was
    // explicitly specified
    //
    // important: in the Aligner class, we rely on sequ_local==true in
    // normalized alignment mode
    if (clp.normalized || clp.penalized) {
        if (!clp.sequ_local_given) {
            clp.sequ_local = true;
        } else {
            if (!clp.sequ_local) {
                std::cerr << "ERROR: Cannot run normalized alignment "
                          << "without --sequ_local on." << std::endl;
                return -1;
            }
        }

        if (clp.struct_local_given && clp.struct_local) {
            std::cerr << "ERROR: Normalized structure local alignment "
                      << "not supported." << std::endl;
            return -1;
        } else {
            clp.struct_local = false;
        }
    }

    // ----------------------------------------
    // Ribosum matrix
    //
    std::unique_ptr<RibosumFreq> ribosum;
    std::unique_ptr<Ribofit> ribofit;

    MainHelper::init_ribo_matrix(clp, ribosum, ribofit);

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    PFoldParams pfparams(clp.no_lonely_pairs, clp.stacking || clp.new_stacking,
                         clp.max_bp_span, 2);

    std::unique_ptr<RnaData> rna_dataA;
    try {
        rna_dataA =
            std::make_unique<RnaData>(clp.fileA, clp.min_prob,
                                      clp.max_bps_length_ratio, pfparams);
    } catch (failure &f) {
        std::cerr << "ERROR:\tfailed to read from file " << clp.fileA
                  << std::endl
                  << "\t" << f.what() << std::endl;
        return -1;
    }

    std::unique_ptr<RnaData> rna_dataB = 0;
    try {
        rna_dataB =
            std::make_unique<RnaData>(clp.fileB, clp.min_prob,
                                      clp.max_bps_length_ratio, pfparams);
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
    // max_diff_pw_alignment/max_diff_alignment_file are incompatible

    // do inconsistency checking for max_diff_pw_alignment and
    // max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment != "" && clp.max_diff_alignment_file != "") {
        std::cerr << "Cannot simultaneously use options --max-diff-pw-alignment"
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
            std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning "
                      << "of alignments." << std::endl;
            return -1;
        }

        std::vector<std::string> alistr;
        split_at_separator(clp.max_diff_pw_alignment, '&', alistr);

        if (alistr.size() != 2) {
            std::cerr << "Invalid argument to --max-diff-pw-alignemnt; require "
                      << "exactly one '&' separating the alignment strings."
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

    TraceController trace_controller(seqA, seqB, multiple_ref_alignment.get(),
                                     clp.max_diff, clp.max_diff_relax);

    // ------------------------------------------------------------
    // Handle constraints (optionally)

    AnchorConstraints seq_constraints(
        lenA,
        seqA.annotation(MultipleAlignment::AnnoType::anchors).single_string(),
        lenB,
        seqB.annotation(MultipleAlignment::AnnoType::anchors).single_string(),
        !clp.relaxed_anchors);

    // not used, since there seems to be no effect
    // trace_controller.restrict_by_anchors(seq_constraints);

    if (clp.verbose) {
        if (!seq_constraints.empty()) {
            std::cout << "Found sequence constraints." << std::endl;
        }
    }

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
            std::make_unique<ArcMatches>(*rna_dataA, *rna_dataB, clp.min_prob,
                                         clp.max_diff_am != -1
                                             ? (size_type)clp.max_diff_am
                                             : std::max(lenA, lenB),
                                         clp.max_diff_at_am != -1
                                             ? (size_type)clp.max_diff_at_am
                                             : std::max(lenA, lenB),
                                         trace_controller, seq_constraints);
    }

    // note: arcmatches has to be assigned ( unless new failed, which
    // we don't catch )

    // ----------------------------------------
    // report on input in verbose mode
    if (clp.verbose)
        MainHelper::report_input(seqA, seqB, *arc_matches.get());

    // ------------------------------------------------------------
    // Sequence match probabilities (for MEA-Alignment)
    //
    // perform parameter consistency checks
    if (clp.read_matchprobs && !clp.mea_alignment) {
        std::cerr << "Warning: clp.read_matchprobs ignored for "
                  << "non-mea alignment.\n";
    }
    if ((clp.write_matchprobs || clp.mea_alignment) && ribosum == nullptr &&
        ribofit == nullptr) {
        std::cerr << "ERROR: Ribosum/fit is required for mea_alignment and "
                  << "computing matchprobs." << std::endl;
        exit(-1);
    }
    //
    std::unique_ptr<MatchProbs> match_probs;

    if (clp.write_matchprobs || clp.mea_alignment) {
        match_probs =
            MainHelper::init_match_probs(clp, rna_dataA.get(), rna_dataB.get(),
                                         ribosum.get(), ribofit.get());
    }
    if (clp.write_matchprobs) {
        MainHelper::write_match_probs(clp, match_probs.get());
        if (!clp.write_arcmatch_scores) {
            return 0;
        } // return from main()
    }
    //

    // ----------------------------------------
    // Scoring Parameter
    //
    double my_exp_probA = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenA);
    double my_exp_probB = clp.exp_prob_given ? clp.exp_prob : prob_exp_f(lenB);

    ScoringParams scoring_params(
        clp.match, clp.mismatch,
        // In true mea alignment gaps are only scored
        // for computing base match probs.
        // Consequently, we set the indel and indel
        // opening cost to 0 for the case of mea
        // alignment!
        (clp.mea_alignment && !clp.mea_gapcost)
            ? 0
            : (clp.indel * (clp.mea_gapcost ? clp.probability_scale / 100 : 1)),
        0, // indel__loop_score, for
        // consistency and least modification to
        // locarna.cc has been set to zero
        (clp.mea_alignment && !clp.mea_gapcost)
            ? 0
            : (clp.indel_opening *
               (clp.mea_gapcost ? clp.probability_scale / 100 : 1)),
        0, // indel_opening_loop_score, for
        // consistency and least modification to
        // locarna.cc has been set to zero
        ribosum.get(), ribofit.get(), clp.unpaired_penalty, clp.struct_weight,
        clp.tau, clp.exclusion, my_exp_probA, my_exp_probB,
        clp.temperature_alipf, clp.stacking, clp.new_stacking,
        clp.mea_alignment, clp.mea_alpha, clp.mea_beta, clp.mea_gamma,
        clp.probability_scale);
    // ------------------------------------------------------------
    // Construct scoring

    Scoring scoring(seqA, seqB, *rna_dataA, *rna_dataB, *arc_matches,
                    match_probs.get(), scoring_params);

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
    auto aligner =
        std::make_unique<Aligner>(Aligner::create()
                                  .seqA(seqA)
                                  .seqB(seqB)
                                  .scoring(scoring)
                                  .no_lonely_pairs(clp.no_lonely_pairs)
                                  .struct_local(clp.struct_local)
                                  .sequ_local(clp.sequ_local)
                                  .free_endgaps(clp.free_endgaps)
                                  .max_diff_am(clp.max_diff_am)
                                  .max_diff_at_am(clp.max_diff_at_am)
                                  .trace_controller(trace_controller)
                                  .stacking(clp.stacking || clp.new_stacking)
                                  .constraints(seq_constraints));

    // enumerate suboptimal alignments (using interval splitting)
    if (clp.subopt) {
        aligner->suboptimal(clp.kbest_k, clp.subopt_threshold, clp.normalized,
                           clp.normalized_L, clp.width, clp.verbose,
                           clp.local_output, clp.pos_output,
                           clp.write_structure);
        return 0;
    }

    infty_score_t score;

    // if option --normalized <L> is given, then do normalized local alignemnt
    if (clp.normalized) {
        score = aligner->normalized_align(clp.normalized_L, clp.verbose);

    } else if (clp.penalized) {
        score = aligner->penalized_align(clp.position_penalty);
    }

    else {
        // ========== STANDARD CASE ==========

        // otherwise compute the best alignment
        score = aligner->align();
    }

    // ----------------------------------------
    // report score
    //
    if (!clp.quiet) {
        std::cout << "Score: " << score << std::endl << std::endl;
    }

    std::unique_ptr<Alignment> alignment;

    // ------------------------------------------------------------
    // Traceback
    //
    if (clp.normalized || clp.penalized || DO_TRACE) {
        if (!(clp.normalized || clp.penalized)) {
            aligner->trace();
        } //otherwise the trace was already computed!

        // copy the alignment from aligner
        alignment = std::make_unique<Alignment>(aligner->get_alignment());
    }

    aligner.reset(); // delete the aligner with DP matrices
    arc_matches.reset(); // delete the arc_matches

    if ((!clp.normalized && !clp.penalized) && DO_TRACE) {

        // if score components should be reported
        if (clp.score_components) {
            // get alignment and compute contributions by sequence
            // similarity and gap cost

            // experimental: I do the calculations here until I
            // decide where to put the code

            std::cout << "WARNING: reporting score components is still "
                      << "experimental." << std::endl;
            std::cout << "         This does not work properly for certain "
                      << "non-standard standard scoring types, "
                      << "e.g. free end gaps." << std::endl;

            const Alignment::edges_t edges = alignment->alignment_edges(false);

            score_t seq_sim = 0;
            score_t gap_cost = 0; // count linear component of gap cost
            score_t gap_numA = 0; // count number of gap openings in A
            score_t gap_numB = 0; // count number of gap openings in B

            bool openA = false; // is a gap open in A
            bool openB = false; // is a gap open in B

            for ( const auto & edge : edges ) {
                if (edge.first.is_pos() && edge.second.is_pos()) {
                    seq_sim +=
                        scoring.basematch(edge.first, edge.second);
                }
                if (edge.first.is_gap()) {
                    if (!openA) {
                        gap_numA++;
                        openA = true;
                    }
                    gap_cost += scoring.gapA(edge.second);
                } else {
                    openA = false;
                }

                if (edge.second.is_gap()) {
                    if (!openB) {
                        gap_numB++;
                        openB = true;
                    }
                    gap_cost += scoring.gapB(edge.first);
                } else {
                    openB = false;
                }
            }

            score_t total_gap_cost =
                gap_cost + (gap_numA + gap_numB) * scoring.indel_opening();

            std::cout << "#SCORE total        " << std::setw(8) << score
                      << std::endl;
            std::cout << "#SCORE seq_sim      " << std::setw(8) << seq_sim
                      << std::endl;
            std::cout << "#SCORE gap_penalty  " << std::setw(8)
                      << total_gap_cost << std::endl;
            std::cout << "#SCORE str_contrib  " << std::setw(8)
                      << score - seq_sim - total_gap_cost << std::endl;
        }
    }

    int return_code = 0;

    if (clp.normalized || clp.penalized || DO_TRACE) {
        // if we did a trace (one way or the other)

        // ----------------------------------------
        // write alignment in different output formats to files
        //

        std::string consensus_structure = "";

        std::unique_ptr<RnaData> consensus =
            MainHelper::consensus(clp, pfparams, my_exp_probA, my_exp_probB,
                                  rna_dataA.get(), rna_dataB.get(), *alignment,
                                  consensus_structure);

        return_code =
            MainHelper::write_alignment(clp, score, consensus_structure,
                                        consensus.get(), *alignment,
                                        multiple_ref_alignment.get());

        // ----------------------------------------
        // write alignment to screen

        if (clp.pos_output) {
            std::cout << "HIT " << score << " " << alignment->local_start().first
                      << " " << alignment->local_start().second << " "
                      << alignment->local_end().first << " " << alignment->local_end().second
                      << " " << std::endl;
            std::cout << std::endl;
        }

        if ((!clp.pos_output && !clp.quiet) || clp.local_output) {
            MultipleAlignment ma(*alignment, clp.local_output);

            if (clp.write_structure) {
                // annotate multiple alignment with structures
                std::string structureA =
                    alignment->dot_bracket_structureA(clp.local_output);
                std::string structureB =
                    alignment->dot_bracket_structureB(clp.local_output);
                ma.prepend(MultipleAlignment::SeqEntry("", structureA));
                ma.append(MultipleAlignment::SeqEntry("", structureB));
            }

            if (clp.pos_output) {
                std::cout << "\t+" << alignment->local_start().first << std::endl
                          << "\t+" << alignment->local_start().second << std::endl
                          << std::endl;
                std::cout << std::endl;
            }

            // recompute consensus structure for screen output, if necessary
            if (consensus_structure != "") {
                if (clp.local_output != clp.local_file_output) {
                    clp.local_file_output = clp.local_output;

                    consensus =
                        MainHelper::consensus(clp, pfparams, my_exp_probA,
                                              my_exp_probB, rna_dataA.get(),
                                              rna_dataB.get(), *alignment,
                                              consensus_structure);
                    clp.local_file_output = !clp.local_output;
                }
                ma.append(MultipleAlignment::SeqEntry(clp.cons_struct_type,
                                                      consensus_structure));
            }
            ma.write(std::cout, clp.width,
                     MultipleAlignment::FormatType::CLUSTAL);

            if (clp.pos_output) {
                std::cout << std::endl
                          << "\t+" << alignment->local_end().second << std::endl
                          << "\t+" << alignment->local_end().second << std::endl
                          << std::endl;
            }
        }

        if (!clp.quiet) {
            std::cout << std::endl;
        }
    }

    stopwatch.stop("total");

    // ----------------------------------------
    // DONE
    return return_code;
}
