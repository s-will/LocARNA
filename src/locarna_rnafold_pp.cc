/************************************************************
 *
 * \file locarna_rnafold_pp.cc
 * \brief Compute and write base pair probabilities in pp-format.
 *
 * Special RNAfold-like program that computes a probability dot plot
 * for a input RNA sequence and writes a file in the
 * LocARNA-internal pp-format.
 *
 * locarna_rnafold_pp folds a sequence or multiple alignment using
 * partition function folding and writes the result in pp 2.0 format.
 *
 * This program is part of the LocARNA package. It is intended for
 * computing pair probabilities of the input sequences.
 *
 ************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <math.h>

#include <string.h>
#include <sstream>
#include <string>

#include <LocARNA/options.hh>
#include <LocARNA/multiple_alignment.hh>
#include <LocARNA/pfold_params.hh>
#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/rna_data.hh>
#include <LocARNA/ext_rna_data.hh>

using namespace LocARNA;

/**
 * \brief Structure for command line parameters
 *
 * Encapsulating all command line parameters in a common structure
 * avoids name conflicts and makes downstream code more informative.
 *
 */
struct command_line_parameters {
    bool help;                   //!< whether to print help
    bool version;                //!< whether to print version
    bool verbose;                //!< whether to print verbose output
    std::string input_file;      //!< input_file
    bool use_struct_constraints; //!< -C use structural constraints
    bool no_lonely_pairs;        //!< no lonely pairs option
    int max_bp_span;             //!< maximum base pair span
    bool stacking;               //!< whether to stacking
    int dangling;                //!< dangling option value
    bool in_loop;                //!< whether to compute in-loop probabilities
    double min_prob;             //!< minimal / cutoff base pair probability
    double prob_unpaired_in_loop_threshold; //!< threshold for
                                            //!prob_unpaired_in_loop
    double prob_basepair_in_loop_threshold; //!< threshold for
                                            //!prob_basepait_in_loop
    std::string output_file;                //!< output file name
    bool force_alifold; //!< use alifold even for single sequences.
};
//! \brief holds command line parameters of locarna
command_line_parameters clp;
// longname,shortname,flag,arg_type,argument,default,argname,description
//! defines command line parameters
option_def my_options[] =
    {{"help", 'h', &clp.help, O_NO_ARG, 0, O_NODEFAULT, "", "Help"},
     {"version", 'V', &clp.version, O_NO_ARG, 0, O_NODEFAULT, "",
      "Version info"},
     {"verbose", 'v', &clp.verbose, O_NO_ARG, 0, O_NODEFAULT, "", "Verbose"},
     {"use-struct-constraints", 'C', &clp.use_struct_constraints, O_NO_ARG, 0,
      O_NODEFAULT, "", "Use structural constraints"},
     {"noLP", 0, &clp.no_lonely_pairs, O_NO_ARG, 0, O_NODEFAULT, "",
      "No lonely pairs"},
     {"maxBPspan", 0, 0, O_ARG_INT, &clp.max_bp_span, "-1", "span",
      "Limit maximum base pair span (default=off)"},
     {"stacking", 0, &clp.stacking, O_NO_ARG, 0, O_NODEFAULT, "",
      "Compute stacking terms"},
     {"dangling", 0, 0, O_ARG_INT, &clp.dangling, "2", "",
      "Dangling option value"},
     {"in-loop", 0, &clp.in_loop, O_NO_ARG, 0, O_NODEFAULT, "",
      "Compute in-loop probabilities"},
     {"min-prob", 'p', 0, O_ARG_DOUBLE, &clp.min_prob, "0.0005", "prob",
      "Minimal probability"},
     {"p_unpaired_in_loop", 0, 0, O_ARG_DOUBLE,
      &clp.prob_unpaired_in_loop_threshold, "0.0005", "threshold",
      "Threshold for prob_unpaired_in_loop"},
     {"p_basepair_in_loop", 0, 0, O_ARG_DOUBLE,
      &clp.prob_basepair_in_loop_threshold, "0.0005", "threshold",
      "Threshold for prob_basepair_in_loop"}, // todo: is the default threshold
                                              // value reasonable?
     {"output", 'o', 0, O_ARG_STRING, &clp.output_file, "", "filename",
      "Output file"},
     {"force-alifold", 0, &clp.force_alifold, O_NO_ARG, 0, O_NODEFAULT, "",
      "Force alifold for single seqeunces"},
     {"", 0, 0, O_ARG_STRING, &clp.input_file, "-", "filename", "Input file"},
     {"", 0, 0, 0, 0, O_NODEFAULT, "", ""}};

/**
 * \brief Main function of locarna_rnafold_pp when Vienna RNA lib is linked
 */
int
main(int argc, char **argv) {
    // ------------------------------------------------------------
    // Process options

    bool process_success = process_options(argc, argv, my_options);

    if (clp.help) {
        std::cout << "locarna_rnafold_pp -- compute RNA pair probabilities and "
                     "write in pp-format"
                  << std::endl;
        //      std::cout << VERSION_STRING<<std::endl<<std::endl;
        print_help(argv[0], my_options);
        return 0;
    }

    if (clp.version || clp.verbose) {
        std::cout << "locarna_rnafold_pp (" << PACKAGE_STRING << ")"
                  << std::endl;
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

    // Reading from stdinput with autodetect of file format works by copying the
    // entire stdinput
    // to memory. Then, autodetection can work on this copy.

    // if we want to get input from stdin, copy content of stdin to
    // stdin_content
    std::string stdin_content;
    if (clp.input_file == "-") {
        std::stringstream stdin;
        stdin << std::cin.rdbuf();
        stdin_content = stdin.str();
    }

    // todo: check that input is proper and catch the wrong inputs
    // MultipleAlignment::FormatType::type input_format;
    MultipleAlignment *mseq = NULL;
    // try fasta format
    bool failed = false;
    try {
        if (clp.input_file.compare("-") != 0) {
            mseq = new MultipleAlignment(clp.input_file,
                                         MultipleAlignment::FormatType::FASTA);
        } else {
            std::istringstream in(stdin_content);
            mseq =
                new MultipleAlignment(in, MultipleAlignment::FormatType::FASTA);
        }
        // even if reading does not fail, we still want to
        // make sure that the result is reasonable. Otherwise,
        // we assume that the file is in a different format.
        if (!mseq->is_proper() || mseq->empty()) {
            failed = true;
        }
    } catch (failure &f) {
        failed = true;
        //      std::cerr << "Not a FASTA file" << std::endl;
        //      std::cerr << f.what() << std::endl;
    }

    if (!failed) {
        // input_format = MultipleAlignment::FormatType::FASTA;
    } else { //
        failed = false;
        try {
            try {
                // std::cerr << "Try reading clustal "<<filename<<"
                // ..."<<std::endl;
                if (clp.input_file.compare("-") != 0) {
                    mseq = new MultipleAlignment(
                        clp.input_file, MultipleAlignment::FormatType::CLUSTAL);
                } else {
                    std::stringstream in(stdin_content);
                    mseq = new MultipleAlignment(
                        in, MultipleAlignment::FormatType::CLUSTAL);

                    // make sure that the result is reasonable. Otherwise,
                }
                // even if reading does not fail, we still want to
                // we assume that the file is in a different format.
                if (!mseq->is_proper() || mseq->empty()) {
                    failed = true;
                }
            } catch (syntax_error_failure &f) {
                throw failure((std::string)"Cannot read input data from clustal file.\n\t"+f.what());
            }
        } catch (failure &f) {
            failed = true;
            //          std::cerr << f.what() << std::endl;
        }
        if (failed) {
            std::cerr << "Error in input format" << std::endl;
            return -1;
        }
    }

    bool use_alifold = true;

    // if the input has only one sequence and forced to do alifold
    if (!clp.force_alifold && mseq->num_of_rows() == 1)
        use_alifold = false;
    if (mseq->has_annotation(MultipleAlignment::AnnoType::structure) &&
        !clp.use_struct_constraints) {
        std::cerr << "Warning locarna_rnafold_pp: structure constraints will "
                     "be ignored"
                  << std::endl;
        mseq->set_annotation(MultipleAlignment::AnnoType::structure,
                             SequenceAnnotation());
    }

    PFoldParams pfoldparams(PFoldParams::args::noLP(clp.no_lonely_pairs),
                            PFoldParams::args::stacking(clp.stacking),
                            PFoldParams::args::max_bp_span(clp.max_bp_span),
                            PFoldParams::args::dangling(clp.dangling));

    RnaEnsemble rna_ensemble(*mseq, pfoldparams, clp.in_loop, use_alifold);

    // write pp file

    // set the appropriate ostream from input_file or std::cout
    std::ofstream of;
    std::streambuf *buff;
    if (clp.output_file.length() == 0) {
        buff = std::cout.rdbuf();
    } else {
        of.open(clp.output_file.c_str());
        buff = of.rdbuf();
    }
    std::ostream out_stream(buff);

    if (clp.in_loop) {
        ExtRnaData ext_rna_data(
            rna_ensemble, clp.min_prob, clp.prob_basepair_in_loop_threshold,
            clp.prob_unpaired_in_loop_threshold,
            0, // don't filter output by max_bps_length_ratio
            0, // don't filter output by max_uil_length_ratio
            0, // don't filter output by max_bpil_length_ratio
            pfoldparams);

        ext_rna_data.write_pp(out_stream); // (no need to filter again => don't
                                           // specify output cutoff)
    } else {
        RnaData rna_data(rna_ensemble, clp.min_prob,
                         0, // don't filter output by max_bps_length_ratio
                         pfoldparams);

        rna_data.write_pp(out_stream); // (no need to filter again => don't
                                       // specify output cutoff)
    }

    return 0;
}
