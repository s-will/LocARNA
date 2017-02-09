/**********************************************************************
 *
 * Exparna-P: Exact Pattern Matching for RNA Structure Ensembles
 * fast structure local exact matching
 *
 * Copyright (C) Christina Otto <schmiedc(@)informatik.uni-freiburg.de>
 *
 **********************************************************************/

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include<limits>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/types.h>
// for gettimeofday()
#include <sys/time.h>
// for setprecision
#include <iomanip>

#include "LocARNA/sequence.hh"
#include "LocARNA/basepairs.hh"

#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/match_probs.hh"

#include "LocARNA/anchor_constraints.hh"
#include "LocARNA/sequence_annotation.hh"
#include "LocARNA/trace_controller.hh"

#include "LocARNA/exact_matcher.hh"
#include "LocARNA/sparsification_mapper.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/global_stopwatch.hh"

#include "LocARNA/main_helper.icc"

using namespace std;

const std::string
VERSION_STRING = (std::string)PACKAGE_STRING;

const bool DO_TRACE=true;
//const bool DO_TRACE=false;

// ------------------------------------------------------------
// Parameter

struct command_line_parameters {
    bool help;
    bool version;
    bool verbose;
    bool quiet;
    bool postscript_output;
    bool subopt;
    bool add_filter;
    bool no_stacking;

    bool stopwatch;

    double min_prob; // only pairs with a probability of at least min_prob are taken into account
    double out_min_prob; // minimal probability for output

    //! maximal ratio of number of base pairs divided by sequence
    //! length. This serves as a second filter on the "significant"
    //! base pairs.
    double max_bps_length_ratio;

    double max_uil_length_ratio; // max unpaired in loop length ratio
    double max_bpil_length_ratio; // max base pairs in loop length ratio

    bool no_lonely_pairs;

    int max_diff; // maximal difference for positions of alignment traces
    // (only used for ends of arcs)
    int max_diff_am; //maximal difference between two arc ends, -1 is off

    //! maximal difference for alignment traces, at arc match
    //! positions
    int max_diff_at_am;

    // only consider arc matchs where
    //   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0)
    //   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)
    // except when there is no additional computation of M matrices necessary,
    // this occurs if arcs are left-incident with larger arcs where 1 and 2 hold

    double prob_unpaired_in_loop_threshold; // threshold for prob_unpaired_in_loop
    double prob_basepair_in_loop_threshold; // threshold for prob_basepair_in_loop

    int alpha_1; //parameter for sequential score
    int alpha_2; //parameter for structural score
    int alpha_3; //parameter for stacking score
    int difference_to_opt_score;
    int min_score;
    int am_threshold;
    long int number_of_EPMs;
    bool inexact_struct_match;
    int struct_mismatch_score;

    std::string seq_constraints_A;
    std::string seq_constraints_B;

    bool no_chaining;

    int max_bp_span;

    bool relaxed_anchors;

    // ------------------------------------------------------------
    // File arguments
    std::string fileA;
    std::string fileB;
    string psFileA;
    string psFileB;
    string locarna_output;
    string output_anchor_pp;
    string clustal_output;
    string epm_list_output;
    string chained_epm_list_output;
};

//! \brief holds command line parameters
command_line_parameters clp;


// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

using namespace LocARNA;

option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","cmd_only"},
    {"help",'h',&clp.help,O_NO_ARG,0,O_NODEFAULT,"","Help"},
    {"version",'V',&clp.version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&clp.verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},
    {"quiet",'q',&clp.quiet,O_NO_ARG,0,O_NODEFAULT,"","Quiet"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},
    {"no-stacking",0,&clp.no_stacking,O_NO_ARG,0,O_NODEFAULT,"stacking",
     "Do not use stacking terms (otherwise needs stacking probs by RNAfold -p2)"},
    {"alpha_1",0,0,O_ARG_INT,&clp.alpha_1,"1","factor","Multiplier for sequential score"},
    {"alpha_2",0,0,O_ARG_INT,&clp.alpha_2,"5","factor","Multiplier for structural score"},
    {"alpha_3",0,0,O_ARG_INT,&clp.alpha_3,"5","factor",
     "Multiplier for stacking score, 0 means no stacking contribution"},
    {"struct-mismatch-score",0,0,O_ARG_INT,&clp.struct_mismatch_score,"-10","score",
     "Score for a structural mismatch (nucleotide mismatch in an arcmatch)"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},
    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,"0.01","prob","Minimal probability"},
    {"max-bps-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_bps_length_ratio,"0.0","factor",
     "Maximal ratio of #base pairs divided by sequence length (default: no effect)"},
    {"max-uil-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_uil_length_ratio,"0.0","factor",
     "Maximal ratio of #unpaired bases in loops divided by sequence length (default: no effect)"},
    {"max-bpil-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_bpil_length_ratio,"0.0","factor",
     "Maximal ratio of #base pairs in loops divided by loop length (default: no effect)"},
    {"max-diff-am",'D',0,O_ARG_INT,&clp.max_diff_am,"30","diff",
     "Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&clp.max_diff,"-1","diff",
     "Maximal difference for alignment traces"},
    {"max-diff-at-am",0,0,O_ARG_INT,&clp.max_diff_at_am,"-1","diff",
     "Maximal difference for alignment traces, only at arc match positions"},
    {"prob_unpaired_in_loop_threshold",0,0,O_ARG_DOUBLE,&clp.prob_unpaired_in_loop_threshold,"0.01","prob",
     "Threshold for prob_unpaired_in_loop"},
    {"prob_basepair_in_loop_threshold",0,0,O_ARG_DOUBLE,&clp.prob_basepair_in_loop_threshold,"0.01","prob",
     "Threshold for prob_basepair_in_loop"},


    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Output"},
    {"out-min-prob",'p',0,O_ARG_DOUBLE,&clp.out_min_prob,"0.0005","prob",
     "Minimal probability for output (min-prob overrides if smaller)"},
    {"output-ps", 0,&clp.postscript_output,O_NO_ARG,0,O_NODEFAULT,"",
     "Output best EPM chain as colored postscript"},
    {"PS_fileA",'a',0,O_ARG_STRING,&clp.psFileA,"","file","Postscript output file for sequence A"},
    {"PS_fileB",'b',0,O_ARG_STRING,&clp.psFileB,"","file","Postscript output file for sequence B"},
    {"output-locarna",'o',0,O_ARG_STRING,&clp.locarna_output,"","file",
     "Fasta file with anchor constraints for locarna"},
    {"output-anchor-pp",0,0,O_ARG_STRING,&clp.output_anchor_pp,"","fileroot",
     "PP files <fileroot>_A.pp and <fileroot>_B.pp, merging input PPs and anchor constraints from chaining"},
    {"output-clustal",0,0,O_ARG_STRING,&clp.clustal_output,"","file",
     "Write file with chain as alignment in clustalw format"},
    {"output-epm-list",0,0,O_ARG_STRING,&clp.epm_list_output,"","file","A list of all found epms"},
    {"output-chained-epm-list",0,0,O_ARG_STRING,&clp.chained_epm_list_output,"","file",
     "A list of all EPMs that are present in the chain"},
    {"inexact-struct-match",0,&clp.inexact_struct_match,O_NO_ARG,0,O_NODEFAULT,"bool",
     "Allow inexact structure matches"},
    {"add-filter",0,&clp.add_filter,O_NO_ARG,0,O_NODEFAULT,"bool",
     "Apply an additional filter to enumerate only EPMs that are maximally extended (only inexact)"},
    {"no-chaining",0,&clp.no_chaining,O_NO_ARG,0,O_NODEFAULT,"bool",
     "Do not use the chaining algorithm to find best overall chain"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Suboptimal traceback"},
    {"subopt",0,&clp.subopt,O_NO_ARG,0,O_NODEFAULT,"subopt_traceback","Use the suboptimal traceback"},
    {"diff-to-opt-score",0,0,O_ARG_INT,&clp.difference_to_opt_score,"-1","threshold",
     "Threshold for suboptimal traceback"},
    {"min-score",0,0,O_ARG_INT,&clp.min_score,"90","score","Minimal score of a traced EPM"},
    {"number-of-EPMs",0,0,O_ARG_INT,&clp.number_of_EPMs,"100","threshold",
     "Maximal number of EPMs for the suboptimal traceback"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},
    {"noLP",0,&clp.no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"bool","use --noLP option for folding"},
    {"maxBPspan",0,0,O_ARG_INT,&clp.max_bp_span,"-1","span","Limit maximum base pair span (default=off)"},
    {"relaxed-anchors",0,&clp.relaxed_anchors,O_NO_ARG,0,O_NODEFAULT,"","Relax anchor constraints (default=off)"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Miscellaneous"},
    {"stopwatch",0,&clp.stopwatch,O_NO_ARG,0,O_NODEFAULT,"","Print run time information."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Input files"},
    {"",0,0,O_ARG_STRING,&clp.fileA,O_NODEFAULT,"Input 1","Input file 1"},
    {"",0,0,O_ARG_STRING,&clp.fileB,O_NODEFAULT,"Input 2","Input file 2"},

    {"",0,0,O_TEXT,0,O_NODEFAULT,"",
     "The two input files <Input 1> and <Input 2> specify the input \
sequences in various, automatically detected formats. Accepted formats are: Fasta, Clustal, Stockholm \
LocARNA PP, ViennaRNA postscript dotplot. Unless in-loop                \
probabilities are provided (only possible in LocARNA PP), base pair     \
probabilities are computed by partition function folding. Clustal, Stockholm, \
and PP input can contain constraints."
    },

    {"",0,0,0,0,O_NODEFAULT,"",""}
};

// ------------------------------------------------------------
// ------------------------------------------------------------
// MAIN

int
main(int argc, char **argv) {
    stopwatch.start("total");

    // ------------------------------------------------------------
    // Process options

    bool process_success=process_options(argc,argv,my_options);

    if (clp.help) {

        cout << "exparna_p: fast exact structure local matching of RNAs."<<endl<<endl;

        //cout << VERSION_STRING<<endl;

        print_help(argv[0],my_options);

        cout << "Report bugs to <schmiedc (at) informatik.uni-freiburg.de>."<<endl<<endl;
        return 0;
    }

    if (clp.quiet) { clp.verbose=false;} // quiet overrides verbose

    if (clp.version || clp.verbose) {
        cout << "exparna_p (" << VERSION_STRING << ")" <<endl;
        if (clp.version) return 0; else cout <<endl;
    }

    if (!process_success) {
        std::cerr << "ERROR --- "
                  <<O_error_msg<<std::endl;
        print_usage(argv[0],my_options);
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

    //if no stacking should be considered, set the parameter for stacking to 0
    if(clp.no_stacking){
        clp.alpha_3 = 0;
    }

    if(clp.no_chaining && clp.chained_epm_list_output.size()>0){
        cout << "Enable chaining in order to output chained epm list " << endl;
        clp.no_chaining = false;
    }

    // no filtering needed if we do exact matching
    if(!clp.inexact_struct_match){
        if(clp.add_filter) cout << "Disable filtering as only exact matches are considered " << endl;
        clp.add_filter = false;
    }

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    PFoldParams pfparams(clp.no_lonely_pairs,(!clp.no_stacking),clp.max_bp_span,2);

    ExtRnaData *rna_dataA=0;
    try {
        rna_dataA = new ExtRnaData(clp.fileA,
                                   std::min(clp.min_prob,clp.out_min_prob),
                                   clp.prob_basepair_in_loop_threshold,
                                   clp.prob_unpaired_in_loop_threshold,
                                   clp.max_bps_length_ratio,
                                   clp.max_uil_length_ratio,
                                   clp.max_bpil_length_ratio,
                                   pfparams);
    } catch (failure &f) {
        std::cerr << "ERROR: failed to read from file "<<clp.fileA <<std::endl
                  << "       "<< f.what() <<std::endl;
        return -1;
    }

    ExtRnaData *rna_dataB=0;
    try {
        rna_dataB = new ExtRnaData(clp.fileB,
                                   std::min(clp.min_prob,clp.out_min_prob),
                                   clp.prob_basepair_in_loop_threshold,
                                   clp.prob_unpaired_in_loop_threshold,
                                   clp.max_bps_length_ratio,
                                   clp.max_uil_length_ratio,
                                   clp.max_bpil_length_ratio,
                                   pfparams);
    } catch (failure &f) {
        std::cerr << "ERROR: failed to read from file "<<clp.fileB <<std::endl
                  << "       "<< f.what() <<std::endl;
        if (rna_dataA) delete rna_dataA;
        return -1;
    }

    const Sequence &seqA=rna_dataA->sequence();
    const Sequence &seqB=rna_dataB->sequence();

    const MultipleAlignment &maA = seqA;
    const MultipleAlignment &maB = seqB;

    // --------------------
    // handle max_diff restriction

    TraceController trace_controller(seqA,seqB,NULL,clp.max_diff);

    // ------------------------------------------------------------
    // Handle constraints (optionally)
    size_t lenA = seqA.length();
    size_t lenB = seqB.length();

    AnchorConstraints seq_constraints(lenA,
                                      seqA.has_annotation(MultipleAlignment::AnnoType::anchors)
                                      ?seqA.annotation(MultipleAlignment::AnnoType::anchors).single_string()
                                      :"",
                                      lenB,
                                      seqB.has_annotation(MultipleAlignment::AnnoType::anchors)
                                      ?seqB.annotation(MultipleAlignment::AnnoType::anchors).single_string()
                                      :"",
                                      !clp.relaxed_anchors
                                      );

    if (clp.verbose) {
        if (! seq_constraints.empty()) {
            std::cout << "Found sequence constraints."<<std::endl;
        }
    }

    // ----------------------------------------
    // construct set of relevant arc matches
    //

    // initialize from RnaData
    ArcMatches *arc_matches = new ArcMatches(*rna_dataA,
                                             *rna_dataB,
                                             clp.min_prob,
                                             (clp.max_diff_am!=-1)
                                             ? (size_type)clp.max_diff_am
                                             : std::max(seqA.length(),seqB.length()),
                                             (clp.max_diff_at_am!=-1)
                                             ? (size_type)clp.max_diff_at_am
                                             : std::max(seqA.length(),seqB.length()),
                                             trace_controller,
                                             seq_constraints
                                             );

    const BasePairs &bpsA = arc_matches->get_base_pairsA();
    const BasePairs &bpsB = arc_matches->get_base_pairsB();


    // ----------------------------------------
    // report on input in verbose mode
    if (clp.verbose) MainHelper::report_input(seqA,seqB,*arc_matches);

    // ------------------------------------------------------------
    // construct datastructures to handle sparse matrices
    //

    //time_t start_mapping = time (NULL);
    SparsificationMapper sparse_mapperA(bpsA,
                                        *rna_dataA,
                                        clp.prob_unpaired_in_loop_threshold,
                                        clp.prob_basepair_in_loop_threshold,
                                        false
                                        );

    SparsificationMapper sparse_mapperB(bpsB,
                                        *rna_dataB,
                                        clp.prob_unpaired_in_loop_threshold,
                                        clp.prob_basepair_in_loop_threshold,
                                        false
                                        );

    SparseTraceController sparse_trace_controller(sparse_mapperA,sparse_mapperB,trace_controller);

    PatternPairMap myEPMs;

    ExactMatcher em(seqA,
                    seqB,
                    *rna_dataA,
                    *rna_dataB,
                    *arc_matches,
                    sparse_trace_controller,
                    myEPMs,
                    clp.alpha_1,
                    clp.alpha_2,
                    clp.alpha_3,
                    clp.difference_to_opt_score,
                    clp.min_score,
                    clp.number_of_EPMs,
                    clp.inexact_struct_match,
                    clp.struct_mismatch_score,
                    clp.add_filter,
                    clp.verbose
                    );

#ifndef NDEBUG
    if(clp.verbose) cout << "test arcmatch score... " << endl;
    em.test_arcmatch_score();
#endif

    // ------------------------------------------------------------
    // Compute Exact Matchings
    //

    stopwatch.start("EPMcomp");
    //compute matrices for finding best and enumerating all matchings
    em.compute_arcmatch_score();

    // ------------------------------------------------------------
    // Traceback
    //
    if (DO_TRACE) {

        if (clp.verbose) {
            if(clp.subopt)  cout << endl << "start suboptimal traceback..." << endl;
            else cout << endl << "start heuristic traceback..." << endl;
        }

        em.trace_EPMs(clp.subopt);

        stopwatch.stop("EPMcomp");

        if(clp.epm_list_output.size()>0){
            if (clp.verbose) { cout << "write list of traced EPMs in file..." << endl;}
            ofstream out_EPM_file (clp.epm_list_output.c_str());
            out_EPM_file << myEPMs.getList() << endl;
            out_EPM_file.close();
        }
    }

    // ------------------------------------------------------------
    // Chaining
    //
    if(!clp.no_chaining){

        if (clp.verbose) {cout << "Start chaining..." << endl;}
        stopwatch.start("chaining");

        PatternPairMap myLCSEPM;
        LCSEPM myChaining(seqA, seqB, myEPMs, myLCSEPM);

        //begin chaining algorithm
        myChaining.calculateLCSEPM(clp.quiet);

        stopwatch.stop("chaining");

        //output chained EPMs to PS files
        if(clp.postscript_output){
            if (clp.verbose) { cout << "write EPM chain as colored postscripts..." << endl;}
            if (clp.psFileA.size()==0){clp.psFileA = seqA.seqentry(0).name()+"_EPMs.ps";}
            if (clp.psFileB.size()==0){clp.psFileB = seqB.seqentry(0).name()+"_EPMs.ps";}

            myChaining.MapToPS(maA.consensus_sequence(), maB.consensus_sequence(), myLCSEPM, clp.psFileA,clp.psFileB);
        }

        if(clp.locarna_output.size()>0){
            if (clp.verbose) { cout << "write locarna anchor constraints..." << endl;}
            myChaining.output_locarna(maA.consensus_sequence(), maB.consensus_sequence(), clp.locarna_output);
        }

        if(clp.output_anchor_pp.size()>0){
            if (clp.verbose) { cout << "write pp files with anchor constraints..." << endl;}

            pair<SequenceAnnotation,SequenceAnnotation> anchors = myChaining.anchor_annotation();

            rna_dataA->set_anchors(anchors.first);
            rna_dataB->set_anchors(anchors.second);

            std::ofstream outA((clp.output_anchor_pp+"_A.pp").c_str());
            std::ofstream outB((clp.output_anchor_pp+"_B.pp").c_str());

            rna_dataA->write_pp(outA);
            rna_dataB->write_pp(outB);

            outA.close();
            outB.close();

        }

        if (clp.clustal_output.size()>0){
            if (clp.verbose) { cout << "write chain as clustal alignment..." << endl;}
            myChaining.output_clustal(clp.clustal_output);
        }

        if(clp.chained_epm_list_output.size()>0){
            if (clp.verbose) { cout << "write list of chained EPMs in file..." << endl;}
            ofstream out_chained_EPM_file (clp.chained_epm_list_output.c_str());
            out_chained_EPM_file << myLCSEPM.getList() << endl;
            out_chained_EPM_file.close();
        }
    }

    // ----------------------------------------
    // DONE

    delete arc_matches;

    delete rna_dataA;
    delete rna_dataB;

    if (clp.verbose) { cout << "... Exparna-P finished!" << endl << endl;}
    stopwatch.stop("total");

    return 0;
}
