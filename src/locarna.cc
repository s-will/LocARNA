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

//using namespace std;
using namespace LocARNA;

//! Version string (from configure.ac via autoconf system)
const std::string 
VERSION_STRING = (std::string)PACKAGE_STRING; 

// ------------------------------------------------------------
// Parameter


// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

//! \brief Switch on/off trace back
//! @note never made it into command line (since no one asked me kindly)
const bool DO_TRACE=true;

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
    
    int kbest_k; //!< kbest_k
    int subopt_threshold; //!< subopt_threshold

    std::string normalized_help =
        "Normalized local alignment with parameter L";
    bool normalized; //!< whether to do normalized alignment

    std::string penalized_help = 
        "Penalized local alignment with penalty PP";
    bool penalized; //!< whether to penalize alignment positions
    int position_penalty; //!< position-wise penalty

    std::string normalized_L_help =
        "Parameter L for normalized local alignment. "
        "Larger values produce larger alignments.";
    int normalized_L; //!< normalized_L

    std::string score_components_help =
        "Output components of the score (experimental)";
    bool score_components; //!< whether to report score components
};

//! \brief holds command line parameters of locarna  
command_line_parameters clp;

//! defines command line parameters
option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","cmd_only"},

    {"help",'h',&clp.help,O_NO_ARG,0,O_NODEFAULT,"",clp.help_help},
    {"galaxy-xml",0,&clp.galaxy_xml,O_NO_ARG,0,
     O_NODEFAULT,"",clp.galaxy_xml_help},
    {"version",'V',&clp.version,O_NO_ARG,0,O_NODEFAULT,"",clp.version_help},
    {"verbose",'v',&clp.verbose,O_NO_ARG,0,O_NODEFAULT,"",clp.verbose_help},
    {"quiet",'q',&clp.quiet,O_NO_ARG,0,O_NODEFAULT,"",clp.quiet_help},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"indel",'i',0,O_ARG_INT,&clp.indel,"-350","score",clp.indel_help},
    {"indel-opening",0,0,O_ARG_INT,&clp.indel_opening,
     "-500","score",clp.indel_opening_help},
    {"ribosum-file",0,0,O_ARG_STRING,&clp.ribosum_file,
     "RIBOSUM85_60","f",clp.ribosum_file_help},
    {"use-ribosum",0,0,O_ARG_BOOL,&clp.use_ribosum,
     "true","bool",clp.use_ribosum_help},
    {"match",'m',0,O_ARG_INT,&clp.match,"50","score",clp.match_help},
    {"mismatch",'M',0,O_ARG_INT,&clp.mismatch,"0","score",clp.mismatch_help},
    {"unpaired-penalty",0,0,O_ARG_INT,&clp.unpaired_penalty,
     "0","score",clp.unpaired_penalty_help},
    {"struct-weight",'s',0,O_ARG_INT,&clp.struct_weight,
     "200","score",clp.struct_weight_help},
    {"exp-prob",'e',&clp.exp_prob_given,O_ARG_DOUBLE,&clp.exp_prob,
     O_NODEFAULT,"prob",clp.exp_prob_help},
    {"tau",'t',0,O_ARG_INT,&clp.tau,"0","factor",clp.tau_help},
    {"exclusion",'E',0,O_ARG_INT,&clp.exclusion,"0","score",clp.exclusion_help},
    {"stacking",0,&clp.stacking,O_NO_ARG,0,O_NODEFAULT,"",clp.stacking_help},
    {"new-stacking",0,&clp.new_stacking,O_NO_ARG,0,
     O_NODEFAULT,"",clp.new_stacking_help},   

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Locality"},

    {"struct-local",0,&clp.struct_local_given,O_ARG_BOOL,&clp.struct_local,
     "false","bool",clp.struct_local_help},
    {"sequ-local",0,&clp.sequ_local_given,O_ARG_BOOL,&clp.sequ_local,
     "false","bool",clp.sequ_local_help},
    {"free-endgaps",0,0,O_ARG_STRING,&clp.free_endgaps,
     "----","spec",clp.free_endgaps_help},
    {"normalized",0,&clp.normalized,O_ARG_INT,&clp.normalized_L,
     "0","L",clp.normalized_help},
    
    {"penalized",0,&clp.penalized,O_ARG_INT,&clp.position_penalty,
     "0","PP",clp.penalized_help},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Output"},

    {"width",'w',0,O_ARG_INT,&clp.width,"120","columns",clp.width_help},
    {"clustal",0,&clp.clustal_given,O_ARG_STRING,&clp.clustal,
     O_NODEFAULT,"file",clp.clustal_help},
    {"stockholm",0,&clp.stockholm_given,O_ARG_STRING,&clp.stockholm,
     O_NODEFAULT,"file",clp.stockholm_help},
    {"pp",0,&clp.pp_given,O_ARG_STRING,&clp.pp,O_NODEFAULT,"file",clp.pp_help},
    {"alifold-consensus-dp",0,&clp.alifold_consensus_dp,O_NO_ARG,0,
     O_NODEFAULT,"",clp.alifold_consensus_dp_help},
    {"consensus-structure",0,0,O_ARG_STRING,&clp.cons_struct_type,
     "alifold","type",clp.cons_struct_type_help},
    {"local-output",'L',&clp.local_output,O_NO_ARG,0,
     O_NODEFAULT,"",clp.local_output_help},
    {"local-file-output",0,&clp.local_file_output,O_NO_ARG,0,
     O_NODEFAULT,"",clp.local_file_output_help},
    {"pos-output",'P',&clp.pos_output,O_NO_ARG,0,
     O_NODEFAULT,"",clp.pos_output_help},
    {"write-structure",0,&clp.write_structure,O_NO_ARG,0,
     O_NODEFAULT,"",clp.write_structure_help},
    {"score-components",0,&clp.score_components,O_NO_ARG,0,
     O_NODEFAULT,"",clp.score_components_help},
    {"stopwatch",0,&clp.stopwatch,O_NO_ARG,0,
     O_NODEFAULT,"",clp.stopwatch_help},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"",
     "Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,
     "0.0005","prob",clp.min_prob_help},
    {"max-bps-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_bps_length_ratio,
     "0.0","factor",clp.max_bps_length_ratio_help},
    {"max-diff-am",'D',0,O_ARG_INT,&clp.max_diff_am,
     "-1","diff",clp.max_diff_am_help},
    {"max-diff",'d',0,O_ARG_INT,&clp.max_diff,"-1","diff",clp.max_diff_help},
    {"max-diff-at-am",0,0,O_ARG_INT,&clp.max_diff_at_am,
     "-1","diff",clp.max_diff_at_am_help},
    {"max-diff-aln",0,0,O_ARG_STRING,&clp.max_diff_alignment_file,
     "","aln file",clp.max_diff_alignment_file_help},
    {"max-diff-pw-aln",0,0,O_ARG_STRING,&clp.max_diff_pw_alignment,
     "","alignment",clp.max_diff_pw_alignment_help},
    {"max-diff-relax",0,&clp.max_diff_relax,O_NO_ARG,0,
     O_NODEFAULT,"",clp.max_diff_relax_help},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Special sauce options"},
    {"kbest",0,&clp.subopt,O_ARG_INT,&clp.kbest_k,
     "-1","k","Enumerate k-best alignments"},
    {"better",0,&clp.subopt,O_ARG_INT,&clp.subopt_threshold,
     "-1000000","t","Enumerate alignments better threshold t"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","MEA score"},

    {"mea-alignment",0,&clp.mea_alignment,O_NO_ARG,0,O_NODEFAULT,"",
     clp.mea_alignment_help},
    {"match-prob-method",0,0,O_ARG_INT,&clp.match_prob_method,
     "0","int",clp.match_prob_method_help},
    {"probcons-file",0,&clp.probcons_file_given,O_ARG_STRING,&clp.probcons_file,
     O_NODEFAULT,"file",
     clp.probcons_file_help},
    {"temperature-alipf",0,0,O_ARG_INT,&clp.temperature_alipf,
     "150","int",clp.temperature_alipf_help},
    {"pf-struct-weight",0,0,O_ARG_INT,&clp.pf_struct_weight,
     "200","weight",clp.pf_struct_weight_help},
    {"mea-gapcost",0,&clp.mea_gapcost,O_NO_ARG,0,
     O_NODEFAULT,"","Use gap cost in mea alignment"},   
    {"mea-alpha",0,0,O_ARG_INT,&clp.mea_alpha,
     "0","weight",clp.mea_alpha_help},
    {"mea-beta",0,0,O_ARG_INT,&clp.mea_beta,
     "200","weight",clp.mea_beta_help},
    {"mea-gamma",0,0,O_ARG_INT,&clp.mea_gamma,
     "100","weight",clp.mea_gamma_help},
    {"probability-scale",0,0,O_ARG_INT,&clp.probability_scale,
     "10000","scale",
     clp.probability_scale_help},
    {"write-match-probs",0,&clp.write_matchprobs,
     O_ARG_STRING,&clp.matchprobs_outfile,
     O_NODEFAULT,"file",clp.write_matchprobs_help},
    {"read-match-probs",0,&clp.read_matchprobs,
     O_ARG_STRING,&clp.matchprobs_infile,
     O_NODEFAULT,"file",clp.read_matchprobs_help},
    {"write-arcmatch-scores",0,&clp.write_arcmatch_scores,
     O_ARG_STRING,&clp.arcmatch_scores_outfile,
     O_NODEFAULT,"file",clp.write_arcmatch_scores_help},
    {"read-arcmatch-scores",0,&clp.read_arcmatch_scores,O_ARG_STRING,
     &clp.arcmatch_scores_infile,O_NODEFAULT,"file",clp.read_arcmatch_scores_help},
    {"read-arcmatch-probs",0,&clp.read_arcmatch_probs,O_ARG_STRING,
     &clp.arcmatch_scores_infile,O_NODEFAULT,"file",
     clp.read_arcmatch_probs_help},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&clp.no_lonely_pairs,O_NO_ARG,0,
     O_NODEFAULT,"",clp.no_lonely_pairs_help},
    {"maxBPspan",0,0,O_ARG_INT,&clp.max_bp_span,
     "-1","span",clp.max_bp_span_help},
    {"relaxed-anchors",0,&clp.relaxed_anchors,O_NO_ARG,0,
     O_NODEFAULT,"",clp.relaxed_anchors_help},
    
    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Hidden Options"},

    {"ribofit",0,0,O_ARG_BOOL,&clp.ribofit,"false","bool",clp.ribofit_help},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"",
     "Input files"},

    {"",0,0,O_ARG_STRING,&clp.fileA,O_NODEFAULT,"Input 1",clp.fileA_help},
    {"",0,0,O_ARG_STRING,&clp.fileB,O_NODEFAULT,"Input 2",clp.fileB_help},

    {"",0,0,O_TEXT,0,O_NODEFAULT,"",clp.files_help},

    {"",0,0,0,0,O_NODEFAULT,"",""}
};


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

    bool process_success=process_options(argc,argv,my_options);

    if (clp.help) {
	std::cout
            << "locarna - pairwise (global and local) alignment of RNA."
            <<std::endl<<std::endl;
	
	//std::cout << VERSION_STRING<<std::endl<<std::endl;

	print_help(argv[0],my_options);

	std::cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."
                  <<std::endl<<std::endl;
	return 0;
    }

    if (clp.quiet) { clp.verbose=false;} // quiet overrides verbose

    if (clp.galaxy_xml) {
    	print_galaxy_xml((char *)"locarna",my_options);
    	return 0;
    }

    if (clp.version || clp.verbose) {
	std::cout << VERSION_STRING<<std::endl;
	if (clp.version) return 0; else std::cout <<std::endl;
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
    if (clp.read_arcmatch_scores && clp.read_arcmatch_probs) {
	std::cerr << "One cannot specify files for arc match scores "
                  << "and probabilities simultaneously."
                  << std::endl;
	return -1;
    }
    
    if (clp.probability_scale<=0) {
	std::cerr << "Probability scale must be greater 0."<<std::endl;
	return -1;
    }
    
    if (clp.struct_weight<0) {
	std::cerr << "Structure weight must be greater equal 0."<<std::endl;
	return -1;
    }

    if (clp.normalized && clp.penalized) {
    	std::cerr << "One cannot specify penalized and normalized "
                  << "simultaneously."<<std::endl;
    	return -1;
    }

    // ----------------------------------------
    // temporarily turn off stacking unless background prob is set
    //
    if (clp.stacking && !clp.exp_prob_given) {
	std::cerr << "WARNING: stacking turned off. "
		  << "Stacking requires setting a background probability "
		  << "explicitely (option --exp-prob)." << std::endl;
	clp.stacking=false;
    }
    
    // ------------------------------------------------------------ 
    // if normalized or penalized alignment shall be computed,
    // automatically turn on sequ_local unless sequ_local mode was
    // explicitly specified
    //
    // important: in the Aligner class, we rely on sequ_local==true in
    // normalized alignment mode
    if (clp.normalized || clp.penalized) {
	if(!clp.sequ_local_given) {
	    clp.sequ_local = true;
	} else {
	    if (!clp.sequ_local) {
		std::cerr << "ERROR: Cannot run normalized alignment "
                          << "without --sequ_local on."
                          << std::endl;
		return -1;
	    }
	}

	if (clp.struct_local_given && clp.struct_local) {
	    std::cerr << "ERROR: Normalized structure local alignment "
                      << "not supported."
                      << std::endl;
	    return -1;
	} else {
	    clp.struct_local=false;
	}
	
    }
    
    // ----------------------------------------  
    // Ribosum matrix
    //
    RibosumFreq *ribosum;
    Ribofit *ribofit;
    
    MainHelper::init_ribo_matrix(clp,&ribosum,&ribofit);
    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    PFoldParams pfparams(clp.no_lonely_pairs, clp.stacking || clp.new_stacking, clp.max_bp_span, 2);
    
    RnaData *rna_dataA=0;
    try {
	rna_dataA = new RnaData(clp.fileA,
				clp.min_prob,
				clp.max_bps_length_ratio,
				pfparams);
    } catch (failure &f) {
	std::cerr << "ERROR:\tfailed to read from file "<<clp.fileA <<std::endl
		  << "\t"<< f.what() <<std::endl;
	return -1;
    }
    
    RnaData *rna_dataB=0;
    try {
	rna_dataB = new RnaData(clp.fileB,
				clp.min_prob,
				clp.max_bps_length_ratio,
				pfparams);
    } catch (failure &f) {
	std::cerr << "ERROR: failed to read from file "<<clp.fileB <<std::endl
		  << "       "<< f.what() <<std::endl;
	if (rna_dataA) delete rna_dataA;
	return -1;
    }
    
    const Sequence &seqA=rna_dataA->sequence();
    const Sequence &seqB=rna_dataB->sequence();
    
    size_type lenA=seqA.length();
    size_type lenB=seqB.length();

    // --------------------
    // handle max_diff restriction  
    
    // missing: proper error handling in case that lenA, lenB, and
    // max_diff_pw_alignment/max_diff_alignment_file are incompatible
    
    // do inconsistency checking for max_diff_pw_alignment and
    // max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment!="" && clp.max_diff_alignment_file!="") {
	std::cerr <<"Cannot simultaneously use options --max-diff-pw-alignment"
                  <<" and --max-diff-alignment-file."<<std::endl;
	return -1;
    }

    // construct TraceController and check inconsistency for with
    // multiplicity of sequences
    //

    MultipleAlignment *multiple_ref_alignment=NULL;
    
    if (clp.max_diff_alignment_file!="") {
	multiple_ref_alignment = new MultipleAlignment(clp.max_diff_alignment_file);
    } else if (clp.max_diff_pw_alignment!="") {
	if ( seqA.num_of_rows()!=1 || seqB.num_of_rows()!=1 ) {
	    std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning "
                      << "of alignments."
                      << std::endl;
	    return -1;
	}
	
	std::vector<std::string> alistr;
	split_at_separator(clp.max_diff_pw_alignment,'&',alistr);
	
	if (alistr.size()!=2) {
	    std::cerr << "Invalid argument to --max-diff-pw-alignemnt; require "
                      << "exactly one '&' separating the alignment strings."
		      << std::endl; 
	    return -1;
	}
    
	if (alistr[0].length() != alistr[1].length()) {
	    std::cerr << "Invalid argument to --max-diff-pw-alignemnt;"
                      <<" alignment strings have unequal lengths."
		      << std::endl; 
	    return -1;
	}
	
	multiple_ref_alignment = new MultipleAlignment(seqA.seqentry(0).name(),
						       seqB.seqentry(0).name(),
						       alistr[0],
						       alistr[1]);
    }
    
    TraceController trace_controller(seqA,seqB,multiple_ref_alignment,
                                     clp.max_diff,clp.max_diff_relax);
    
    
    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    AnchorConstraints 
        seq_constraints(lenA,
                        seqA.annotation(MultipleAlignment::AnnoType::anchors)
                        .single_string(),
                        lenB,
                        seqB.annotation(MultipleAlignment::AnnoType::anchors)
                        .single_string(),
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
    ArcMatches *arc_matches;
    
    // ------------------------------------------------------------
    // handle reading and writing of arcmatch_scores
    //
    // (needed for mea alignment with probabilistic consistency
    // transformation of arc match scores)
    //
    if (clp.read_arcmatch_scores || clp.read_arcmatch_probs) {
	if (clp.verbose) {
	    std::cout << "Read arcmatch scores from file "
                      << clp.arcmatch_scores_infile << "." <<std::endl;
	}
	arc_matches = new ArcMatches(seqA,
				     seqB,
				     clp.arcmatch_scores_infile,
				     clp.read_arcmatch_probs
				     ? ((clp.mea_beta*clp.probability_scale)
                                        /100)
				     : -1,
				     clp.max_diff_am!=-1
				     ? (size_type)clp.max_diff_am
				     : std::max(lenA,lenB),
				     clp.max_diff_at_am!=-1
				     ? (size_type)clp.max_diff_at_am
				     : std::max(lenA,lenB),
				     trace_controller,
				     seq_constraints
				     );
    } else {
	// initialize from RnaData
	arc_matches = new ArcMatches(*rna_dataA,
				     *rna_dataB,
				     clp.min_prob,
				     clp.max_diff_am!=-1
				     ? (size_type)clp.max_diff_am
				     : std::max(lenA,lenB),
				     clp.max_diff_at_am!=-1
				     ? (size_type)clp.max_diff_at_am
				     : std::max(lenA,lenB),
				     trace_controller,
				     seq_constraints
				     );
    }
    
    // note: arcmatches has to be assigned ( unless new failed, which
    // we don't catch )
    
    // ----------------------------------------
    // report on input in verbose mode
    if (clp.verbose) MainHelper::report_input(seqA,seqB,*arc_matches);


    // ------------------------------------------------------------
    // Sequence match probabilities (for MEA-Alignment)
    //
    // perform parameter consistency checks
    if (clp.read_matchprobs && !clp.mea_alignment) {
        std::cerr << "Warning: clp.read_matchprobs ignored for "<<
            "non-mea alignment.\n"; 
    }
    if ( (clp.write_matchprobs || clp.mea_alignment)
         && ribosum==NULL && ribofit==NULL
         ) {
        std::cerr << "ERROR: Ribosum/fit is required for mea_alignment and "
                  << "computing matchprobs."<<std::endl;
        exit(-1);
    }
    //
    MatchProbs *match_probs=0L;
    if (clp.write_matchprobs || clp.mea_alignment) {
        match_probs = MainHelper::init_match_probs(clp,
                                                   rna_dataA,rna_dataB,
                                                   ribosum,ribofit);
    }
    if (clp.write_matchprobs) {
        MainHelper::write_match_probs(clp,match_probs);
        if (!clp.write_arcmatch_scores) { return 0; } // return from main()
    }
    //

    // ----------------------------------------
    // Scoring Parameter
    //
    double my_exp_probA = clp.exp_prob_given?clp.exp_prob:prob_exp_f(lenA);
    double my_exp_probB = clp.exp_prob_given?clp.exp_prob:prob_exp_f(lenB);

    ScoringParams
        scoring_params(clp.match,
                       clp.mismatch,
                       // In true mea alignment gaps are only scored
                       // for computing base match probs.
                       // Consequently, we set the indel and indel
                       // opening cost to 0 for the case of mea
                       // alignment!
                       (clp.mea_alignment && !clp.mea_gapcost)
                       ?0
                       :(clp.indel
                         * (clp.mea_gapcost?clp.probability_scale/100:1)),
                       0, // indel__loop_score, for
                       // consistency and least modification to
                       // locarna.cc has been set to zero
                       (clp.mea_alignment && !clp.mea_gapcost)
                       ?0
                       :(clp.indel_opening
                         * (clp.mea_gapcost?clp.probability_scale/100:1)),
                       0, // indel_opening_loop_score, for
                       // consistency and least modification to
                       // locarna.cc has been set to zero
                       ribosum,
                       ribofit,
                       clp.unpaired_penalty,
                       clp.struct_weight,
                       clp.tau,
                       clp.exclusion,
                       my_exp_probA,
                       my_exp_probB,
                       clp.temperature_alipf,
                       clp.stacking,
                       clp.new_stacking,
                       clp.mea_alignment,
                       clp.mea_alpha,
                       clp.mea_beta,
                       clp.mea_gamma,
                       clp.probability_scale
                       );
    // ------------------------------------------------------------
    // Construct scoring
    
    Scoring scoring(seqA,
		    seqB,
		    *rna_dataA,
		    *rna_dataB,
		    *arc_matches,
		    match_probs,
		    scoring_params,
		    false // no Boltzmann weights (as required for LocARNA-P)
		    );    

    if (clp.write_arcmatch_scores) {
	if (clp.verbose) {
	    std::cout << "Write arcmatch scores to file "
                      << clp.arcmatch_scores_outfile<<" and exit."<<std::endl;
	}
	arc_matches->write_arcmatch_scores(clp.arcmatch_scores_outfile,scoring);
	return 0;
    }
        

    // ------------------------------------------------------------
    // Computation of the alignment score
    //
    
    // initialize aligner object, which does the alignment computation
    Aligner aligner = Aligner::create()
	. seqA(seqA)
	. seqB(seqB)
	. arc_matches(*arc_matches)
	. scoring(scoring)
	. no_lonely_pairs(clp.no_lonely_pairs)
	. struct_local(clp.struct_local)
	. sequ_local(clp.sequ_local)
	. free_endgaps(clp.free_endgaps)
	. max_diff_am(clp.max_diff_am)
	. max_diff_at_am(clp.max_diff_at_am)
	. trace_controller(trace_controller)
	. stacking(clp.stacking || clp.new_stacking)
	. constraints(seq_constraints);

    // enumerate suboptimal alignments (using interval splitting)
    if (clp.subopt) {
	aligner.suboptimal(clp.kbest_k,
			   clp.subopt_threshold,
			   clp.normalized,
			   clp.normalized_L,
			   clp.width,
			   clp.verbose,
			   clp.local_output,
			   clp.pos_output,
			   clp.write_structure
			   );
	return 0;
    }
    
    infty_score_t score;

    // if option --normalized <L> is given, then do normalized local alignemnt
    if (clp.normalized) {
		
	score = aligner.normalized_align(clp.normalized_L,clp.verbose);
	
    } else if (clp.penalized) {
    	score = aligner.penalized_align(clp.position_penalty);
    }

    else {
	
	// ========== STANDARD CASE ==========
	
	// otherwise compute the best alignment
	score = aligner.align();
	
    }
    
    // ----------------------------------------
    // report score
    //
    if (!clp.quiet) {
        std::cout << "Score: "<<score<<std::endl<<std::endl;
    }

    // ------------------------------------------------------------
    // Traceback
    //
    if ((!clp.normalized && !clp.penalized) && DO_TRACE) {
	
	aligner.trace();
	
	// for debugging:
	// aligner.get_alignment().write_debug(std::cout);


	// if score components should be reported
	if (clp.score_components) {
	    // get alignment and compute contributions by sequence
	    // similarity and gap cost
	    const Alignment &alignment = aligner.get_alignment();
	    
	    // experimental: I do the calculations here until I
	    // decide where to put the code
	    
	    std::cout << "WARNING: reporting score components is still "
                      << "experimental."
                      << std::endl;
	    std::cout << "         This does not work properly for certain "
                      << "non-standard standard scoring types, "
                      << "e.g. free end gaps." << std::endl;

	    const Alignment::edges_t edges = alignment.alignment_edges(false);

	    score_t seq_sim=0;
	    score_t gap_cost=0; // count linear component of gap cost
	    score_t gap_numA=0; // count number of gap openings in A
	    score_t gap_numB=0; // count number of gap openings in B
	    
	    bool openA=false; // is a gap open in A
	    bool openB=false; // is a gap open in B
	    
	    for (size_t k=0; k< edges.size(); k++) {
		if (edges.first[k].is_pos() && edges.second[k].is_pos()) {
		    seq_sim += scoring.basematch(edges.first[k],
                                                 edges.second[k]);
		}
		if (edges.first[k].is_gap()) {
		    if (!openA) {
			gap_numA++;
			openA=true;
		    }
		    gap_cost += scoring.gapA(edges.second[k]);
		} else {
		    openA=false;
		}

		if (edges.second[k].is_gap()) {
		    if (!openB) {
			gap_numB++;
			openB=true;
		    }
		    gap_cost += scoring.gapB(edges.first[k]);
		} else {
		    openB=false;
		}
	    }
	    
	    
	    score_t total_gap_cost = gap_cost
                + (gap_numA+gap_numB)*scoring.indel_opening();
	    
	    std::cout << "#SCORE total        "
                      << std::setw(8) << score<<std::endl;
	    std::cout << "#SCORE seq_sim      "
                      << std::setw(8) << seq_sim<<std::endl;
	    std::cout << "#SCORE gap_penalty  "
		      << std::setw(8) << total_gap_cost<<std::endl;
	    std::cout << "#SCORE str_contrib  "
                      << std::setw(8) << score-seq_sim-total_gap_cost
                      << std::endl;
	}
    }
    
    int return_code=0;

    if (clp.normalized || clp.penalized || DO_TRACE) {
	// if we did a trace (one way or the other)
	
	const Alignment &alignment = aligner.get_alignment();

        // ----------------------------------------
	// write alignment in different output formats to files
	//
        
        std::string consensus_structure=""; 
        
        RnaData *consensus =
            MainHelper::consensus(clp,
                                  pfparams,
                                  my_exp_probA, my_exp_probB,
                                  rna_dataA, rna_dataB,
                                  alignment,
                                  consensus_structure);
        
        return_code = MainHelper::write_alignment(clp,
                                                  score,
                                                  consensus_structure,
                                                  consensus,
                                                  alignment,
                                                  multiple_ref_alignment);
        
        // ----------------------------------------
        // write alignment to screen
                
	if (clp.pos_output) {
	    std::cout << "HIT "<<score<<" "
		      <<alignment.local_startA()<<" "
		      <<alignment.local_startB()<<" "
		      <<alignment.local_endA()<<" "
		      <<alignment.local_endB()<<" "
		      <<std::endl;
            std::cout << std::endl;
	}
        
	if ((!clp.pos_output && !clp.quiet) || clp.local_output) {
	    MultipleAlignment ma(alignment,clp.local_output);
	    
	    if (clp.write_structure) {
		// annotate multiple alignment with structures
                std::string structureA =
                    alignment.dot_bracket_structureA(clp.local_output);
                std::string structureB =
                    alignment.dot_bracket_structureB(clp.local_output);
                ma.prepend(MultipleAlignment::SeqEntry("",structureA));
		ma.append(MultipleAlignment::SeqEntry("",structureB));
	    }
	    
	    if (clp.pos_output) {
		std::cout  << "\t+" << alignment.local_startA() << std::endl
			   << "\t+" << alignment.local_startB() << std::endl
			   << std::endl;
                std::cout << std::endl;
            }
            
            if (consensus_structure!="") {
                if (clp.local_output!=clp.local_file_output) {
                    // recompute consensus structure
                    clp.local_file_output=clp.local_output;
                    
                    if (consensus) { delete consensus; }
                    consensus =
                        MainHelper::consensus(clp,
                                              pfparams,
                                              my_exp_probA, my_exp_probB,
                                              rna_dataA, rna_dataB,
                                              alignment,
                                              consensus_structure);
                    clp.local_file_output=!clp.local_output;
                }
                ma.append(MultipleAlignment::SeqEntry(clp.cons_struct_type,
                                                      consensus_structure));
            }
	    ma.write(std::cout,clp.width,
                     MultipleAlignment::FormatType::CLUSTAL);

	    if (clp.pos_output) {
		std::cout  << std::endl 
			   << "\t+" << alignment.local_endA() << std::endl
			   << "\t+" << alignment.local_endB() << std::endl
			   << std::endl;
	    }
            
        }

	if (!clp.quiet) {
            std::cout<<std::endl;
	}
        
        if (consensus) { delete consensus; }

    }

    // ----------------------------------------
    //  clean up
    //
    if (match_probs) delete match_probs;
    delete arc_matches;
    if (multiple_ref_alignment) delete multiple_ref_alignment;
    if (ribofit) delete ribofit;
    if (ribosum) delete ribosum;

    delete rna_dataA;
    delete rna_dataB;
        
    stopwatch.stop("total");

    // ----------------------------------------
    // DONE
    return return_code;
}
