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
struct command_line_parameters : public MainHelper::std_command_line_parameters {
    bool opt_subopt; //!< find suboptimal solution (either k-best or all solutions better than a threshold)
    
    int kbest_k; //!< kbest_k
    int subopt_threshold; //!< subopt_threshold

    bool opt_mea_alignment; //!< whether to perform mea alignment
    bool opt_mea_gapcost; //!< whether to use mea gapcost
    int mea_alpha; //!< mea alpha
    int mea_beta; //!< mea beta
    int mea_gamma; //!< mea gamma
    int probability_scale; //!< probability scale

    bool opt_normalized; //!< whether to do normalized alignment

    bool opt_penalized; //!< whether to penalize alignment positions

    int normalized_L; //!< normalized_L

    bool opt_score_components; //!< whether to report score components
};

//! \brief holds command line parameters of locarna  
command_line_parameters clp;

//! defines command line parameters
option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","cmd_only"},

    {"help",'h',&clp.opt_help,O_NO_ARG,0,O_NODEFAULT,"","Help"},
    {"galaxy-xml",0,&clp.opt_galaxy_xml,O_NO_ARG,0,O_NODEFAULT,"","Galaxy xml wrapper"},
    {"version",'V',&clp.opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&clp.opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring_parameters"},

    {"match",'m',0,O_ARG_INT,&clp.match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&clp.mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&clp.ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&clp.use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&clp.indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&clp.indel_opening_score,"-500","score","Indel opening score"},
    {"unpaired_penalty",0,0,O_ARG_INT,&clp.unpaired_penalty,"0","score","Penalty for unpaired bases"},
//    {"position_penalty",0,0,O_ARG_INT,&clp.position_penalty,"0","score","Penalty for each alignment position"},
    {"struct-weight",'s',0,O_ARG_INT,&clp.struct_weight,"200","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&clp.opt_exp_prob,O_ARG_DOUBLE,&clp.exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&clp.tau_factor,"0","factor","Tau factor in percent"},
    {"exclusion",'E',0,O_ARG_INT,&clp.exclusion_score,"0","score","Exclusion weight"},
    {"stacking",0,&clp.opt_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use stacking terms (needs stack-probs by RNAfold -p2)"},
    {"new-stacking",0,&clp.opt_new_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use new stacking terms (needs stack-probs by RNAfold -p2)"},   

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Locality_type"},

    {"struct-local",0,&clp.struct_local_given,O_ARG_BOOL,&clp.struct_local,"false","bool","Structure local"},
    {"sequ-local",0,&clp.sequ_local_given,O_ARG_BOOL,&clp.sequ_local,"false","bool","Sequence local"},
    {"free-endgaps",0,0,O_ARG_STRING,&clp.free_endgaps,"----","spec","Whether and which end gaps are free. order: L1,R1,L2,R2"},
    {"normalized",0,&clp.opt_normalized,O_ARG_INT,&clp.normalized_L,"0","L","Normalized local alignment with parameter L"},
    {"penalized",0,&clp.opt_penalized,O_ARG_INT,&clp.position_penalty,"0","PP","Penalized local alignment with penalty PP"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling_output"},

    {"width",'w',0,O_ARG_INT,&clp.output_width,"120","columns","Output width"},
    {"clustal",0,&clp.opt_clustal_out,O_ARG_STRING,&clp.clustal_out,O_NODEFAULT,"file","ClustalW (aln) output"},
    {"stockholm",0,&clp.opt_stockholm_out,O_ARG_STRING,&clp.stockholm_out,O_NODEFAULT,"file","Stockholm output"},
    {"pp",0,&clp.opt_pp_out,O_ARG_STRING,&clp.pp_out,O_NODEFAULT,"file","PP output"},
    {"alifold-consensus-dp",0,&clp.opt_alifold_consensus_dp,O_NO_ARG,0,O_NODEFAULT,"","Compute consensus dot plot by alifold"},
    {"local-output",'L',&clp.opt_local_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment (to standard out)"},
    {"local-file-output",0,&clp.opt_local_file_output,O_NO_ARG,0,O_NODEFAULT,"","Write only local sub-alignment to output files"},
    {"pos-output",'P',&clp.opt_pos_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment positions"},
    {"write-structure",0,&clp.opt_write_structure,O_NO_ARG,0,O_NODEFAULT,"","Write guidance structure in output"},
    {"score-components",0,&clp.opt_score_components,O_NO_ARG,0,O_NODEFAULT,"","Output components of the score (experimental)"},
    {"stopwatch",0,&clp.opt_stopwatch,O_NO_ARG,0,O_NODEFAULT,"","Print run time information."},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,"0.0005","prob","Minimal probability"},
    {"max-bps-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_bps_length_ratio,"0.0","factor","Maximal ratio of #base pairs divided by sequence length (default: no effect)"},
    {"max-diff-am",'D',0,O_ARG_INT,&clp.max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&clp.max_diff,"-1","diff","Maximal difference for alignment traces"},
    {"max-diff-at-am",0,0,O_ARG_INT,&clp.max_diff_at_am,"-1","diff","Maximal difference for alignment traces, only at arc match positions"},
    {"max-diff-aln",0,0,O_ARG_STRING,&clp.max_diff_alignment_file,"","aln file","Maximal difference relative to given alignment (file in clustalw format))"},
    {"max-diff-pw-aln",0,0,O_ARG_STRING,&clp.max_diff_pw_alignment,"","alignment","Maximal difference relative to given alignment (string, delim=AMPERSAND)"},
    {"max-diff-relax",0,&clp.opt_max_diff_relax,O_NO_ARG,0,O_NODEFAULT,"","Relax deviation constraints in multiple aligmnent"},
    {"min-am-prob",'a',0,O_ARG_DOUBLE,&clp.min_am_prob,"0.0005","amprob","Minimal Arc-match probability"},
    {"min-bm-prob",'b',0,O_ARG_DOUBLE,&clp.min_bm_prob,"0.0005","bmprob","Minimal Base-match probability"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Special sauce options"},
    {"kbest",0,&clp.opt_subopt,O_ARG_INT,&clp.kbest_k,"-1","k","Enumerate k-best alignments"},
    {"better",0,&clp.opt_subopt,O_ARG_INT,&clp.subopt_threshold,"-1000000","t","Enumerate alignments better threshold t"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","MEA_score controlling options"},

    {"mea-alignment",0,&clp.opt_mea_alignment,O_NO_ARG,0,O_NODEFAULT,"","Do MEA alignment"},
    {"probcons-file",0,&clp.opt_probcons_file,O_ARG_STRING,&clp.probcons_file,O_NODEFAULT,"file","Probcons parameter file"},

    {"match-prob-method",0,0,O_ARG_INT,&clp.match_prob_method,"0","int","Method for computation of match probs"},
    {"temperature",0,0,O_ARG_INT,&clp.temperature,"150","int","Temperature for PF-computation"},
    {"pf-struct-weight",0,0,O_ARG_INT,&clp.pf_struct_weight,"200","weight","Structure weight in PF-computation"},

    {"mea-gapcost",0,&clp.opt_mea_gapcost,O_NO_ARG,0,O_NODEFAULT,"","Use gap cost in mea alignment"},   
    {"mea-alpha",0,0,O_ARG_INT,&clp.mea_alpha,"0","weight","Weight alpha for MEA"},
    {"mea-beta",0,0,O_ARG_INT,&clp.mea_beta,"200","weight","Weight beta for MEA"},
    {"mea-gamma",0,0,O_ARG_INT,&clp.mea_gamma,"100","weight","Weight gamma for MEA"},
    {"probability-scale",0,0,O_ARG_INT,&clp.probability_scale,"10000","scale","Scale for probabilities/resolution of mea score"},

    {"write-match-probs",0,&clp.opt_write_matchprobs,O_ARG_STRING,&clp.matchprobs_file,O_NODEFAULT,"file","Write match probs to file (don't align!)"},
    {"read-match-probs",0,&clp.opt_read_matchprobs,O_ARG_STRING,&clp.matchprobs_file,O_NODEFAULT,"file","Read match probabilities from file"},

    {"write-arcmatch-scores",0,&clp.opt_write_arcmatch_scores,O_ARG_STRING,&clp.arcmatch_scores_file,O_NODEFAULT,"file","Write arcmatch scores (don't align!)"},
    {"read-arcmatch-scores",0,&clp.opt_read_arcmatch_scores,O_ARG_STRING,&clp.arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch scores"},
    {"read-arcmatch-probs",0,&clp.opt_read_arcmatch_probs,O_ARG_STRING,&clp.arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch probabilities (weight by mea_beta/100)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&clp.no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs"},
    // {"anchorA",0,0,O_ARG_STRING,&clp.seq_anchors_A,"","string","Anchor constraints sequence A"},
    // {"anchorB",0,0,O_ARG_STRING,&clp.seq_anchors_B,"","string","Anchor constraints sequence B"},
    //{"ignore-constraints",0,&clp.opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints input files"},

    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Hidden Options"},
    {"ribofit",0,0,O_ARG_BOOL,&clp.opt_ribofit,"false","bool","Use Ribofit base and arc match scores (overrides ribosum)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Input_files RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&clp.fileA,O_NODEFAULT,"input1","Input file 1"},
    {"",0,0,O_ARG_STRING,&clp.fileB,O_NODEFAULT,"input2","Input file 2"},
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

    if (clp.opt_help) {
	std::cout << "locarna - a tool for pairwise (global and local) alignment of RNA."<<std::endl<<std::endl;
	
	//std::cout << VERSION_STRING<<std::endl<<std::endl;

	print_help(argv[0],my_options);

	std::cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<std::endl<<std::endl;
	return 0;
    }

    if (clp.opt_galaxy_xml) {
    	print_galaxy_xml((char *)"locarna",my_options);
    	return 0;
    }

    if (clp.opt_version || clp.opt_verbose) {
	std::cout << VERSION_STRING<<std::endl;
	if (clp.opt_version) return 0; else std::cout <<std::endl;
    }

    if (!process_success) {
	std::cerr << "ERROR --- "
		  <<O_error_msg<<std::endl;
	printf("USAGE: ");
	print_usage(argv[0],my_options);
	printf("\n");
	return -1;
    }

    if (clp.opt_stopwatch) {
	stopwatch.set_print_on_exit(true);
    }
    
    if (clp.opt_verbose) {
	print_options(my_options);
    }
    
    // ------------------------------------------------------------
    // parameter consistency
    if (clp.opt_read_arcmatch_scores && clp.opt_read_arcmatch_probs) {
	std::cerr << "One cannot specify arc match score and probabilities file simultaneously."<<std::endl;
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

    if (clp.opt_normalized && clp.opt_penalized) {
    	std::cerr << "One cannot specify penalized and normalized simultaneously."<<std::endl;
    	return -1;
    }

    // ----------------------------------------
    // temporarily turn off stacking unless background prob is set
    //
    if (clp.opt_stacking && !clp.opt_exp_prob) {
	std::cerr << "WARNING: stacking turned off. "
		  << "Stacking requires setting a background probability "
		  << "explicitely (option --exp-prob)." << std::endl;
	clp.opt_stacking=false;
    }
    
    // ------------------------------------------------------------
    // if normalized or penalized alignment shall be computed, automatically turn on
    // sequ_local unless sequ_local mode was explicitly specified
    //
    // important: in the Aligner class, we rely on sequ_local==true in normalized alignment mode
    if (clp.opt_normalized || clp.opt_penalized) {
	if(!clp.sequ_local_given) {
	    clp.sequ_local = true;
	} else {
	    if (!clp.sequ_local) {
		std::cerr 
		    << "ERROR: Cannot run normalized alignment with sequence local alignment turned off."<<std::endl;
		return -1;
	    }
	}

	if (clp.struct_local_given && clp.struct_local) {
	    std::cerr 
		<< "ERROR: Normalized structure local alignment not supported."
		<<std::endl;
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

    PFoldParams pfparams(clp.no_lonely_pairs, clp.opt_stacking || clp.opt_new_stacking);
    
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
    
    // missing: proper error handling in case that lenA, lenB, and max_diff_pw_alignment/max_diff_alignment_file are incompatible 
    
    // do inconsistency checking for max_diff_pw_alignment and max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment!="" && clp.max_diff_alignment_file!="") {
	std::cerr <<"Cannot simultaneously use options --max-diff-pw-alignment"
                  <<" and --max-diff-alignment-file."<<std::endl;
	return -1;
    }

    // construct TraceController and check inconsistency for with multiplicity of sequences
    //

    MultipleAlignment *multiple_ref_alignment=NULL;
    
    if (clp.max_diff_alignment_file!="") {
	multiple_ref_alignment = new MultipleAlignment(clp.max_diff_alignment_file);
    } else if (clp.max_diff_pw_alignment!="") {
	if ( seqA.num_of_rows()!=1 || seqB.num_of_rows()!=1 ) {
	    std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning of alignments." << std::endl;
	    return -1;
	}
	
	std::vector<std::string> alistr;
	split_at_separator(clp.max_diff_pw_alignment,'&',alistr);
	
	if (alistr.size()!=2) {
	    std::cerr << "Invalid argument to --max-diff-pw-alignemnt; require exactly one '&' separating the alignment strings."
		      << std::endl; 
	    return -1;
	}
    
	if (alistr[0].length() != alistr[1].length()) {
	    std::cerr << "Invalid argument to --max-diff-pw-alignemnt; alignment strings have unequal lengths."
		      << std::endl; 
	    return -1;
	}
	
	multiple_ref_alignment = new MultipleAlignment(seqA.seqentry(0).name(),
						       seqB.seqentry(0).name(),
						       alistr[0],
						       alistr[1]);
    }
    
    TraceController trace_controller(seqA,seqB,multiple_ref_alignment,clp.max_diff,clp.opt_max_diff_relax);
    
    
    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    AnchorConstraints seq_constraints(lenA,
				      seqA.annotation(MultipleAlignment::AnnoType::anchors).single_string(),
				      lenB,
				      seqB.annotation(MultipleAlignment::AnnoType::anchors).single_string());
    
    if (clp.opt_verbose) {
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
    // (needed for mea alignment with probabilistic consistency transformation of arc match scores)
    //
    if (clp.opt_read_arcmatch_scores || clp.opt_read_arcmatch_probs) {
	if (clp.opt_verbose) {
	    std::cout << "Read arcmatch scores from file " << clp.arcmatch_scores_file << "." <<std::endl;
	}
	arc_matches = new ArcMatches(seqA,
				     seqB,
				     clp.arcmatch_scores_file,
				     clp.opt_read_arcmatch_probs
				     ? ((clp.mea_beta*clp.probability_scale)/100)
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
    
    // note: arcmatches has to be assigned ( unless new failed, which we don't catch )
    
    // ----------------------------------------
    // report on input in verbose mode
    if (clp.opt_verbose) MainHelper::report_input(seqA,seqB,*arc_matches);


    // ------------------------------------------------------------
    // Sequence match probabilities (for MEA-Alignment)
    //
    // perform parameter consistency checks
    if (clp.opt_read_matchprobs && !clp.opt_mea_alignment) {
        std::cerr << "Warning: clp.opt_read_matchprobs ignored for non-mea alignment.\n"; 
    }
    if ( (clp.opt_write_matchprobs || clp.opt_mea_alignment)
         && ribosum==NULL && ribofit==NULL
         ) {
        std::cerr << "ERROR: Ribosum/fit is required for mea_alignment and computing matchprobs."<<std::endl;
        exit(-1);
    }
    //
    MatchProbs *match_probs=0L;
    if (clp.opt_write_matchprobs || clp.opt_mea_alignment) {
        match_probs = MainHelper::init_match_probs(clp,rna_dataA,rna_dataB,ribosum,ribofit);
    }
    if (clp.opt_write_matchprobs) {
        MainHelper::write_match_probs(clp,match_probs);
        if (!clp.opt_write_arcmatch_scores) { return 0; } // return from main()
    }
    //

    // ----------------------------------------
    // Scoring Parameter
    //
    double my_exp_probA = clp.opt_exp_prob?clp.exp_prob:prob_exp_f(lenA);
    double my_exp_probB = clp.opt_exp_prob?clp.exp_prob:prob_exp_f(lenB);

    ScoringParams scoring_params(clp.match_score,
				 clp.mismatch_score,
				 // In true mea alignment gaps are only 
				 // scored for computing base match probs.
				 // Consequently, we set the indel and indel opening cost to 0
				 // for the case of mea alignment!
				 (clp.opt_mea_alignment && !clp.opt_mea_gapcost)
				 ?0
				 :clp.indel_score * (clp.opt_mea_gapcost?clp.probability_scale/100:1),
				 0, // indel__loop_score, for
                                    // consistency and least
                                    // modification to locarna.cc has
                                    // been set to zero
				 (clp.opt_mea_alignment && !clp.opt_mea_gapcost)
				 ?0
				 :clp.indel_opening_score * (clp.opt_mea_gapcost?clp.probability_scale/100:1),
				 0, // indel_opening_loop_score, for
                                    // consistency and least
                                    // modification to locarna.cc has
                                    // been set to zero
				 ribosum,
				 ribofit,
				 clp.unpaired_penalty,
				 clp.struct_weight,
				 clp.tau_factor,
				 clp.exclusion_score,
				 my_exp_probA,
				 my_exp_probB,
				 clp.temperature,
				 clp.opt_stacking,
				 clp.opt_new_stacking,
				 clp.opt_mea_alignment,
				 clp.mea_alpha,
				 clp.mea_beta,
				 clp.mea_gamma,
				 clp.probability_scale
				 );
    // --------------------------------------------------------------------------------
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

    if (clp.opt_write_arcmatch_scores) {
	if (clp.opt_verbose) {
	    std::cout << "Write arcmatch scores to file "<< clp.arcmatch_scores_file<<" and exit."<<std::endl;
	}
	arc_matches->write_arcmatch_scores(clp.arcmatch_scores_file,scoring);
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
	. min_am_prob(clp.min_am_prob)
	. min_bm_prob(clp.min_bm_prob)
	. stacking(clp.opt_stacking || clp.opt_new_stacking)
	. constraints(seq_constraints);

    // enumerate suboptimal alignments (using interval splitting)
    if (clp.opt_subopt) {
	aligner.suboptimal(clp.kbest_k,
			   clp.subopt_threshold,
			   clp.opt_normalized,
			   clp.normalized_L,
			   clp.output_width,
			   clp.opt_verbose,
			   clp.opt_local_output,
			   clp.opt_pos_output,
			   clp.opt_write_structure
			   );
	return 0;
    }
    
    infty_score_t score;

    // if option --normalized <L> is given, then do normalized local alignemnt
    if (clp.opt_normalized) {
		
	score = aligner.normalized_align(clp.normalized_L,clp.opt_verbose);
	
    } else if (clp.opt_penalized) {
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
    std::cout << "Score: "<<score<<std::endl;


    // ------------------------------------------------------------
    // Traceback
    //
    if ((!clp.opt_normalized && !clp.opt_penalized) && DO_TRACE) {
	    
	if (clp.opt_verbose) {
	    std::cout << "Traceback."<<std::endl;
	}
	
	aligner.trace();
	
	// for debugging:
	// aligner.get_alignment().write_debug(std::cout);


	// if score components should be reported
	if (clp.opt_score_components) {
	    // get alignment and compute contributions by sequence similarity and gap cost
	    const Alignment &alignment = aligner.get_alignment();
	    
	    // under development: I do the calculations here until I decide where to put the code
	    
	    std::cout << "WARNING: reporting score components is still experimental." << std::endl;
	    std::cout << "         This does not work properly for some less standard scoring types, e.g. free end gaps." << std::endl;

	    const Alignment::edges_t edges = alignment.alignment_edges(false);
	    
	    score_t seq_sim=0;
	    score_t gap_cost=0; // count linear component of gap cost
	    score_t gap_numA=0; // count number of gap openings in A
	    score_t gap_numB=0; // count number of gap openings in B
	    
	    bool openA=false; // is a gap open in A
	    bool openB=false; // is a gap open in B
	    
	    for (size_t k=0; k< edges.size(); k++) {
		if (edges.first[k].is_pos() && edges.second[k].is_pos()) {
		    seq_sim += scoring.basematch(edges.first[k],edges.second[k]);
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
	    
	    
	    score_t total_gap_cost = gap_cost + (gap_numA+gap_numB)*scoring.indel_opening();
	    
	    std::cout << "#SCORE total        " << std::setw(8) << score<<std::endl;
	    std::cout << "#SCORE seq_sim      " << std::setw(8) << seq_sim<<std::endl;
	    std::cout << "#SCORE gap_penalty  "
		//<<gap_cost<<"+"<<(gap_numA+gap_numB)<<"*("<<scoring.indel_opening()<<") = "
		      << std::setw(8) << total_gap_cost<<std::endl;
	    std::cout << "#SCORE str_contrib  " << std::setw(8) << score-seq_sim-total_gap_cost<<std::endl;
	}
    }
    
    int return_code=0;

    if (clp.opt_normalized || clp.opt_penalized || DO_TRACE) {
	// if we did a trace (one way or the other)
	
	const Alignment &alignment = aligner.get_alignment();
	
	if (clp.opt_pos_output) {
	    std::cout << "HIT "<<score<<" "
		      <<alignment.local_startA()<<" "
		      <<alignment.local_startB()<<" "
		      <<alignment.local_endA()<<" "
		      <<alignment.local_endB()<<" "
		      <<std::endl;
	} 
	if (!clp.opt_pos_output || clp.opt_local_output) {
	    MultipleAlignment ma(alignment,clp.opt_local_output);
	    
	    if (clp.opt_write_structure) {
		// annotate multiple alignment with structures
		ma.prepend(MultipleAlignment::SeqEntry("",
						       alignment.dot_bracket_structureA(clp.opt_local_output)));
		ma.append(MultipleAlignment::SeqEntry("",
						      alignment.dot_bracket_structureB(clp.opt_local_output)));
	    }
	    
	    if (clp.opt_pos_output) {
		std::cout  << std::endl 
			   << "\t+" << alignment.local_startA() << std::endl
			   << "\t+" << alignment.local_startB() << std::endl
			   << std::endl;
	    }
	    
	    std::cout << std::endl;
	    ma.write(std::cout,clp.output_width,MultipleAlignment::FormatType::CLUSTAL);

	    if (clp.opt_pos_output) {
		std::cout  << std::endl 
			   << "\t+" << alignment.local_endA() << std::endl
			   << "\t+" << alignment.local_endB() << std::endl
			   << std::endl;
	    }

	}

	std::cout<<std::endl;
	
        // ----------------------------------------
	// write alignment in different output formats
	//
        
        RnaData *consensus=0L;
        std::string consensus_structure=""; 
        
	consensus = MainHelper::consensus(clp,
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
