/**
 * \file locarna.cc
 *
 * \brief Defines main function of locarna
 *
 * LocARNA: global and LOCal Alignment of RNA
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de> 
 *               2005-2011
 * 
 */


#include <iostream>
#include <fstream>
#include <vector>

#include <memory> // for auto_ptr

//#include <math.h>

#include <LocARNA/sequence.hh>
#include <LocARNA/basepairs.hh>
#include <LocARNA/alignment.hh>
#include <LocARNA/aligner_n.hh>
#include <LocARNA/evaluator.hh>
#include <LocARNA/rna_data.hh>
#include <LocARNA/arc_matches.hh>
#include <LocARNA/match_probs.hh>
#include <LocARNA/ribosum.hh>
#include <LocARNA/anchor_constraints.hh>
#include <LocARNA/trace_controller.hh>
#include <LocARNA/ribosum85_60.icc>


using namespace std;
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
#include <LocARNA/options.hh>

//! \brief Switch on/off trace back
//! @note never made it into command line
const bool DO_TRACE=true;

//! \brief Structure for command line parameters of locarna
//!
//! Encapsulating all command line parameters in a common structure
//! avoids name conflicts and makes downstream code more informative.
//!
struct command_line_parameters {
    //! only pairs with a probability of at least min_prob are taken into account
    double min_prob; 
    
    int match_score; //!< match score
    
    int mismatch_score; //!< mismatch score
    
    int indel_score; //!< indel extension score
    
    int indel_opening_score; //!< indel opening score
    
    int temperature; //!< temperature
    
    int struct_weight; //!< structure weight

    //! contribution of sequence similarity in an arc match (in percent)
    int tau_factor;

    bool no_lonely_pairs; //!< no lonely pairs option

    //! allow exclusions for maximizing alignment of connected substructures
    bool struct_local;

    bool sequ_local; //!< maximize alignment of subsequences

    //! specification of free end gaps, order left end sequence 1,
    //! right 1, left 2, right 2 e.g. "+---" allows free end gaps at
    //! the left end of the first alignment string ; "----" forbids
    //! free end gaps
    std::string free_endgaps; 
    
    //! maximal difference for positions of alignment
    //! traces (only used for ends of arcs)
    int max_diff; 
    
    //! maximal difference between two arc ends, -1 is off
    int max_diff_am;

    //! pairwise reference alignment for max-diff heuristic,
    //!separator &
    std::string max_diff_pw_alignment;
    
    //! reference alignment for max-diff heuristic, name of clustalw
    //! format file
    std::string max_diff_alignment_file;

    //! use relaxed variant of max diff with reference alignment
    bool opt_max_diff_relax; 

    //! Score contribution per exclusion
    //! set to zero for unrestricted structure locality
    int exclusion_score; 
    
    //! expected probability of a base pair (null-model)
    double exp_prob;
    
    //! expected probability given?
    bool opt_exp_prob;

    //! width of alignment output
    int output_width;

    // ------------------------------------------------------------
    // File arguments
    
    //! first input file 
    std::string file1;
    
    //! second input file
    std::string file2;

    std::string clustal_out; //!< name of clustal output file

    bool opt_clustal_out; //!< whether to write clustal output to file

    std::string pp_out; //!< name of pp output file
    
    bool opt_pp_out; //!< whether to write pp output to file
    
    bool opt_alifold_consensus_dp; //!< whether to compute consensus dp by alifold

    bool opt_help; //!< whether to print help
    bool opt_version; //!< whether to print version
    bool opt_verbose; //!< whether to print verbose output
    bool opt_local_output; //!< whether to use local output
    bool opt_pos_output; //!< whether to output positions

    bool opt_write_structure; //!< whether to write structure


    bool opt_stacking; //!< whether to stacking

    std::string ribosum_file; //!< ribosum_file
    bool use_ribosum; //!< use_ribosum

    bool opt_probcons_file; //!< whether to probcons_file
    std::string probcons_file; //!< probcons_file

    bool opt_mea_alignment; //!< whether to mea_alignment

    bool opt_write_matchprobs; //!< whether to write_matchprobs
    bool opt_read_matchprobs; //!< whether to read_matchprobs
    std::string matchprobs_file; //!< matchprobs_file

    bool opt_write_arcmatch_scores; //!< whether to write_arcmatch_scores

    bool opt_read_arcmatch_scores; //!< whether to read arcmatch scores
    bool opt_read_arcmatch_probs; //!< whether to read arcmatch probabilities
    std::string arcmatch_scores_file; //!< arcmatch scores file

    int match_prob_method; //!< method for computing match probabilities

    double min_am_prob; //!< only matched arc-pair with a probability of at least min_am_prob are taken into account
    double min_bm_prob; //!< only matched base-pair with a probability of at least min_bm_prob are taken into account

    bool opt_subopt; //!< find suboptimal solution (either k-best or all solutions better than a threshold)

    int kbest_k; //!< kbest_k
    int subopt_threshold; //!< subopt_threshold

    std::string seq_constraints_A; //!< seq_constraints_A
    std::string seq_constraints_B; //!< seq_constraints_B

    bool opt_ignore_constraints; //!< whether to ignore_constraints

    int pf_struct_weight; //!< pf_struct_weight

    bool opt_mea_gapcost; //!< whether to use mea gapcost
    int mea_alpha; //!< mea alpha
    int mea_beta; //!< mea beta
    int mea_gamma; //!< mea gamma
    int probability_scale; //!< probability scale

    bool opt_eval; //!< whether to evaluate

    bool opt_normalized; //!< whether to do normalized alignment
    int normalized_L; //!< normalized_L
};


//! \brief holds command line parameters of locarna  
command_line_parameters clp;


//! defines command line parameters
option_def my_options[] = {
    {"help",'h',&clp.opt_help,O_NO_ARG,0,O_NODEFAULT,"","Help"},
    {"version",'V',&clp.opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&clp.opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"match",'m',0,O_ARG_INT,&clp.match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&clp.mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&clp.ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&clp.use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&clp.indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&clp.indel_opening_score,"-500","score","Indel opening score"},
    {"struct-weight",'s',0,O_ARG_INT,&clp.struct_weight,"200","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&clp.opt_exp_prob,O_ARG_DOUBLE,&clp.exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&clp.tau_factor,"0","factor","Tau factor in percent"},
    {"exclusion",'E',0,O_ARG_INT,&clp.exclusion_score,"0","score","Exclusion weight"},
    {"stacking",0,&clp.opt_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use stacking terms (needs stack-probs by RNAfold -p2)"},   

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Type of locality"},

    {"struct-local",0,0,O_ARG_BOOL,&clp.struct_local,"false","bool","Structure local"},
    {"sequ-local",0,0,O_ARG_BOOL,&clp.sequ_local,"false","bool","Sequence local"},
    {"free-endgaps",0,0,O_ARG_STRING,&clp.free_endgaps,"----","spec","Whether and which end gaps are free. order: L1,R1,L2,R2"},
    {"normalized",0,&clp.opt_normalized,O_ARG_INT,&clp.normalized_L,"0","L","Normalized local alignment with parameter L"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling output"},

    {"width",'w',0,O_ARG_INT,&clp.output_width,"120","columns","Output width"},
    {"clustal",0,&clp.opt_clustal_out,O_ARG_STRING,&clp.clustal_out,O_NODEFAULT,"file","Clustal output"},
    {"pp",0,&clp.opt_pp_out,O_ARG_STRING,&clp.pp_out,O_NODEFAULT,"file","PP output"},
    
#ifdef HAVE_LIBRNA
    {"alifold-consensus-dp",0,&clp.opt_alifold_consensus_dp,O_NO_ARG,0,O_NODEFAULT,"","Compute consensus dot plot by alifold"},
#endif

    {"local-output",'L',&clp.opt_local_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment"},
    {"pos-output",'P',&clp.opt_pos_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment positions"},
    {"write-structure",0,&clp.opt_write_structure,O_NO_ARG,0,O_NODEFAULT,"","Write guidance structure in output"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&clp.max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&clp.max_diff,"-1","diff","Maximal difference for alignment traces"},
    {"max-diff-aln",0,0,O_ARG_STRING,&clp.max_diff_alignment_file,"","aln file","Maximal difference relative to given alignment (file in clustalw format))"},
    {"max-diff-pw-aln",0,0,O_ARG_STRING,&clp.max_diff_pw_alignment,"","alignment","Maximal difference relative to given alignment (string, delim=&)"},
    {"max-diff-relax",0,&clp.opt_max_diff_relax,O_NO_ARG,0,O_NODEFAULT,"","Relax deviation constraints in multiple aligmnent"},
    {"min-am-prob",'a',0,O_ARG_DOUBLE,&clp.min_am_prob,"0.0005","amprob","Minimal Arc-match probability"},
    {"min-bm-prob",'b',0,O_ARG_DOUBLE,&clp.min_bm_prob,"0.0005","bmprob","Minimal Base-match probability"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Special sauce options"},
    {"kbest",0,&clp.opt_subopt,O_ARG_INT,&clp.kbest_k,"-1","k","Enumerate k-best alignments"},
    {"better",0,&clp.opt_subopt,O_ARG_INT,&clp.subopt_threshold,"-1000000","t","Enumerate alignments better threshold t"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Options for controlling MEA score"},

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
    {"anchorA",0,0,O_ARG_STRING,&clp.seq_constraints_A,"","string","Anchor constraints sequence A"},
    {"anchorB",0,0,O_ARG_STRING,&clp.seq_constraints_B,"","string","Anchor constraints sequence B"},
    {"ignore-constraints",0,&clp.opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Mode of operation"},
    {"eval",0,&clp.opt_eval,O_NO_ARG,0,O_NODEFAULT,"","Turn on evaluation mode"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&clp.file1,O_NODEFAULT,"file 1","Basepairs input file 1 (alignment in eval mode)"},
    {"",0,0,O_ARG_STRING,&clp.file2,O_NODEFAULT,"file 2","Basepairs input file 2 (dp dir in eval mode)"},
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
	cout << "locarna - a tool for pairwise (global and local) alignment of RNA!"<<endl<<endl;
	
	cout << VERSION_STRING<<endl<<endl;

	print_help(argv[0],my_options);

	cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (clp.opt_version || clp.opt_verbose) {
	cout << "locarna ("<< VERSION_STRING<<")"<<endl;
	if (clp.opt_version) exit(0); else cout <<endl;
    }

    if (!process_success) {
	std::cerr << "ERROR --- "
		  <<O_error_msg<<std::endl;
	printf("USAGE: ");
	print_usage(argv[0],my_options);
	printf("\n");
	exit(-1);
    }
    
    if (clp.opt_verbose) {
	print_options(my_options);
    }
    
    // ------------------------------------------------------------
    // parameter consistency
    if (clp.opt_read_arcmatch_scores && clp.opt_read_arcmatch_probs) {
	std::cerr << "You cannot specify arc match score and probabilities file simultaneously."<<std::endl;
	exit(-1);
    }
    
    if (clp.probability_scale<=0) {
	std::cerr << "Probability scale must be greater 0."<<std::endl;
	exit(-1);
    }
    
    if (clp.struct_weight<0) {
	std::cerr << "Structure weight must be greater equal 0."<<std::endl;
	exit(-1);
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


    // ----------------------------------------  
    // Ribosum matrix
    //
    std::auto_ptr<RibosumFreq> ribosum(NULL);
	
    if (clp.use_ribosum) {
	if (clp.ribosum_file == "RIBOSUM85_60") {
	    if (clp.opt_verbose) {
		std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	  	    ribosum = auto_ptr<RibosumFreq>(new Ribosum85_60);
	} else {
	    ribosum = auto_ptr<RibosumFreq>(new RibosumFreq(clp.ribosum_file));
	}	
	/*
	  std::cout <<" A: "<< ribosum->base_nonstruct_prob('A')
	  <<" C: "<< ribosum->base_nonstruct_prob('C')
	  <<" G: "<< ribosum->base_nonstruct_prob('G')
	  <<" U: "<< ribosum->base_nonstruct_prob('U')
	  << std::endl;
	
	  ribosum->print_basematch_scores_corrected();
	*/
    }
    
    // ----------------------------------------
    // Scoring Parameter
    //
    ScoringParams scoring_params(clp.match_score,
				 clp.mismatch_score,
				 // In true mea alignment gaps are only 
				 // scored for computing base match probs.
				 // Consequently, we set the indel and indel opening cost to 0
				 // for the case of mea alignment!
				 (clp.opt_mea_alignment && !clp.opt_mea_gapcost)
				 ?0
				 :clp.indel_score * (clp.opt_mea_gapcost?clp.probability_scale/100:1),
				 (clp.opt_mea_alignment && !clp.opt_mea_gapcost)
				 ?0
				 :clp.indel_opening_score * (clp.opt_mea_gapcost?clp.probability_scale/100:1),
				 ribosum.get(),
				 clp.struct_weight,
				 clp.tau_factor,
				 clp.exclusion_score,
				 clp.opt_exp_prob?clp.exp_prob:-1,
				 clp.temperature,
				 clp.opt_stacking,
				 clp.opt_mea_alignment,
				 clp.mea_alpha,
				 clp.mea_beta,
				 clp.mea_gamma,
				 clp.probability_scale
				 );

    
    // ----------------------------------------
    // CHOOSE MODE OF OPERATION
    //
    // whether to do evaluation
    //
    
    if (clp.opt_eval) {
	std::cout <<"Evaluation Mode"<<std::endl;
	std::cout <<"==============="<<std::endl;
	
	Evaluator e(clp.file1,clp.file2,&scoring_params);
	std::cout <<"SCORE: "<< e.eval() << std::endl;
	
	std::exit(0);
    }
    
    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
    
    PFoldParams params(clp.no_lonely_pairs,clp.opt_stacking);

    RnaData rnadataA(clp.file1,true,clp.opt_stacking,true);
    if (!rnadataA.pairProbsAvailable() || !rnadataA.inLoopProbsAvailable()) {
	if (clp.opt_verbose) {
	    std::cout << "Compute ensemble probabilities for first input sequence."
		      << std::endl;
	}
	rnadataA.computeEnsembleProbs(params,true);
    }

    RnaData rnadataB(clp.file2,true,clp.opt_stacking,true);
    if (!rnadataB.pairProbsAvailable() || !rnadataB.inLoopProbsAvailable()) {
	if (clp.opt_verbose) {
	    std::cout << "Compute ensemble probabilities for second input sequence."
		      << std::endl;
	}
	rnadataB.computeEnsembleProbs(params,true);
    }

    Sequence seqA=rnadataA.get_sequence();
    Sequence seqB=rnadataB.get_sequence();
    
    size_type lenA=seqA.length();
    size_type lenB=seqB.length();

    // --------------------
    //Forbid unsupported option of LocARNA_N
    if ( clp.struct_local )
    {
    	std::cerr << "Exclusions is not supported" << std::endl;
    	exit (-1);
    }
    if ( clp.indel_opening_score != 0 )
    {
    	std::cerr << "Affine gap cost is not supported" << std::endl;
    	exit (-1);
    }


    if( clp.no_lonely_pairs )
    {
    	std::cerr << "No lonely pairs option is not supported" << std::endl;
    	exit (-1);
    }
    if( clp.sequ_local )
    {
    	std::cerr << "Local sequence alignment is not supported" << std::endl;
    	exit (-1);
    }


    // --------------------
    // handle max_diff restriction  
    
    // missing: proper error handling in case that lenA, lenB, and max_diff_pw_alignment/max_diff_alignment_file are incompatible 
    
    // do inconsistency checking for max_diff_pw_alignment and max_diff_alignment_file
    //
    if (clp.max_diff_pw_alignment!="" && clp.max_diff_alignment_file!="") {
	std::cerr <<"Cannot simultaneously use both options --max-diff-pw-alignemnt and --max-diff-alignment-file."<<std::endl;
	exit(-1);
    }

    // construct TraceController and check inconsistency for with multiplicity of sequences
    //

    MultipleAlignment *multiple_ref_alignment=NULL;
    
    if (clp.max_diff_alignment_file!="") {
	multiple_ref_alignment = new MultipleAlignment(clp.max_diff_alignment_file);
    } else if (clp.max_diff_pw_alignment!="") {
	if ( seqA.row_number()!=1 || seqB.row_number()!=1 ) {
	    std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning of alignments." << std::endl;
	    exit(-1);
	}
	
	multiple_ref_alignment = new MultipleAlignment(seqA.names()[0],seqB.names()[0],clp.max_diff_pw_alignment);
    }

    // if (multiple_ref_alignment) {
    // 	std::cout<<"Reference aligment:"<<std::endl;
    // 	multiple_ref_alignment->print_debug(std::cout);
    // 	std::cout << std::flush;
    // }
    
    TraceController trace_controller(seqA,seqB,multiple_ref_alignment,clp.max_diff,clp.opt_max_diff_relax);
    
    
    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    std::string seqCA = clp.seq_constraints_A;
    std::string seqCB = clp.seq_constraints_B;

    if (!clp.opt_ignore_constraints) {
	if ( seqCA=="" ) seqCA = rnadataA.get_seq_constraints();
	if ( seqCB=="" ) seqCB = rnadataB.get_seq_constraints();
    }

    AnchorConstraints seq_constraints(seqA.length(),seqCA,
				      seqB.length(),seqCB);
    
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
				     (clp.opt_read_arcmatch_probs
				      ?((clp.mea_beta*clp.probability_scale)/100)
				      :-1),
				     ((clp.max_diff_am!=-1)
				      ?(size_type)clp.max_diff_am
				      :std::max(lenA,lenB)),
				     trace_controller,
				     seq_constraints
				     );
    } else {
	// initialize from RnaData
	arc_matches = new ArcMatches(rnadataA,
				     rnadataB,
				     clp.min_prob,
				     (clp.max_diff_am!=-1)?(size_type)clp.max_diff_am:std::max(lenA,lenB),
				     trace_controller,
				     seq_constraints
				     );
    }
    
    BasePairs bpsA = arc_matches->get_base_pairsA();
    BasePairs bpsB = arc_matches->get_base_pairsB();
    
    // ----------------------------------------
    // report on input in verbose mode
    if (clp.opt_verbose) {
	std::cout << "Sequence A: "<<std::endl;
	seqA.write(cout);
	std::cout<<" (Length:"<< seqA.length()<<", Basepairs:"<<bpsA.num_bps() << ")" <<std::endl;

	std::cout << "Sequence B: "<<std::endl;
	seqB.write(cout);
	std::cout<<" (Length:"<< seqB.length()<<", Basepairs:"<<bpsB.num_bps() << ")" <<std::endl;

	cout <<std::endl 
	     <<"Base Pair Matches: "<<arc_matches->num_arc_matches() << "." <<std::endl;
	// cout << "Base Identity: "<<(seq_identity(seqA,seqB)*100)<<endl; 
    }


    // ------------------------------------------------------------
    // Sequence match probabilities (for MEA-Alignment)
    //
    MatchProbs *match_probs=0L;

    if (clp.opt_read_matchprobs && !clp.opt_mea_alignment) {
	std::cerr << "Warning: clp.opt_read_matchprobs ignored for non-mea alignment.\n"; 
    }

    if (clp.opt_write_matchprobs || clp.opt_mea_alignment) {
	match_probs = new MatchProbs;

	if (!clp.use_ribosum) {
	    std::cerr << "WARNING: Ribosum scoring used for mea_alignment and computing matchprobs."<<std::endl;
	}

	if (clp.opt_read_matchprobs) {
	    match_probs->read_sparse(clp.matchprobs_file,seqA.length(),seqB.length());
	} else {
	    if (clp.match_prob_method==1) {
		if (!clp.opt_probcons_file) {
		    std::cerr << "Probcons parameter file required for pairHMM-style computation"
			      <<" of basematch probabilities."<<std::endl;
		    print_usage(argv[0],my_options);
		    std::cerr << std::endl;
		    exit(-1);
		}
		if (clp.opt_verbose) {
		    std::cout << "Compute match probabilities using pairHMM."<<std::endl; 
		}

		match_probs->pairHMM_probs(seqA,seqB,clp.probcons_file);
	    } else {
		bool sl=clp.sequ_local;
		if (clp.match_prob_method==2) sl=true;
		if (clp.match_prob_method==3) sl=false;

		if (clp.opt_verbose) {
		    std::cout << "Compute match probabilities using PF sequence alignment."<<std::endl; 
		}

		match_probs->pf_probs(rnadataA,rnadataB,
				      ribosum->get_basematch_scores(),
				      ribosum->alphabet(),
				      clp.indel_opening_score/100.0,
				      clp.indel_score/100.0,
				      clp.pf_struct_weight/100.0,
				      clp.temperature/100.0,
				      sl);
	    }
	}

	if (clp.opt_write_matchprobs) {
	    if (clp.opt_verbose) {
		std::cout << "Write match probabilities to file "<<clp.matchprobs_file<<"."<<std::endl; 
	    }

	    match_probs->write_sparse(clp.matchprobs_file,1.0/clp.probability_scale);
	    if (!clp.opt_write_arcmatch_scores) exit(0); // else we exit there!
	}
    }
   

    // ----------------------------------------
    // construct scoring
   
    Scoring scoring(seqA,seqB,arc_matches,match_probs,&scoring_params);    

    if (clp.opt_write_arcmatch_scores) {
	if (clp.opt_verbose) {
	    std::cout << "Write arcmatch scores to file "<< clp.arcmatch_scores_file<<" and exit."<<std::endl;
	}
	arc_matches->write_arcmatch_scores(clp.arcmatch_scores_file,scoring);
	exit(0);
    }
        

    // ------------------------------------------------------------
    // Computation of the alignment score
    //

    // parameter for the alignment
    AlignerParams aligner_params(clp.no_lonely_pairs,
				 clp.struct_local,
				 clp.sequ_local,
				 clp.free_endgaps,
				 trace_controller,
				 clp.max_diff_am,
				 clp.min_am_prob,
				 clp.min_bm_prob,
				 clp.opt_stacking,
				 seq_constraints
				 );
    
    // initialize aligner object, which does the alignment computation
    AlignerN aligner(seqA,seqB,*arc_matches,&aligner_params,&scoring);

    
    // enumerate suboptimal alignments (using interval splitting)
    if (clp.opt_subopt) {
    	std::cerr
      			<< "ERROR: suboptimal alignment not supported."
      			<<std::endl;
        exit(-1);

	/*aligner.suboptimal(clp.kbest_k,
			   clp.subopt_threshold,
			   clp.opt_normalized,
			   clp.normalized_L,
			   clp.output_width,
			   clp.opt_verbose,
			   clp.opt_local_output,
			   clp.opt_pos_output,
			   clp.opt_write_structure
			   );
	exit(0);
	*/
    }
    
    infty_score_t score;

    // if option --normalized <L> is given, then do normalized local alignemnt
    if (clp.opt_normalized) {
    	  std::cerr
    			<< "ERROR: Normalized alignment not supported."
    			<<std::endl;
      exit(-1);


	// do some option consistency checks and output errors
	if (clp.struct_local) {
	    std::cerr 
		<< "ERROR: Normalized structure local alignment not supported."
		<<std::endl
		<< "LocARNA ignores struct_local option."<<std::endl;
	    exit(-1);
	}
	if (!clp.sequ_local) { // important: in the Aligner class, we rely on this
	    std::cerr 
		<< "ERROR: Normalized alignment requires option --sequ_local."<<std::endl;
	    exit(-1);
	}
//	score = aligner.normalized_align(clp.normalized_L,clp.opt_verbose);
	
    } else {
	
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
    if ((!clp.opt_normalized) && DO_TRACE) {
	    
	if (clp.opt_verbose) {
	    std::cout << "Traceback."<<std::endl;
	}
	
	aligner.trace();
	
	// for debugging:
	//if (clp.opt_verbose)
	//    aligner.get_alignment().write_debug(std::cout);
    }
    
    if (clp.opt_normalized || DO_TRACE) { // if we did a trace (one way or
				      // the other)

	aligner.get_alignment().write(std::cout, 
				      clp.output_width,
				      score,
				      clp.opt_local_output,
				      clp.opt_pos_output,
				      clp.opt_write_structure
				      );
	std::cout<<endl;
	
	// test MultipleAlignment
	if (clp.opt_verbose) {
	    MultipleAlignment resultMA(aligner.get_alignment());
	    //std::cout << "MultipleAlignment"<<std::endl; 
	    //resultMA.print_debug(cout);
	    if (multiple_ref_alignment) {
		std::cout << "Deviation to reference: "<< multiple_ref_alignment->deviation(resultMA)<<std::endl;
	    }
	}
	
	// ----------------------------------------
	// optionally write output formats
	//
	if (clp.opt_clustal_out) {
	    ofstream out(clp.clustal_out.c_str());
	    if (out.good()) {
		aligner.get_alignment().write_clustal(out, clp.output_width,
						      score,
						      clp.opt_local_output,
						      clp.opt_pos_output,
						      true,
						      clp.opt_write_structure
						      );
	    } else {
		cerr << "Cannot write to "<<clp.clustal_out<<endl<<"! Exit.";
		exit(-1);
	    }
	}
	if (clp.opt_pp_out) {

	    // if compiled without vienna rna lib, deactivate alifold
	    // consensus dot plot option
#         ifndef HAVE_LIBRNA
	    clp.opt_alifold_consensus_dp=false;
#         endif

	    ofstream out(clp.pp_out.c_str());
	    if (out.good()) {
		aligner.get_alignment().
		    write_pp(out,bpsA,bpsB,
			     scoring,
			     seq_constraints,
			     clp.output_width,
			     clp.opt_alifold_consensus_dp);
	    } else {
		cerr << "Cannot write to "<<clp.pp_out<<endl<<"! Exit.";
		exit(-1);
	    }
	}
    }
    
    stopwatch.stop("total");
    
    // ----------------------------------------
    // DONE
    exit(0);
}
