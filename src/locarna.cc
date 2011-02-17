/**********************************************************************
 *
 * LocARNA: LOCal Alignment of RNA
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de> 
 *               2005-2010
 *
 **********************************************************************/


#include <iostream>
#include <fstream>
#include <vector>

#include <memory> // for auto_ptr

//#include <math.h>

#include <LocARNA/sequence.hh>
#include <LocARNA/basepairs.hh>
#include <LocARNA/alignment.hh>
#include <LocARNA/aligner.hh>
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

const std::string 
VERSION_STRING = (std::string)PACKAGE_STRING; 

// ------------------------------------------------------------
// Parameter

double min_prob; // only pairs with a probability of at least min_prob
		 // are taken into account

int match_score;
int mismatch_score;
int indel_score;
int indel_opening_score;
int temperature;
int struct_weight;

int tau_factor; // contribution of sequence similarity in an arc match
		// (in percent)

bool no_lonely_pairs; // no lonely pairs option

bool struct_local; // allow exclusions for maximizing alignment of
		   // connected substructures
bool sequ_local; // maximize alignment of subsequences

std::string free_endgaps; //!< specification of free end gaps,
// order left end sequence 1, right 1, left 2, right 2
// e.g. "+---" allows free end gaps at the left end of the first alignment string
// ; "----" forbids free end gaps

const bool DO_TRACE=true;


int max_diff; // maximal difference for positions of alignment traces
// (only used for ends of arcs)
int max_diff_am; //maximal difference between two arc ends, -1 is off

std::string max_diff_pw_alignment; // pairwise reference alignment for max-diff heuristic, separator &
std::string max_diff_alignment_file; // reference alignment for max-diff heuristic, name of clustalw format file

bool opt_max_diff_relax=false;

// only consider arc matchs where
//   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0) !!! changed due to TraceController/reference alignment
//   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)

int exclusion_score; // Score contribution per exclusion
// set to zero for unrestricted structure locality

// double prob_exp_f(int seqlen) {return 1.0/(2*seqlen);} // expected probability of a base pair (null-model)

double exp_prob;
bool opt_exp_prob;

int output_width;

// ------------------------------------------------------------
// File arguments
std::string file1;
std::string file2;

std::string clustal_out;
bool opt_clustal_out;

std::string pp_out;
bool opt_pp_out;

// ------------------------------------------------------------
//
// Options
//
#include <LocARNA/options.hh>


bool opt_help;
bool opt_version;
bool opt_verbose;
bool opt_local_output;
bool opt_pos_output;

bool opt_write_structure;


bool opt_stacking;

std::string ribosum_file;
bool use_ribosum;

bool opt_probcons_file;
std::string probcons_file;

bool opt_mea_alignment;

bool opt_write_matchprobs;
bool opt_read_matchprobs;
std::string matchprobs_file;

bool opt_write_arcmatch_scores;

bool opt_read_arcmatch_scores;
bool opt_read_arcmatch_probs;
std::string arcmatch_scores_file;

int match_prob_method;

double min_am_prob; // only matched arc-pair with a probability of at least min_am_prob are taken into account
double min_bm_prob; // only matched base-pair with a probability of at least min_bm_prob are taken into account

bool opt_subopt; // find suboptimal solution (either k-best or all solutions better than a threshold)

int kbest_k;
int subopt_threshold;

std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;

int pf_struct_weight;

bool opt_mea_gapcost;
int mea_alpha;
int mea_beta;
int mea_gamma;
int probability_scale;

bool opt_eval;

bool opt_normalized;
int normalized_L;

option_def my_options[] = {
    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","Help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"match",'m',0,O_ARG_INT,&match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&indel_opening_score,"-500","score","Indel opening score"},
    {"struct-weight",'s',0,O_ARG_INT,&struct_weight,"200","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&opt_exp_prob,O_ARG_DOUBLE,&exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&tau_factor,"0","factor","Tau factor in percent"},
    {"exclusion",'E',0,O_ARG_INT,&exclusion_score,"0","score","Exclusion weight"},
    {"stacking",0,&opt_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use stacking terms (needs stack-probs by RNAfold -p2)"},   

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Type of locality"},

    {"struct-local",0,0,O_ARG_BOOL,&struct_local,"false","bool","Structure local"},
    {"sequ-local",0,0,O_ARG_BOOL,&sequ_local,"false","bool","Sequence local"},
    {"free-endgaps",0,0,O_ARG_STRING,&free_endgaps,"----","spec","Whether and which end gaps are free. order: L1,R1,L2,R2"},
    {"normalized",0,&opt_normalized,O_ARG_INT,&normalized_L,"0","L","Normalized local alignment with parameter L"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling output"},

    {"width",'w',0,O_ARG_INT,&output_width,"120","columns","Output width"},
    {"clustal",0,&opt_clustal_out,O_ARG_STRING,&clustal_out,O_NODEFAULT,"file","Clustal output"},
    {"pp",0,&opt_pp_out,O_ARG_STRING,&pp_out,O_NODEFAULT,"file","PP output"},
    {"local-output",'L',&opt_local_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment"},
    {"pos-output",'P',&opt_pos_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment positions"},
    {"write-structure",0,&opt_write_structure,O_NO_ARG,0,O_NODEFAULT,"","Write guidance structure in output"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},
    {"max-diff-aln",0,0,O_ARG_STRING,&max_diff_alignment_file,"","aln file","Maximal difference relative to given alignment (file in clustalw format))"},
    {"max-diff-pw-aln",0,0,O_ARG_STRING,&max_diff_pw_alignment,"","alignment","Maximal difference relative to given alignment (string, delim=&)"},
    {"max-diff-relax",0,&opt_max_diff_relax,O_NO_ARG,0,O_NODEFAULT,"","Relax deviation constraints in multiple aligmnent"},
    {"min-am-prob",'a',0,O_ARG_DOUBLE,&min_am_prob,"0.0005","amprob","Minimal Arc-match probability"},
    {"min-bm-prob",'b',0,O_ARG_DOUBLE,&min_bm_prob,"0.0005","bmprob","Minimal Base-match probability"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Special sauce options"},
    {"kbest",0,&opt_subopt,O_ARG_INT,&kbest_k,"-1","k","Enumerate k-best alignments"},
    {"better",0,&opt_subopt,O_ARG_INT,&subopt_threshold,"-1000000","t","Enumerate alignments better threshold t"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Options for controlling MEA score"},

    {"mea-alignment",0,&opt_mea_alignment,O_NO_ARG,0,O_NODEFAULT,"","Do MEA alignment"},
    {"probcons-file",0,&opt_probcons_file,O_ARG_STRING,&probcons_file,O_NODEFAULT,"file","Probcons parameter file"},

    {"match-prob-method",0,0,O_ARG_INT,&match_prob_method,"0","int","Method for computation of match probs"},
    {"temperature",0,0,O_ARG_INT,&temperature,"150","int","Temperature for PF-computation"},
    {"pf-struct-weight",0,0,O_ARG_INT,&pf_struct_weight,"200","weight","Structure weight in PF-computation"},

    {"mea-gapcost",0,&opt_mea_gapcost,O_NO_ARG,0,O_NODEFAULT,"","Use gap cost in mea alignment"},   
    {"mea-alpha",0,0,O_ARG_INT,&mea_alpha,"0","weight","Weight alpha for MEA"},
    {"mea-beta",0,0,O_ARG_INT,&mea_beta,"200","weight","Weight beta for MEA"},
    {"mea-gamma",0,0,O_ARG_INT,&mea_gamma,"100","weight","Weight gamma for MEA"},
    {"probability-scale",0,0,O_ARG_INT,&probability_scale,"10000","scale","Scale for probabilities/resolution of mea score"},

    {"write-match-probs",0,&opt_write_matchprobs,O_ARG_STRING,&matchprobs_file,O_NODEFAULT,"file","Write match probs to file (don't align!)"},
    {"read-match-probs",0,&opt_read_matchprobs,O_ARG_STRING,&matchprobs_file,O_NODEFAULT,"file","Read match probabilities from file"},

    {"write-arcmatch-scores",0,&opt_write_arcmatch_scores,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Write arcmatch scores (don't align!)"},
    {"read-arcmatch-scores",0,&opt_read_arcmatch_scores,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch scores"},
    {"read-arcmatch-probs",0,&opt_read_arcmatch_probs,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch probabilities (weight by mea_beta/100)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs"},
    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A"},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B"},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Mode of operation"},
    {"eval",0,&opt_eval,O_NO_ARG,0,O_NODEFAULT,"","Turn on evaluation mode"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&file1,O_NODEFAULT,"file 1","Basepairs input file 1 (alignment in eval mode)"},
    {"",0,0,O_ARG_STRING,&file2,O_NODEFAULT,"file 2","Basepairs input file 2 (dp dir in eval mode)"},
    {"",0,0,0,0,O_NODEFAULT,"",""}
};


// ------------------------------------------------------------



// ------------------------------------------------------------
// MAIN

int
main(int argc, char **argv) {

    typedef std::vector<int>::size_type size_type;

    // ------------------------------------------------------------
    // Process options

    bool process_success=process_options(argc,argv,my_options);

    if (opt_help) {
	cout << "locarna - a tool for pairwise (global and local) alignment of RNA."<<endl<<endl;
	
	cout << VERSION_STRING<<endl<<endl;

	print_help(argv[0],my_options);

	cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (opt_version || opt_verbose) {
	cout << "locarna ("<< VERSION_STRING<<")"<<endl;
	if (opt_version) exit(0); else cout <<endl;
    }

    if (!process_success) {
      std::cerr << "ERROR --- "
		<<O_error_msg<<std::endl;
      printf("USAGE: ");
      print_usage(argv[0],my_options);
      printf("\n");
      exit(-1);
    }
    
    if (opt_verbose) {
      print_options(my_options);
    }
    
    // ------------------------------------------------------------
    // parameter consistency
    if (opt_read_arcmatch_scores && opt_read_arcmatch_probs) {
	std::cerr << "You cannot specify arc match score and probabilities file simultaneously."<<std::endl;
	exit(-1);
    }
    
    if (probability_scale<=0) {
	std::cerr << "Probability scale must be greater 0."<<std::endl;
	exit(-1);
    }
    
    if (struct_weight<0) {
	std::cerr << "Structure weight must be greater equal 0."<<std::endl;
	exit(-1);
    }

    // ----------------------------------------
    // temporarily turn off stacking unless background prob is set
    //
    if (opt_stacking && !opt_exp_prob) {
	std::cerr << "WARNING: stacking turned off. "
		  << "Stacking requires setting a background probability "
		  << "explicitely (option --exp-prob)." << std::endl;
	opt_stacking=false;
    }


    // ----------------------------------------  
    // Ribosum matrix
    //
    std::auto_ptr<RibosumFreq> ribosum(NULL);
	
    if (use_ribosum) {
	if (ribosum_file == "RIBOSUM85_60") {
	    if (opt_verbose) {
		std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	    ribosum = auto_ptr<RibosumFreq>(new Ribosum85_60);
	} else {
	    ribosum = auto_ptr<RibosumFreq>(new RibosumFreq(ribosum_file));
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
    ScoringParams scoring_params(match_score,
				 mismatch_score,
				 // In true mea alignment gaps are only 
				 // scored for computing base match probs.
				 // Consequently, we set the indel and indel opening cost to 0
				 // for the case of mea alignment!
				 (opt_mea_alignment && !opt_mea_gapcost)?0:indel_score * (opt_mea_gapcost?probability_scale/100:1),
				 (opt_mea_alignment && !opt_mea_gapcost)?0:indel_opening_score * (opt_mea_gapcost?probability_scale/100:1),
				 ribosum.get(),
				 struct_weight,
				 tau_factor,
				 exclusion_score,
				 opt_exp_prob?exp_prob:-1,
				 temperature,
				 opt_stacking,
				 opt_mea_alignment,
				 mea_alpha,
				 mea_beta,
				 mea_gamma,
				 probability_scale
				 );

    
    // ----------------------------------------
    // CHOOSE MODE OF OPERATION
    //
    // whether to do evaluation
    //
    
    if (opt_eval) {
	std::cout <<"Evaluation Mode"<<std::endl;
	std::cout <<"==============="<<std::endl;
	
	Evaluator e(file1,file2,&scoring_params);
	std::cout <<"SCORE: "<< e.eval() << std::endl;
	
	std::exit(0);
    }
    
    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
    
    RnaData rnadataA(file1,opt_stacking);
    RnaData rnadataB(file2,opt_stacking);

    Sequence seqA=rnadataA.get_sequence();
    Sequence seqB=rnadataB.get_sequence();
    
    size_type lenA=seqA.length();
    size_type lenB=seqB.length();

    // --------------------
    // handle max_diff restriction  
    
    // missing: proper error handling in case that lenA, lenB, and max_diff_pw_alignment/max_diff_alignment_file are incompatible 
    
    // do inconsistency checking for max_diff_pw_alignment and max_diff_alignment_file
    //
    if (max_diff_pw_alignment!="" && max_diff_alignment_file!="") {
	std::cerr <<"Cannot simultaneously use both options --max-diff-pw-alignemnt and --max-diff-alignment-file."<<std::endl;
	exit(-1);
    }

    // construct TraceController and check inconsistency for with multiplicity of sequences
    //

    MultipleAlignment *multiple_ref_alignment=NULL;
    
    if (max_diff_alignment_file!="") {
	multiple_ref_alignment = new MultipleAlignment(max_diff_alignment_file);
    } else if (max_diff_pw_alignment!="") {
	if ( seqA.get_rows()!=1 || seqB.get_rows()!=1 ) {
	    std::cerr << "Cannot use --max-diff-pw-alignemnt for aligning of alignments." << std::endl;
	    exit(-1);
	}
	
	multiple_ref_alignment = new MultipleAlignment(seqA.names()[0],seqB.names()[0],max_diff_pw_alignment);
    }

    // if (multiple_ref_alignment) {
    // 	std::cout<<"Reference aligment:"<<std::endl;
    // 	multiple_ref_alignment->print_debug(std::cout);
    // 	std::cout << std::flush;
    // }
    
    TraceController trace_controller(seqA,seqB,multiple_ref_alignment,max_diff,opt_max_diff_relax);
    
    
    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    std::string seqCA = seq_constraints_A;
    std::string seqCB = seq_constraints_B;

    if (!opt_ignore_constraints) {
	if ( seqCA=="" ) seqCA = rnadataA.get_seq_constraints();
	if ( seqCB=="" ) seqCB = rnadataB.get_seq_constraints();
    }

    AnchorConstraints seq_constraints(seqA.length(),seqCA,
				      seqB.length(),seqCB);
    
    if (opt_verbose) {
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
    if (opt_read_arcmatch_scores || opt_read_arcmatch_probs) {
	if (opt_verbose) {
	    std::cout << "Read arcmatch scores from file " << arcmatch_scores_file << "." <<std::endl;
	}
	arc_matches = new ArcMatches(seqA,
				     seqB,
				     arcmatch_scores_file,
				     opt_read_arcmatch_probs?((mea_beta*probability_scale)/100):-1,
				     (max_diff_am!=-1)?(size_type)max_diff_am:std::max(lenA,lenB),
				     trace_controller,
				     seq_constraints
				     );
    } else {
	// initialize from RnaData
	arc_matches = new ArcMatches(rnadataA,
				     rnadataB,
				     min_prob,
				     (max_diff_am!=-1)?(size_type)max_diff_am:std::max(lenA,lenB),
				     trace_controller,
				     seq_constraints
				     );
    }
    
    BasePairs bpsA = arc_matches->get_base_pairsA();
    BasePairs bpsB = arc_matches->get_base_pairsB();
    
    // ----------------------------------------
    // report on input in verbose mode
    if (opt_verbose) {
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

    if (opt_read_matchprobs && !opt_mea_alignment) {
	std::cerr << "Warning: opt_read_matchprobs ignored for non-mea alignment.\n"; 
    }

    if (opt_write_matchprobs || opt_mea_alignment) {
	match_probs = new MatchProbs;

	if (!use_ribosum) {
	    std::cerr << "WARNING: Ribosum scoring used for mea_alignment and computing matchprobs."<<std::endl;
	}

	if (opt_read_matchprobs) {
	    match_probs->read_sparse(matchprobs_file,seqA.length(),seqB.length());
	} else {
	    if (match_prob_method==1) {
		if (!opt_probcons_file) {
		    std::cerr << "Probcons parameter file required for pairHMM-style computation"
			      <<" of basematch probabilities."<<std::endl;
		    print_usage(argv[0],my_options);
		    std::cerr << std::endl;
		    exit(-1);
		}
		if (opt_verbose) {
		    std::cout << "Compute match probabilities using pairHMM."<<std::endl; 
		}

		match_probs->pairHMM_probs(seqA,seqB,probcons_file);
	    } else {
		bool sl=sequ_local;
		if (match_prob_method==2) sl=true;
		if (match_prob_method==3) sl=false;

		if (opt_verbose) {
		    std::cout << "Compute match probabilities using PF sequence alignment."<<std::endl; 
		}

		match_probs->pf_probs(rnadataA,rnadataB,
				      ribosum->get_basematch_scores(),
				      ribosum->alphabet(),
				      indel_opening_score/100.0,
				      indel_score/100.0,
				      pf_struct_weight/100.0,
				      temperature/100.0,
				      sl);
	    }
	}

	if (opt_write_matchprobs) {
	    if (opt_verbose) {
		std::cout << "Write match probabilities to file "<<matchprobs_file<<"."<<std::endl; 
	    }

	    match_probs->write_sparse(matchprobs_file,1.0/probability_scale);
	    if (!opt_write_arcmatch_scores) exit(0); // else we exit there!
	}
    }
   

    // ----------------------------------------
    // construct scoring
   
    Scoring scoring(seqA,seqB,arc_matches,match_probs,&scoring_params);    

    if (opt_write_arcmatch_scores) {
	if (opt_verbose) {
	    std::cout << "Write arcmatch scores to file "<< arcmatch_scores_file<<" and exit."<<std::endl;
	}
	arc_matches->write_arcmatch_scores(arcmatch_scores_file,scoring);
	exit(0);
    }
        

    // ------------------------------------------------------------
    // Computation of the alignment score
    //

    // parameter for the alignment
    AlignerParams aligner_params(no_lonely_pairs,
				 struct_local,
				 sequ_local,
				 free_endgaps,
				 trace_controller,
				 max_diff_am,
				 min_am_prob,
				 min_bm_prob,
				 opt_stacking,
				 seq_constraints
				 );
    
    // initialize aligner object, which does the alignment computation
    Aligner aligner(seqA,seqB,*arc_matches,&aligner_params,&scoring);

    
    // enumerate suboptimal alignments (using interval splitting)
    if (opt_subopt) {
	aligner.suboptimal(kbest_k,
			   subopt_threshold,
			   opt_normalized,
			   normalized_L,
			   output_width,
			   opt_verbose,
			   opt_local_output,
			   opt_pos_output,
			   opt_write_structure
			   );
	exit(0);
    }
    
    infty_score_t score;

    // if option --normalized <L> is given, then do normalized local alignemnt
    if (opt_normalized) {
	
	// do some option consistency checks and output errors
	if (struct_local) {
	    std::cerr 
		<< "ERROR: Normalized structure local alignment not supported."
		<<std::endl
		<< "LocARNA ignores struct_local option."<<std::endl;
	    exit(-1);
	}
	if (!sequ_local) { // important: in the Aligner class, we rely on this
	    std::cerr 
		<< "ERROR: Normalized alignment requires option --sequ_local."<<std::endl;
	    exit(-1);
	}
	
	score = aligner.normalized_align(normalized_L,opt_verbose);
	
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
    if ((!opt_normalized) && DO_TRACE) {
	    
	if (opt_verbose) {
	    std::cout << "Traceback."<<std::endl;
	}
	
	aligner.trace();
	
	// for debugging:
	// aligner.get_alignment().write_debug(std::cout);
    }
    
    if (opt_normalized || DO_TRACE) { // if we did a trace (one way or
				      // the other)

	aligner.get_alignment().write(std::cout, 
				      output_width,
				      score,
				      opt_local_output,
				      opt_pos_output,
				      opt_write_structure
				      );
	std::cout<<endl;
	
	// test MultipleAlignment
	if (opt_verbose) {
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
	if (opt_clustal_out) {
	    ofstream out(clustal_out.c_str());
	    if (out.good()) {
		aligner.get_alignment().write_clustal(out, output_width,
						      score,
						      opt_local_output,
						      opt_pos_output,
						      true,
						      opt_write_structure
						      );
	    } else {
		cerr << "Cannot write to "<<clustal_out<<endl<<"! Exit.";
		exit(-1);
	    }
	}
	if (opt_pp_out) {
	    ofstream out(pp_out.c_str());
	    if (out.good()) {
		aligner.get_alignment().
		    write_pp(out,bpsA,bpsB,scoring,seq_constraints,output_width);
	    } else {
		cerr << "Cannot write to "<<pp_out<<endl<<"! Exit.";
		exit(-1);
	    }
	}
    }
    
    // ----------------------------------------
    // DONE
    exit(0);
}
