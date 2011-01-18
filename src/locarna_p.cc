/**********************************************************************
 *
 * LocARNA-P: LOCal Alignment of RNA - Partitition Function Variant
 *
 * Copyright (C) Sebastian Will <will(@@)informatik.uni-freiburg.de> 
 *               2005-2010
 *
 **********************************************************************/


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <math.h>

#include "sequence.hh"
#include "basepairs.hh"
#include "aligner_p.hh"

#include "rna_data.hh"
#include "arc_matches.hh"

#include "match_probs.hh"
#include "ribosum.hh"

#include "ribosum85_60.icc"

#include "anchor_constraints.hh"

#include "trace_controller.hh"
#include "multiple_alignment.hh"

using namespace std;

const std::string 
VERSION_STRING = (std::string)PACKAGE_STRING; 

// ------------------------------------------------------------
// Parameter

double min_prob; // only pairs with a probability of at least min_prob are taken into account

int mismatch_score;
int indel_score;
int indel_opening_score=0;
int temperature;
int struct_weight;

int tau_factor; // contribution of sequence similarity in an arc match (in percent)


bool struct_local; // allow exclusions for maximizing alignment of connected substructures 
bool sequ_local; // maximize alignment of subsequences

const bool DO_TRACE=true;


int max_diff; // maximal difference for positions of alignment traces
// (only used for ends of arcs)
int max_diff_am; //maximal difference between two arc ends, -1 is off

std::string max_diff_pw_alignment; // pairwise reference alignment for max-diff heuristic, separator &
std::string max_diff_alignment_file; // reference alignment for max-diff heuristic, name of clustalw format file

double min_am_prob; // only matched arc-pair with a probability of at least min_am_prob are taken into account
double min_bm_prob; // only matched base-pair with a probability of at least min_bm_prob are taken into account
int match_score;

// only consider arc matchs where
//   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0) !!!!! this changed due to TraceController
//   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)

int exclusion_score; // Score contribution per exclusion
// set to zero for unrestricted structure locality

//double prob_exp_f(int seqlen) {return 1.0/(2*seqlen);} // expected probability of a base pair (null-model)

double exp_prob;
bool opt_exp_prob;
int output_width=70;
bool opt_write_arcmatch_probs;
bool opt_write_basematch_probs;
// ------------------------------------------------------------
// File arguments

std::string arcmatch_probs_file;
std::string basematch_probs_file;

std::string bpsfile1;
std::string bpsfile2;

std::string clustal_out;
bool opt_clustal_out;

std::string pp_out;
bool opt_pp_out;

// ------------------------------------------------------------
//
// Options
//
#include "options.hh"


bool opt_help;
bool opt_version;
bool opt_verbose;
bool opt_local_output;
bool opt_pos_output;

std::string ribosum_file;
bool use_ribosum;

bool basematch_probs_include_arcmatch;

std::string fragment_match_probs;

// IMPORTANT: here use double (not pf_score_t), since option parser cannot interpret pf_score_t, when using long double
// later this will be converted to pf_score_t
double pf_scale;

option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"match",'m',0,O_ARG_INT,&match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&indel_opening_score,"-500","score","Indel opening score"},
    {"struct-weight",'s',0,O_ARG_INT,&struct_weight,"180","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&opt_exp_prob,O_ARG_DOUBLE,&exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&tau_factor,"0","factor","Tau factor in percent"},
    {"temperature",0,0,O_ARG_INT,&temperature,"150","int","Temperature for PF-computation"},
    
    {"pf-scale",0,0,O_ARG_DOUBLE,&pf_scale,"1.0","scale","Scaling of the partition function. Use in order to avoid overflow."},
    
    {"min-prob",'p',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"min-am-prob",'a',0,O_ARG_DOUBLE,&min_am_prob,"0.0005","amprob","Minimal Arc-match probability"},
    {"min-bm-prob",'b',0,O_ARG_DOUBLE,&min_bm_prob,"0.0005","bmprob","Minimal Base-match probability"},

    {"include-am-in-bm",0,&basematch_probs_include_arcmatch,O_NO_ARG,0,O_NODEFAULT,"","Include arc match cases in computation of base match probabilities"},

    // {"kbest",'k',0,O_ARG_INT,&kbest_k,"-1","k","Find k-best alignments"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling output"},
    {"write-arcmatch-probs",0,&opt_write_arcmatch_probs,O_ARG_STRING,&arcmatch_probs_file,O_NODEFAULT,"file","Write arcmatch probabilities"},
    {"write-basematch-probs",0,&opt_write_basematch_probs,O_ARG_STRING,&basematch_probs_file,O_NODEFAULT,"file","Write basematch probabilities"},
    {"width",'w',0,O_ARG_INT,&output_width,"120","columns","Output width"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    
    {"max-diff-aln",0,0,O_ARG_STRING,&max_diff_alignment_file,"","aln file","Maximal difference relative to given alignment (file in clustalw format))"},
    {"max-diff-pw-aln",0,0,O_ARG_STRING,&max_diff_pw_alignment,"","alignment","Maximal difference relative to given alignment (string, delim=&)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Computed probabilities"},
    {"fragment-match-probs",0,0,O_ARG_STRING,&fragment_match_probs,"","\"i j k l\"",
     "Requests probabilities for the match of fragments [i..j] and [k..l]. Accepts a ';' separated list of ranges."},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Standard options"},

    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&bpsfile1,O_NODEFAULT,"bps-file 1","Basepairs input file 1"},
    {"",0,0,O_ARG_STRING,&bpsfile2,O_NODEFAULT,"bps-file 2","Basepairs input file 2"},


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
	cout << VERSION_STRING<<endl;

	cout << "Copyright Sebastian Will, 2004-2009"<<endl<<endl;

	cout << "A tool for pairwise Local (and global) Alignment of RNA."<<endl<<endl;
	cout << "LocARNA-P: computes partitition function and edge probabilities of alignment."<<endl<<endl;
		
	print_help(argv[0],my_options);

	cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (opt_version || opt_verbose) {
	cout << VERSION_STRING<<endl;
	if (opt_version) exit(0); else cout <<endl;
    }

    if (!process_success) {
	puts(O_error_msg);
	printf("USAGE: ");
	print_usage(argv[0],my_options);
	printf("\n");
	exit(-1);
    }

    if (opt_verbose)
	print_options(my_options);


    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    RnaData rnadataA(bpsfile1);
    RnaData rnadataB(bpsfile2);

    Sequence seqA=rnadataA.get_sequence();
    Sequence seqB=rnadataB.get_sequence();

    size_type lenA=seqA.length();
    size_type lenB=seqB.length();

    AnchorConstraints seq_constraints(seqA.length(),"",seqB.length(),"");

    // --------------------
    // handle max_diff restriction
    
    // missing: proper error handling in case that lenA, lenB, and max_diff_alignment are incompatible 

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
    
    TraceController trace_controller(seqA,seqB,multiple_ref_alignment,max_diff);
    
    if (multiple_ref_alignment) {
	delete multiple_ref_alignment;
    }
    
    // ----------------------------------------
    // construct set of relevant arc matches
    //
    ArcMatches *arc_matches;
    
    // initialize from RnaData
    arc_matches = new ArcMatches(rnadataA,
				 rnadataB,
				 min_prob,
				 (max_diff_am!=-1)?(size_type)max_diff_am:std::max(lenA,lenB),
				 trace_controller,
				 seq_constraints
				 );
	
    BasePairs bpsA = arc_matches->get_base_pairsA();
    BasePairs bpsB = arc_matches->get_base_pairsB();
	
    // ----------------------------------------
    // report on input in verbose mode
    if (opt_verbose) {
	cout << "Sequence A: "<<endl;
	seqA.write(cout);
	cout<<" (Length:"<< seqA.length()<<", Basepairs:"<<bpsA.num_bps() << ")" <<endl;

	cout << "Sequence B: "<<endl;
	seqB.write(cout);
	cout<<" (Length:"<< seqB.length()<<", Basepairs:"<<bpsB.num_bps() << ")" <<endl;

	// cout << "Base Identity: "<<(seq_identity(seqA,seqB)*100)<<endl; 
    }
	
    // ----------------------------------------  
    // Ribosum matrix
    //
    RibosumFreq *ribosum=NULL;

    if (use_ribosum) {
	if (ribosum_file == "RIBOSUM85_60") {
	    if (opt_verbose) {
		std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	    ribosum = new Ribosum85_60();
	} else {
	    ribosum = new RibosumFreq(ribosum_file);
	}
    }

    //  MatchProbs *match_probs=0L;

    // ----------------------------------------
    // construct scoring
    if (opt_verbose) {
	std::cout << "Construct scoring."<<std::endl;
    }

	
    ScoringParams scoring_params(match_score,
				 mismatch_score,
				 indel_score,
				 indel_opening_score,
				 ribosum,
				 struct_weight,
				 tau_factor,
				 0, // exclusion score
				 opt_exp_prob?exp_prob:-1,
				 temperature,//temperature,
				 false,//opt_stacking,
				 false,//opt_mea_alignment,
				 0,//mea_alpha,
				 0,//mea_beta,
				 0,//mea_gamma,
				 0 //probability_scale
				 );
	
    Scoring scoring(seqA,seqB,arc_matches,0L,&scoring_params,true);
    
    // ------------------------------------------------------------
    // Computation of the alignment score
    //
    // parameter for the alignment
    AlignerParams aligner_params(false,//no_lonely_pairs,
				 false,//struct_local,
				 false,//sequ_local,
				 "",
				 trace_controller,
				 max_diff_am,
				 min_am_prob,
				 min_bm_prob,
				 false,//opt_stacking,
				 seq_constraints
				 );
    
    
    // initialize aligner object, which does the alignment computation
    AlignerP aligner(seqA,seqB,*arc_matches,&aligner_params,&scoring,pf_scale);

    if (opt_verbose) {
	std::cout << "Run inside algorithm."<<std::endl;
#       ifdef LARGE_PF
	std::cout << "Use large partition function type."<<std::endl;
#       endif
    }

    pf_score_t pf=aligner.align_inside();
    
    std::cout << "Partition function: "<<pf<<std::endl;
    
    if (opt_verbose) {
	std::cout << "Run outside algorithm."<<std::endl;
    }
	
    aligner.align_outside();

    if (opt_verbose) {
	std::cout << "Compute probabilities."<<std::endl;
    }
	
    aligner.compute_arcmatch_probabilities();

    
//     if (opt_verbose) {
// 	std::cout << "Arc match probabilities:"<<std::endl;
// 	aligner.write_arcmatch_probabilities(std::cout);
// 	std::cout << std::endl;
//     }
	
    if (opt_write_arcmatch_probs) {
	if (opt_verbose) {
	    std::cout << "Write Arc-match probabilities to file "<<arcmatch_probs_file<<"."<<std::endl; 
	}
	ofstream out(arcmatch_probs_file.c_str());
	if (out.good()) {
	    aligner.write_arcmatch_probabilities(out);
	} else {
	    cerr << "Cannot write to "<<arcmatch_probs_file<<"! Exit."<<endl;
	    exit(-1);
	}
    }
	
    
    aligner.compute_basematch_probabilities(basematch_probs_include_arcmatch);
	
//     if (opt_verbose) {
// 	std::cout << "Base match probabilities:"<<std::endl;
// 	aligner.write_basematch_probabilities(std::cout);
// 	std::cout << std::endl;
//     }


    if (opt_write_basematch_probs) {
	if (opt_verbose) {
	    std::cout << "Write Base-match probabilities to file "<<basematch_probs_file<<"."<<std::endl; 
	}
	ofstream out(basematch_probs_file.c_str());
	if (out.good()) {
	    aligner.write_basematch_probabilities(out);
	} else {
	    cerr << "Cannot write to "<<basematch_probs_file<<"! Exit."<<endl;
	    exit(-1);
	}
    }

    // ----------------------------------------
    // optionally, compute probabilities that certain ranges are matched
    //
    
    if (fragment_match_probs != "") {
	// parse the string fragment_match_probs
	// the string must have the format 
	// fragment_list ::= empty | pos pos pos pos {; pos pos pos pos}
	
	std::istringstream fragments(fragment_match_probs);
	
	size_type i;
	size_type j;
	size_type k;
	size_type l;

	while (! fragments.eof()) {
	    
	    if (!(fragments >> i >> j >> k >> l)) {
		std::cerr << "WARNING: expected 4 positions or end when parsing argument fragment-match-probs"
			  <<"\""<<fragment_match_probs<<"\""<<std::endl;
		break;
	    }
	    
	    std::cout <<"FRAGMENT_MATCH_PROB "<<i<<" "<<j<<" "<<k<<" "<<l<<" : ";
	    std::cout.flush();
	    std::cout << aligner.compute_fragment_match_prob(i,j,k,l) << std::endl;
	    
	    if (! fragments.eof()) {
		char sep;
		fragments >> sep;
		if (sep != ';') {
		    std::cerr << "WARNING: expected ';' or end when parsing argument fragment-match-probs"
			  <<std::endl;
		    break;
		}
	    }
	}
    }
    
    // ----------------------------------------
    // for debugging: print all arcmatch scores
    
    //for (ArcMatches::const_iterator it=arc_matches->begin(); it!=arc_matches->end(); ++it) {
    // std::cout << it->arcA() <<" " << it->arcB() <<" " << it->idx() << " " << scoring.exp_arcmatch(*it)<< std::endl;
    //}

    
    // DONE
    exit(0);
}
