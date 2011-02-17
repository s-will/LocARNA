/**********************************************************************
 *
 * LocARNA X: LOCal Alignment of RNA eXact:
 * fast structure local exact matching
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de> 
 *               2005-2009
 *
 **********************************************************************/


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//#include <math.h>

#include <LocARNA/sequence.hh>
#include <LocARNA/basepairs.hh>

//#include <LocARNA/exact_matcher.hh>

#include <LocARNA/rna_data.hh>
#include <LocARNA/arc_matches.hh>
#include <LocARNA/match_probs.hh>

#include <LocARNA/anchor_constraints.hh>
#include <LocARNA/trace_controller.hh>

#include <LocARNA/exact_matcher.hh>

using namespace std;


const std::string
VERSION_STRING = (std::string)PACKAGE_STRING;

// ------------------------------------------------------------
// Parameter

double min_prob; // only pairs with a probability of at least min_prob are taken into account

// bool no_lonely_pairs; // no lonely pairs option

const bool DO_TRACE=true;

int max_diff; // maximal difference for positions of alignment traces
// (only used for ends of arcs)
int max_diff_am; //maximal difference between two arc ends, -1 is off

// only consider arc matchs where
//   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0)
//   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)
// except when there is no additional computation of M matrices necessary,
// this occurs if arcs are left-incident with larger arcs where 1 and 2 hold



std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;


// ------------------------------------------------------------
// File arguments
std::string file1;
std::string file2;

// ------------------------------------------------------------
//
// Options
//
#include <LocARNA/options.hh>

using namespace LocARNA;

bool opt_help;
bool opt_version;
bool opt_verbose;


option_def my_options[] = {
    {"min-prob",'p',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff-match",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},

    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A."},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B."},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    
    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    
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
	cout << VERSION_STRING<<endl;

	cout << "Copyright Sebastian Will, 2005-2009"<<endl<<endl;

	cout << "A tool for pairwise Local (and global) Alignment of RNA: Exact Local Matchings."<<endl<<endl;

	print_help(argv[0],my_options);

	cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (opt_version || opt_verbose) {
	cout << VERSION_STRING<<endl;
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

    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
    
    RnaData rnadataA(file1,false);
    RnaData rnadataB(file2,false);

    Sequence seqA=rnadataA.get_sequence();
    Sequence seqB=rnadataB.get_sequence();
    
    size_type lenA=seqA.length();
    size_type lenB=seqB.length();

    // --------------------
    // handle max_diff restriction
    
    TraceController trace_controller(seqA,seqB,NULL,max_diff);
    
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
    // Compute Exact Matchings
    //

    ExactMatcher em(seqA,seqB,*arc_matches);
    
    // compute matrices for finding best and enumerating all matchings
    infty_score_t score = em.compute_matchings();
    
    // ----------------------------------------
    // report score
    //
    std::cout << "Size of best exact matching: "<<score<<std::endl;
    
    
    // ------------------------------------------------------------
    // Traceback
    //
    if (DO_TRACE) {
	
	if (opt_verbose) {
	    std::cout << "Traceback."<<std::endl;
	}
	
	em.write_matchings(cout);
    }
    
    // ----------------------------------------
    // DONE
    exit(0);
}
