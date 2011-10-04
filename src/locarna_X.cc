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
#include <sstream>
#include<limits>
#include <stdio.h>
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

int EPM_threshold; //threshold for Exact Pattern Matches
int EPM_min_size; //minimum size for Exact Pattern Matches
double prob_unpaired_threshold; // threshold for prob_unpaired_in_loop
int alpha_1; //parameter for sequential score
int alpha_2; //parameter for structural score
int alpha_3; //parameter for stacking score
int easier_scoring_par;
int subopt_score;


std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;


// ------------------------------------------------------------
// File arguments
std::string file1;
std::string file2;
string psFile1;
string psFile2;

// ------------------------------------------------------------
//
// Options
//
#include <LocARNA/options.hh>

using namespace LocARNA;

bool opt_help;
bool opt_version;
bool opt_verbose;
bool opt_stacking;
bool opt_locarna_output;
bool opt_suboptimal;

option_def my_options[] = {
    {"min-prob",'P',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff-match",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},

    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A."},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B."},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    
    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"stacking",'S',&opt_stacking,O_NO_ARG,0,O_NODEFAULT,"stacking","Use stacking terms (needs stack-probs by RNAfold -p2)"},
    {"EPM_threshold",'t',0,O_ARG_INT,&EPM_threshold,"5","threshold","User-defined threshold for Exact Pattern Matches"},
    {"EPM_minimum_size",'s',0,O_ARG_INT,&EPM_min_size,"3","min_size","User-defined minimum size for Exact Pattern Matches"},
    {"prob_unpaired_threshold",'p',0,O_ARG_DOUBLE,&prob_unpaired_threshold,"0.001","threshold","Threshold for prob_unpaired_in_loop"},
    {"alpha_1",0,0,O_ARG_INT,&alpha_1,"1","alpha_1","Parameter for sequential score"},
    {"alpha_2",0,0,O_ARG_INT,&alpha_2,"1","alpha_2","Parameter for structural score"},
    {"alpha_3",0,0,O_ARG_INT,&alpha_3,"1","alpha_3","Parameter for stacking score, 0 means no stacking contribution"},
    {"suboptimal",0,&opt_suboptimal,O_NO_ARG,0,O_NODEFAULT,"suboptimal_traceback","Use a suboptimal traceback for the computation of the exact pattern matchings"},
    {"suboptimal_score",0,0,O_ARG_INT,&subopt_score,"5","alpha_1","Threshold for suboptimal traceback"},
    {"easier_scoring_par",'e',0,O_ARG_INT,&easier_scoring_par,"0","alpha","use only sequential and a constant structural score alpha (easier_scoring_par) for each matched base of a basepair"},
    
    {"",0,0,O_ARG_STRING,&file1,O_NODEFAULT,"file 1","Basepairs input file 1 (alignment in eval mode)"},
    {"",0,0,O_ARG_STRING,&file2,O_NODEFAULT,"file 2","Basepairs input file 2 (dp dir in eval mode)"},
    {"PS_file1",'a',0,O_ARG_STRING,&psFile1,"SequenceA","psFile1","Postscript output file for sequence A"},
    {"PS_file2",'b',0,O_ARG_STRING,&psFile2,"SequenceB","psFile2","Postscript output file for sequence B"},
    {"output_locarna", 'i',&opt_locarna_output,O_NO_ARG,0,O_NODEFAULT,"","Output with anchor constraints for locarna"},
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
    
    //if no stacking should be considered, set the parameter for stacking to 0
    if(!opt_stacking){
    	alpha_3 = 0;
    }

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
   
    time_t start_preprocessing = time (NULL);
    
    MultipleAlignment maA(file1);
    MultipleAlignment maB(file2);
    
    
    
    Sequence seqA(maA);
    Sequence seqB(maB);

    time_t start_RNAdataA = time (NULL);
    RnaData rnadataA(seqA,true,opt_stacking);
    time_t stop_RNAdataA = time (NULL);
    cout << "time for RNAdataA: " << stop_RNAdataA - start_RNAdataA << "sec " << endl;
    time_t start_RNAdataB = time (NULL);
    RnaData rnadataB(seqB,true,opt_stacking);
    time_t stop_RNAdataB = time (NULL);
    cout << "time for RNAdataB: " << stop_RNAdataB - start_RNAdataB << "sec " << endl;

   
    //Sequence seqA=rnadataA.get_sequence();
    //Sequence seqB=rnadataB.get_sequence();
    
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
    
    time_t start_mapping = time (NULL);
    Mapping mappingA(bpsA,rnadataA,prob_unpaired_threshold);
    Mapping mappingB(bpsB,rnadataB,prob_unpaired_threshold);
    time_t stop_mapping = time (NULL);
    cout << "time for mapping: " << stop_mapping - start_mapping << "sec " << endl;
    
    string sequenceA= MultipleAlignment(seqA).seqentry(0).seq().to_string();
    string sequenceB= MultipleAlignment(seqB).seqentry(0).seq().to_string();
    ExactMatcher em(seqA,
		    seqB,
		    *arc_matches,
		    mappingA,
		    mappingB,
		    EPM_threshold,
		    EPM_min_size,
		    alpha_1,
		    alpha_2,
		    alpha_3,
		    subopt_score,
		    easier_scoring_par,
		    sequenceA,
		    sequenceB,
		    psFile1,
		    psFile2
 		  );

    time_t stop_preprocessing = time (NULL);

    cout << "time for preprocessing: " << stop_preprocessing - start_preprocessing << "sec " << endl;

    time_t start_computeMatrices = time (NULL);
    //compute matrices for finding best and enumerating all matchings
    em.compute_matrices();
    time_t stop_computeMatrices = time (NULL);
    cout << "time for computing Matrices : " << stop_computeMatrices - start_computeMatrices << "sec " << endl;
    
    
    // ------------------------------------------------------------
    // Traceback
    //
   if (DO_TRACE) {
	
	if (opt_verbose) {
	    std::cout << "Traceback."<<std::endl;
	}
	
	 if(opt_suboptimal){
		 cout << "suboptimal " << endl;
		 em.compute_EPMs_suboptimal();

	 }

	 else{
		 cout << "heuristic " << endl;
		 time_t start_traceback = time (NULL);
		 em.compute_EPMs_heuristic();
		 time_t stop_traceback = time(NULL);
		 cout << "time for traceback : " << stop_traceback - start_traceback << "sec " << endl;
	 }


    }

    if(opt_locarna_output){
    	em.output_locarna();
    }
    // ----------------------------------------
    // DONE
    delete arc_matches;
    return 0;
}
