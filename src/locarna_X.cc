/**********************************************************************
 *
 * LocARNA X: LOCal Alignment of RNA eXact:
 * fast structure local exact matching
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de> 
 *               2005-2009
 *
 **********************************************************************/

// compile only when libRNA is available for linking
#ifdef HAVE_LIBRNA


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include<limits>
#include <stdio.h>
//#include <math.h>
// for getrusage()
#include <sys/resource.h>
#include <sys/types.h>
// for gettimeofday()
#include <sys/time.h>
// for setprecision
#include <iomanip>

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

bool no_lonely_pairs=false; // no lonely pairs option (currently not supported)

int max_diff; // maximal difference for positions of alignment traces
// (only used for ends of arcs)
int max_diff_am; //maximal difference between two arc ends, -1 is off

// only consider arc matchs where
//   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0)
//   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)
// except when there is no additional computation of M matrices necessary,
// this occurs if arcs are left-incident with larger arcs where 1 and 2 hold

int EPM_min_size; //minimum size for Exact Pattern Matches
double prob_unpaired_in_loop_threshold; // threshold for prob_unpaired_in_loop
double prob_unpaired_in_F_threshold;
//double prob_basepair_external_threshold;
int alpha_1; //parameter for sequential score
int alpha_2; //parameter for structural score
int alpha_3; //parameter for stacking score
int easier_scoring_par;
int difference_to_opt_score;
int min_subopt_score;
double subopt_range;
int am_threshold;
double coverage_cutoff;
std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;


// ------------------------------------------------------------
// File arguments
std::string file1;
std::string file2;
string psFile1;
string psFile2;
string locarna_output;
string clustal_output;
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
//bool opt_locarna_output;
bool opt_postscript_output;
bool opt_suboptimal;

option_def my_options[] = {
    {"min-prob",'P',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},

    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A."},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B."},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    
    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"stacking",'S',&opt_stacking,O_NO_ARG,0,O_NODEFAULT,"stacking","Use stacking terms (needs stack-probs by RNAfold -p2)"},
    {"EPM_minimum_size",'s',0,O_ARG_INT,&EPM_min_size,"2","min_size","User-defined minimum size for Exact Pattern Matches (chaining only)"},
    {"prob_unpaired_in_loop_threshold",'p',0,O_ARG_DOUBLE,&prob_unpaired_in_loop_threshold,"0.001","threshold","Threshold for prob_unpaired_in_loop"},
    {"prob_unpaired_in_F_threshold",0,0,O_ARG_DOUBLE,&prob_unpaired_in_F_threshold,"0.1","threshold","Threshold for prob_unpaired_in_F"},
    //{"prob_basepair_external_threshold",0,0,O_ARG_DOUBLE,&prob_basepair_external_threshold,"0.00001","threshold","Threshold for prob_basepair_external"},
    {"alpha_1",0,0,O_ARG_INT,&alpha_1,"1","alpha_1","Parameter for sequential score"},
    {"alpha_2",0,0,O_ARG_INT,&alpha_2,"1","alpha_2","Parameter for structural score"},
    {"alpha_3",0,0,O_ARG_INT,&alpha_3,"1","alpha_3","Parameter for stacking score, 0 means no stacking contribution"},
    {"suboptimal",0,&opt_suboptimal,O_NO_ARG,0,O_NODEFAULT,"suboptimal_traceback","Use a suboptimal traceback for the computation of the exact pattern matchings"},
    {"difference_to_optimal_score",0,0,O_ARG_INT,&difference_to_opt_score,"10","threshold","Threshold for suboptimal traceback"},
    {"min_subopt_score",0,0,O_ARG_INT,&min_subopt_score,"3","min","Minimal suboptimal score"},
    {"am-threshold",0,0,O_ARG_INT,&am_threshold,"3","am","Minimal arcmatch score in F matrix"},
    {"suboptimal-range",0,0,O_ARG_DOUBLE,&subopt_range,"0.0","range","trace EPMs within that range of best EPM score"},
    {"coverage-cutoff",0,0,O_ARG_DOUBLE,&coverage_cutoff,"0.5","cov","Skip chaining if best EPM has larger coverage on shortest seq"},
    {"easier_scoring_par",'e',0,O_ARG_INT,&easier_scoring_par,"0","alpha","use only sequential and a constant structural score alpha (easier_scoring_par) for each matched base of a basepair"},
    
    {"",0,0,O_ARG_STRING,&file1,O_NODEFAULT,"file 1","Basepairs input file 1 (alignment in eval mode)"},
    {"",0,0,O_ARG_STRING,&file2,O_NODEFAULT,"file 2","Basepairs input file 2 (dp dir in eval mode)"},
    {"PS_file1",'a',0,O_ARG_STRING,&psFile1,"","psFile1","Postscript output file for sequence A"},
    {"PS_file2",'b',0,O_ARG_STRING,&psFile2,"","psFile2","Postscript output file for sequence B"},
    {"output-ps", 0,&opt_postscript_output,O_NO_ARG,0,O_NODEFAULT,"","Output best EPM chain as colored postscript"},
//    {"output-locarna", 'o',&opt_locarna_output,O_NO_ARG,0,O_NODEFAULT,"","Output fasta file with anchor constraints for locarna"},
    {"output-locarna",'o',0,O_ARG_STRING,&locarna_output,"locarna_constraints_input.txt","constraintsFile","Fasta file with anchor constraints for locarna"},
    {"output-clustal",0,0,O_ARG_STRING,&clustal_output,"","filename","Write file with chain as alignment in clustalw format"},
    {"",0,0,0,0,O_NODEFAULT,"",""}
};


// ------------------------------------------------------------
// ------------------------------------------------------------
// MAIN

int
main(int argc, char **argv) {

	struct timeval tp;
	struct rusage ruse;

	gettimeofday( &tp, NULL );
	double start = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

	getrusage( RUSAGE_SELF, &ruse );
	double startR = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

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

    gettimeofday( &tp, NULL );
    double start_preproc = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
    getrusage( RUSAGE_SELF, &ruse );
    double startR_preproc = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
    
    gettimeofday( &tp, NULL );
    double start_fold = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
    getrusage( RUSAGE_SELF, &ruse );
    double start_foldR = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

    PFoldParams params(no_lonely_pairs,opt_stacking);

    RnaData rnadataA(file1,true,opt_stacking,true);
    if (!rnadataA.pairProbsAvailable() || !rnadataA.inLoopProbsAvailable()) {
	rnadataA.computeEnsembleProbs(params,true);
    }

    RnaData rnadataB(file2,true,opt_stacking,true);
    if (!rnadataB.pairProbsAvailable() || !rnadataB.inLoopProbsAvailable()) {
	rnadataB.computeEnsembleProbs(params,true);
    }

    Sequence seqA=rnadataA.get_sequence();
    Sequence seqB=rnadataB.get_sequence();
    
    MultipleAlignment maA(seqA);
    MultipleAlignment maB(seqB);

    gettimeofday( &tp, NULL );
    double end_fold = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
    getrusage( RUSAGE_SELF, &ruse );
    double end_foldR = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

    cout << "time_wall McCaskill_all = " << setprecision(3) << end_fold - start_fold << " sec" << endl;
    cout << "time_cpu McCaskill_all = " << setprecision(3) << end_foldR - start_foldR << " sec" << endl;

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
				 (max_diff_am!=-1)?(size_type)max_diff_am:std::max(seqA.length(),seqB.length()),
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
    
    //time_t start_mapping = time (NULL);
    Mapping mappingA(bpsA,
    		rnadataA,
    		prob_unpaired_in_loop_threshold,
    		//prob_unpaired_external_threshold,
    		//prob_basepair_external_threshold
    		prob_unpaired_in_F_threshold
    		);

    Mapping mappingB(bpsB,
    		rnadataB,
    		prob_unpaired_in_loop_threshold,
    		//prob_unpaired_external_threshold,
    		//prob_basepair_external_threshold
    		prob_unpaired_in_F_threshold
    		);
    //time_t stop_mapping = time (NULL);
    //cout << "time for mapping: " << stop_mapping - start_mapping << "sec " << endl;
    
    //string sequenceA= MultipleAlignment(seqA).seqentry(0).seq().to_string();
    //string sequenceB= MultipleAlignment(seqB).seqentry(0).seq().to_string();

    PatternPairMap myEPMs;

    ExactMatcher em(seqA,
		    seqB,
		    *arc_matches,
		    mappingA,
		    mappingB,
		    myEPMs,
		    alpha_1,
		    alpha_2,
		    alpha_3,
		    difference_to_opt_score,
		    min_subopt_score,
		    easier_scoring_par,
		    subopt_range,
		    am_threshold,
		    coverage_cutoff
		    );


    gettimeofday( &tp, NULL );
    double end_preproc = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
    getrusage( RUSAGE_SELF, &ruse );
    double endR_preproc = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

    cout << endl << "time_wall preprocessing = " << setprecision(3) << end_preproc - start_preproc << " sec" << endl;
    cout << "time_cpu preprocessing = " << setprecision(3) << endR_preproc - startR_preproc << " sec" << endl << endl;

    time_t start_computeMatrices = time (NULL);
    //compute matrices for finding best and enumerating all matchings
    em.compute_matrices();
    time_t stop_computeMatrices = time (NULL);
    cout << "time for computing EPM matrices : " << stop_computeMatrices - start_computeMatrices << "sec " << endl;
    
    
    // ------------------------------------------------------------
    // Traceback
    //
    const bool DO_TRACE=true;
    if (DO_TRACE) {
	
//	if (opt_verbose) {
//	    std::cout << "Traceback."<<std::endl;
//	}
	
	gettimeofday( &tp, NULL );
	double start_trace = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
	getrusage( RUSAGE_SELF, &ruse );
	double startR_trace = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;


	 if(opt_suboptimal){
		 cout << endl << "start suboptimal traceback..." << endl;
		 em.compute_EPMs_suboptimal();
	 }

	 else{
		 cout << endl << "start heuristic traceback..." << endl;
		 em.compute_EPMs_heuristic();
	 }

	 gettimeofday( &tp, NULL );
	 double end_trace = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
	 getrusage( RUSAGE_SELF, &ruse );
	 double endR_trace = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
	 cout << "time_wall traceback = " << setprecision(3) << end_trace - start_trace << " sec" << endl;
	 cout << "time_cpu traceback = " << setprecision(3) << endR_trace - startR_trace << " sec" << endl << endl;

   }
	 // chaining
	 //time_t start_chaining = time (	epm.sort_patVec();NULL);

    	 cout << "Start chaining..." << endl;
	 gettimeofday( &tp, NULL );
	 double start_chain = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

	 getrusage( RUSAGE_SELF, &ruse );
	 double startR_chain = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

	 PatternPairMap myLCSEPM;
	 LCSEPM myChaining(seqA, seqB, myEPMs, myLCSEPM, EPM_min_size);

	 //begin chaining algorithm
	 myChaining.calculateLCSEPM();

	 //time_t stop_chaining = time (NULL);

	 gettimeofday( &tp, NULL );
	 double end_chain = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

	 getrusage( RUSAGE_SELF, &ruse );
	 double endR_chain = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
	 cout << endl<< "time_wall chaining = " << end_chain - start_chain << " sec" << endl;
	 cout << "time_cpu chaining = "  << endR_chain - startR_chain << " sec" << endl << endl;

	 //
//	 cout << "time for chaining : " << stop_chaining - start_chaining << "sec " << endl;

	 //output chained EPMs to PS files
	 if(opt_postscript_output){
	//	 time_t start_ps = time (NULL);
	     assert(seqA.rows()==0);
	     assert(seqB.rows()==0);

		 if (opt_verbose) { cout << "write EPM chain as colored postscripts..." << endl;}
		 if (psFile1.size()==0){psFile1 = seqA.names()[0]+"_EPMs.ps";}
		 if (psFile2.size()==0){psFile2 = seqB.names()[0]+"_EPMs.ps";}

		 myChaining.MapToPS(maA.consensus_sequence(), maB.consensus_sequence(), myLCSEPM, psFile1,psFile2);
//		 time_t stop_ps = time (NULL);
		// cout << "time for map to ps : " << stop_ps - start_ps << "sec " << endl;
	 }

	 //if(opt_locarna_output){
		 if (opt_verbose) { cout << "write locarna anchor constraints..." << endl;}
		 myChaining.output_locarna(maA.consensus_sequence(), maB.consensus_sequence(), locarna_output);
	 //}

	 if (clustal_output.size()>0){
		 if (opt_verbose) { cout << "write chain as clustal alignment..." << endl;}
		 myChaining.output_clustal(clustal_output);
	 }

    gettimeofday( &tp, NULL );
    double end = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

    getrusage( RUSAGE_SELF, &ruse );
    double endR = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
    cout << endl << "time_wall main = " << setprecision(3) << end - start << " sec" << endl;
    cout << "time_cpu main = " << setprecision(3) << endR - startR << " sec" << endl << endl;

    // ----------------------------------------
    // DONE
    delete arc_matches;
    cout << "... locarna_X finished!" << endl << endl;
    return 0;
}

#else // HAVE_LIBRNA
#include "iostream"
int
main() {
    std::cerr
	<< "Functionality of locarna_X is not available,"
	<< "since LocARNA was compiled without libRNA support." << std::endl
	<< "Requires recompilation with configure option --enable-librna." <<std::endl; 
    return -1;
}

#endif // HAVE_LIBRNA
