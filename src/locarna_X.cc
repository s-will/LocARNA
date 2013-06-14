/**********************************************************************
 *
 * LocARNA X: LOCal Alignment of RNA eXact:
 * fast structure local exact matching
 *
 * Copyright (C) Sebastian Will <will(@)informatik.uni-freiburg.de> 
 *
 **********************************************************************/

// need to include config.h already here
// because of following #ifdef HAVE_LIBRNA
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

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

#include "LocARNA/sequence.hh"
#include "LocARNA/basepairs.hh"

//#include "LocARNA/exact_matcher.hh"

#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/match_probs.hh"

#include "LocARNA/anchor_constraints.hh"
#include "LocARNA/sequence_annotations.hh"
#include "LocARNA/trace_controller.hh"

#include "LocARNA/exact_matcher.hh"
#include "LocARNA/sparsification_mapper.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/global_stopwatch.hh"


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
double prob_basepair_in_loop_threshold; // threshold for prob_basepair_in_loop
double prob_unpaired_in_F_threshold;
//double prob_basepair_external_threshold;
int alpha_1; //parameter for sequential score
int alpha_2; //parameter for structural score
int alpha_3; //parameter for stacking score
LocARNA::score_t easier_scoring_par;
int difference_to_opt_score;
int min_subopt_score;
double subopt_range;
int am_threshold;
double coverage_cutoff;
std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;

bool chaining;


// ------------------------------------------------------------
// File arguments
std::string fileA;
std::string fileB;
string psFileA;
string psFileB;
string locarna_output;
string clustal_output;
string epm_list_output;
string chained_epm_list_output;

// ------------------------------------------------------------
//
// Options
//
#include "LocARNA/options.hh"

using namespace LocARNA;

bool opt_help;
bool opt_version;
bool opt_verbose;
bool opt_stacking;
//bool opt_locarna_output;
bool opt_postscript_output;
bool opt_suboptimal;

bool opt_stopwatch;

option_def my_options[] = {
    {"min-prob",'P',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment traces"},

    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    
    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"stacking",'S',&opt_stacking,O_NO_ARG,0,O_NODEFAULT,"stacking","Use stacking terms (needs stack-probs by RNAfold -p2)"},
    {"EPM_minimum_size",'s',0,O_ARG_INT,&EPM_min_size,"2","min_size","User-defined minimum size for Exact Pattern Matches (chaining only)"},
    {"prob_unpaired_in_loop_threshold",'p',0,O_ARG_DOUBLE,&prob_unpaired_in_loop_threshold,"0.001","threshold","Threshold for prob_unpaired_in_loop"},
    {"prob_basepair_in_loop_threshold",'q',0,O_ARG_DOUBLE,&prob_basepair_in_loop_threshold,"0.001","threshold","Threshold for prob_basepair_in_loop"},
    {"alpha_1",0,0,O_ARG_INT,&alpha_1,"1","alpha_1","Multiplier for sequential score"},
    {"alpha_2",0,0,O_ARG_INT,&alpha_2,"1","alpha_2","Multiplier for structural score"},
    {"alpha_3",0,0,O_ARG_INT,&alpha_3,"1","alpha_3","Multiplier for stacking score, 0 means no stacking contribution"},
    {"suboptimal",0,&opt_suboptimal,O_NO_ARG,0,O_NODEFAULT,"suboptimal_traceback","Use a suboptimal traceback for the computation of the exact pattern matchings"},
    {"difference_to_optimal_score",0,0,O_ARG_INT,&difference_to_opt_score,"10","threshold","Threshold for suboptimal traceback"},
    {"min_subopt_score",0,0,O_ARG_INT,&min_subopt_score,"3","min","Minimal suboptimal score"},
    {"am-threshold",0,0,O_ARG_INT,&am_threshold,"3","am","Minimal arcmatch score in F matrix"},
    {"suboptimal-range",0,0,O_ARG_DOUBLE,&subopt_range,"0.0","range","trace EPMs within that range of best EPM score"},
    {"coverage-cutoff",0,0,O_ARG_DOUBLE,&coverage_cutoff,"0.5","cov","Skip chaining if best EPM has larger coverage on shortest seq"},
    {"easier_scoring_par",'e',0,O_ARG_INT,&easier_scoring_par,"0","alpha","use only sequential and a constant structural score alpha (easier_scoring_par) for each matched base of a basepair"},
    
    {"stopwatch",0,&opt_stopwatch,O_NO_ARG,0,O_NODEFAULT,"","Print run time information."},

    {"",0,0,O_ARG_STRING,&fileA,O_NODEFAULT,"file A","input file A"},
    {"",0,0,O_ARG_STRING,&fileB,O_NODEFAULT,"file B","input file B"},
    {"PS_fileA",'a',0,O_ARG_STRING,&psFileA,"","psFileA","Postscript output file for sequence A"},
    {"PS_fileB",'b',0,O_ARG_STRING,&psFileB,"","psFileB","Postscript output file for sequence B"},
    {"output-ps", 0,&opt_postscript_output,O_NO_ARG,0,O_NODEFAULT,"","Output best EPM chain as colored postscript"},
//    {"output-locarna", 'o',&opt_locarna_output,O_NO_ARG,0,O_NODEFAULT,"","Output fasta file with anchor constraints for locarna"},
    {"output-locarna",'o',0,O_ARG_STRING,&locarna_output,"locarna_constraints_input.txt","constraintsFile","Fasta file with anchor constraints for locarna"},
    {"output-clustal",0,0,O_ARG_STRING,&clustal_output,"","filename","Write file with chain as alignment in clustalw format"},
    {"output-epm-list",0,0,O_ARG_STRING,&epm_list_output,"","epm list","A list of all found epms"},
    {"output-chained-epm-list",0,0,O_ARG_STRING,&chained_epm_list_output,"","chained epm list","A list of all EPMs that are present in the chain"},
    {"chaining",0,&chaining,O_NO_ARG,0,O_NODEFAULT,"chaining","use the chaining algorithm to find best overall chain"},
    {"",0,0,0,0,O_NODEFAULT,"",""}


};

// ------------------------------------------------------------
// ------------------------------------------------------------
// MAIN

int
main(int argc, char **argv) {
    stopwatch.start("total");

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

	cout << "(C) Sebastian Will"<<endl<<endl;

	cout << "A tool for pairwise Local (and global) Alignment of RNA: Exact Local Matchings."<<endl<<endl;

	print_help(argv[0],my_options);

	cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	return 0;
    }

    if (opt_version || opt_verbose) {
	cout << VERSION_STRING<<endl;
	if (opt_version) return 0; else cout <<endl;
    }

    if (!process_success) {
      std::cerr << "ERROR --- "
		<<O_error_msg<<std::endl;
      printf("USAGE: ");
      print_usage(argv[0],my_options);
      printf("\n");
      return -1;
    }
    
    if (opt_stopwatch) {
	stopwatch.set_print_on_exit(true);
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

    if(!chaining && chained_epm_list_output.size()>0){
    	cout << "Enable chaining in order to output chained epm list " << endl;
    	chaining = true;
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

    PFoldParams pfparams(no_lonely_pairs,opt_stacking);

    ExtRnaData *rna_dataA=0;
    try {
	rna_dataA = new ExtRnaData(fileA,
				   min_prob,
				   prob_basepair_in_loop_threshold,
				   prob_unpaired_in_loop_threshold,
				   pfparams);
    } catch (failure &f) {
	std::cerr << "ERROR: failed to read from file "<<fileA <<std::endl
		  << "       "<< f.what() <<std::endl;
	return -1;
    }
    
    ExtRnaData *rna_dataB=0;
    try {
	rna_dataB = new ExtRnaData(fileB,
				   min_prob,
				   prob_basepair_in_loop_threshold,
				   prob_unpaired_in_loop_threshold,
				   pfparams);
    } catch (failure &f) {
	std::cerr << "ERROR: failed to read from file "<<fileB <<std::endl
		  << "       "<< f.what() <<std::endl;
	if (rna_dataA) delete rna_dataA;
	return -1;
    }
    
    const Sequence &seqA=rna_dataA->sequence();
    const Sequence &seqB=rna_dataB->sequence();
    
    const MultipleAlignment &maA = seqA;
    const MultipleAlignment &maB = seqB;
    
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

    AnchorConstraints seq_constraints(seqA.length(),
				      seqA.sequence_anchors().single_string(),
				      seqB.length(),
				      seqB.sequence_anchors().single_string());
    
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
    arc_matches = new ArcMatches(*rna_dataA,
				 *rna_dataB,
				 min_prob,
				 (max_diff_am!=-1)?(size_type)max_diff_am:std::max(seqA.length(),seqB.length()),
				 trace_controller,
				 seq_constraints
				 );
    
    const BasePairs &bpsA = arc_matches->get_base_pairsA();
    const BasePairs &bpsB = arc_matches->get_base_pairsB();
    
    
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
    SparsificationMapper sparse_mapperA(bpsA,
    		*rna_dataA,
    		prob_unpaired_in_loop_threshold,
    		prob_basepair_in_loop_threshold,
    		false
    		);

    SparsificationMapper sparse_mapperB(bpsB,
    		*rna_dataB,
    		prob_unpaired_in_loop_threshold,
    		prob_basepair_in_loop_threshold,
    		false
    		);

    SparseTraceController sparse_trace_controller(sparse_mapperA,sparse_mapperB,trace_controller);

    //time_t stop_mapping = time (NULL);
    //cout << "time for mapping: " << stop_mapping - start_mapping << "sec " << endl;
    
    //string sequenceA= MultipleAlignment(seqA).seqentry(0).seq().to_string();
    //string sequenceB= MultipleAlignment(seqB).seqentry(0).seq().to_string();


    PatternPairMap myEPMs;

    ExactMatcher em(seqA,
		    seqB,
		    *rna_dataA,
		    *rna_dataB,
		    *arc_matches,
		    sparse_trace_controller,
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
    em.compute_arcmatch_score();
    time_t stop_computeMatrices = time (NULL);
    cout << "time for computing EPM matrices : " << stop_computeMatrices - start_computeMatrices << "sec " << endl;
    
    
    // ------------------------------------------------------------
    // Traceback
    //
    const bool DO_TRACE=true;
    //const bool DO_TRACE=false;
    if (DO_TRACE) {
	
//	if (opt_verbose) {
//	    std::cout << "Traceback."<<std::endl;
//	}
	
	gettimeofday( &tp, NULL );
	double start_trace = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
	getrusage( RUSAGE_SELF, &ruse );
	double startR_trace = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;


	if(opt_suboptimal)  cout << endl << "start suboptimal traceback..." << endl;
	else cout << endl << "start heuristic traceback..." << endl;

	//trace_controller not yet implemented for suboptimal traceback!
	if(opt_suboptimal && max_diff !=-1){
		cerr << "suboptimal traceback is not implemented yet with max_diff option " << endl;
		return 0;
	}

	em.trace_EPMs(opt_suboptimal);

	 gettimeofday( &tp, NULL );
	 double end_trace = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;
	 getrusage( RUSAGE_SELF, &ruse );
	 double endR_trace = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
	 cout << "time_wall traceback = " << setprecision(3) << end_trace - start_trace << " sec" << endl;
	 cout << "time_cpu traceback = " << setprecision(3) << endR_trace - startR_trace << " sec" << endl << endl;

	 if(epm_list_output.size()>0){
		 if (opt_verbose) { cout << "write list of traced EPMs in file..." << endl;}
		 ofstream out_EPM_file (epm_list_output.c_str());
		 out_EPM_file << myEPMs.getList() << endl;
		 out_EPM_file.close();
	 }
   }

    if(chaining){
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
    		if (opt_verbose) { cout << "write EPM chain as colored postscripts..." << endl;}
    		if (psFileA.size()==0){psFileA = seqA.seqentry(0).name()+"_EPMs.ps";}
    		if (psFileB.size()==0){psFileB = seqB.seqentry(0).name()+"_EPMs.ps";}

    		myChaining.MapToPS(maA.consensus_sequence(), maB.consensus_sequence(), myLCSEPM, psFileA,psFileB);
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

    	if(chained_epm_list_output.size()>0){
    		if (opt_verbose) { cout << "write list of chained EPMs in file..." << endl;}
    		ofstream out_chained_EPM_file (chained_epm_list_output.c_str());
    		out_chained_EPM_file << myLCSEPM.getList() << endl;
    		out_chained_EPM_file.close();
    	}
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

    stopwatch.stop("total");

    return 0;
}

#else // HAVE_LIBRNA
#include <LocARNA/aux.hh>
int
main() {
    write_errormsg_rnalib_unvailable();
    return -1;
}

#endif // HAVE_LIBRNA
