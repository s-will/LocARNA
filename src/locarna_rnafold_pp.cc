/************************************************************/
/**
 * \file locarna_rnafold_pp.cc
 * \brief Compute and write base pair probabilities in pp-format.
 *
 * Special RNAfold-like program that computes a probability dot plot
 * for a input RNA sequence and writes a file in the
 * LocARNA-internal pp-format.
 * 
 * This program is part of the LocARNA package. It is intended for
 * computing pair probabilities of the input sequences.
 *
 * Reads sequence in fasta from cin and writes pp-files to cout
 *
 * command line argument --TEST provides a way to test for linking to
 * the ViennaLib. (This should be eventually replaced by a less
 * idiosyncratic mechanism.)
 *
 * @todo use RnaEnsemble to avoid code duplication (caution: at the
 * time of writing, RnaEnsemble does not support constraints yet.)
 */
/************************************************************/

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <stdlib.h>
#include <iostream>
#include <fstream>

#ifdef HAVE_LIBRNA

#include <math.h>

#include <string.h>
#include <sstream>
#include <string>

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils.h>
}

#include <LocARNA/options.hh>
#include <LocARNA/multiple_alignment.hh>
#include <LocARNA/pfold_params.hh>
#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/rna_data.hh>
#include <LocARNA/ext_rna_data.hh>



using namespace LocARNA;

    /*
	 * Concerning rnafold_pp: we need to support at least the following
    options/arguments

    input file: Unlike RNAfold, we want to
    interpret multiple sequences in fasta format as one alignment. We need
    to check that the input alignment is proper.

    -C           use structural constraints. [Actually, this becomes
    obsolete; the tool could use constraints, if they are present in the
    input. However, the tool could print a warning to std::cerr if
    constraints are found but -C is missing. We could also add an option
    like --ignore-constraints. Alternatively, the behavior could be to use
    structure constraints only when -C is given and otherwise ignore them.
    The latter would be closer to the behavior of RNAfold; maybe a good idea:)].



     input file: input format should be either FASTA or CLUSTAL;
    if no filename is provided or filename=="-" then the input is read from stdin.

    -C           use structural constraints. if it is intended to use structural constraints(annotations), this parameter should be set, otherwise the given constraints will be ignored.

    --noLP       forbid lonely base pairs

    --stacking   compute stacking terms

    --in-loop  compute in-loop probabilities

    -p,--min-prob    set cutoff probability (default 0.0005)

    --prob_unpaired_in_loop_threshold   Threshold for prob_unpaired_in_loop

    --prob_basepair_in_loop_threshold    Threshold for prob_basepair_in_loop


    -o,--output=<file>                       PP output [default: output to
    std::cout; if file is give output only to file]

    --force-alifold    use alifold even for single sequences. [The standard
    behavior could be to use alifold unless the input contains only a single
    sequence.]

    -w,--width=<columns>(120)         Output width


    We need mechanisms to catch wrong input; some of these things are
    already handled in the lib classes, which raise exceptions on errors.
    You should catch such exceptions and print errors on wrong input.
    Currently the constraint string is not  checked for allowed characters
    (see RNAfold man page) and balanced parentheses; probably this could
    crash the fold routine of RNAlib.  Can you please suggest a solution? I
    think, we should check this in the set_annotation method of
    MultipleAlignment, which would make it impossible to put
    MultipleAlignment objects into some inconsistent state.
	 * */

//! \brief Structure for command line parameters of locarna
//!
//! Encapsulating all command line parameters in a common structure
//! avoids name conflicts and makes downstream code more informative.
//!
struct command_line_parameters {
    bool opt_help; //!< whether to print help
    bool opt_version; //!< whether to print version
    bool opt_verbose; //!< whether to print verbose output
    std::string input_file; //!< input_file
    bool use_struct_constraints; //!< -C use structural constraints
    bool no_lonely_pairs; //!< no lonely pairs option
    bool opt_stacking; //!< whether to stacking
    bool opt_in_loop; //!< whether to compute in-loop probabilities
    double min_prob;  //! < only pairs with a probability of at least min_prob are taken into account
    double prob_unpaired_in_loop_threshold; //!< threshold for prob_unpaired_in_loop
    double prob_basepair_in_loop_threshold; //!< threshold for prob_basepait_in_loop
    std::string output_file; //!< output file name
    bool force_alifold; //!< use alifold even for single sequences.
    uint width; //!< Output width
    bool test; //!< Test the binary and return 1 if it the functionality is available
};
//! \brief holds command line parameters of locarna
command_line_parameters clp;
//longname,shortname,flag,arg_type,argument,default,argname,description
//! defines command line parameters
option_def my_options[] = {
    {"help",'h',&clp.opt_help,O_NO_ARG,0,O_NODEFAULT,"","Help"},
    {"version",'V',&clp.opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&clp.opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},
    {"input",0,0,O_ARG_STRING,&clp.input_file,"-","f","Input file"},
    {"use-struct-constraints",'C',&clp.use_struct_constraints, O_NO_ARG, 0, O_NODEFAULT, "","Use structural constraints"},
    {"noLP",0,&clp.no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs"},
    {"stacking",0,&clp.opt_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use stacking"},
    {"in-loop",0,&clp.opt_in_loop,O_NO_ARG,0,O_NODEFAULT,"","Compute in-loop probabilities"},
    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,"0.0005","prob","Minimal probability"},
    {"prob_unpaired_in_loop_threshold",0,0,O_ARG_DOUBLE,&clp.prob_unpaired_in_loop_threshold,"0.0005","threshold","Threshold for prob_unpaired_in_loop"},
    {"prob_basepair_in_loop_threshold",0,0,O_ARG_DOUBLE,&clp.prob_basepair_in_loop_threshold,"0.0005","threshold","Threshold for prob_basepair_in_loop"}, //todo: is the default threshold value reasonable?
    {"output",'o',0,O_ARG_STRING,&clp.output_file,"","f","Output file"},
    {"force-alifold",0,&clp.force_alifold,O_NO_ARG,0,O_NODEFAULT,"","Force alifold for single sequnces"},
    {"width",'w',0,O_ARG_INT,&clp.width,"120","size","Output width"},
    {"test",0,&clp.test,O_NO_ARG,0,O_NODEFAULT,"","Test avialability"},
    {"",0,0,0,0,O_NODEFAULT,"",""}

};





/** 
 * \brief Main function of locarna_rnafold_pp when Vienna RNA lib is linked
 */
int
main(int argc, char **argv) {
    
    // ------------------------------------------------------------
    // Process options

    bool process_success=process_options(argc,argv,my_options);

    if (clp.opt_help) {
	std::cout << "locarna_rnafold_pp -- compute RNA pair probabilities and write in pp-format" << std::endl;
//	std::cout << VERSION_STRING<<std::endl<<std::endl;
	print_help(argv[0],my_options);
	return 0;
    }

    if (clp.opt_version || clp.opt_verbose) {
	std::cout << "locarna_rnafold_pp "/*<< VERSION_STRING<<")"*/<<std::endl;
	if (clp.opt_version) return 0; else std::cout <<std::endl;
    }
    //test whether vrna is available
    if ( clp.test ) {
	std::cout << "1";
	return 0;
    }

    if (!process_success) {
	std::cerr << "ERROR --- "
		<<O_error_msg<<std::endl;
	std::cout << "USAGE: " << std::endl;
	print_usage(argv[0],my_options);
	std::cout << std::endl;
	return -1;
    }
    //todo: 	check that input is proper and catch the wrong inputs
    MultipleAlignment::FormatType::type input_format;
    MultipleAlignment* mseq = NULL;
    // try fasta format
    bool failed = false;
    try {

	if (clp.input_file.compare("-") != 0)
	{
	    mseq = new MultipleAlignment(clp.input_file, MultipleAlignment::FormatType::FASTA);
	}
	else
	{
	    mseq = new MultipleAlignment(std::cin, MultipleAlignment::FormatType::FASTA);
	}
	// even if reading does not fail, we still want to
	// make sure that the result is reasonable. Otherwise,
	// we assume that the file is in a different format.
	if (! mseq->is_proper() || mseq->empty() ) {
	    failed=true;
	}
    } catch (failure &f) {
	failed=true;
//	std::cerr << "Not a FASTA file" << std::endl;
//	std::cerr << f.what() << std::endl;
    }

    if (!failed)
    {
	input_format = MultipleAlignment::FormatType::FASTA;
    }
    else
    { //
	failed=false;
	try{

	    try
	    {
		// std::cerr << "Try reading clustal "<<filename<<" ..."<<std::endl;
		if (clp.input_file.compare("-") != 0)
		{
		    mseq = new MultipleAlignment (clp.input_file, MultipleAlignment::FormatType::CLUSTAL);
		}
		else
		{
		    mseq = new MultipleAlignment (std::cin, MultipleAlignment::FormatType::CLUSTAL);

		// make sure that the result is reasonable. Otherwise,
		}
		// even if reading does not fail, we still want to
		// we assume that the file is in a different format.
		if (!mseq->is_proper() || mseq->empty() ) {
		    failed=true;
		}
	    } catch (syntax_error_failure &f) {
		throw failure((std::string)"Cannot read input data from clustal file.\n\t"+f.what());
	    }
	}
	catch (failure &f) {
	    failed=true;
	    std::cerr << f.what() << std::endl;
	}
	if (failed)
	{
	    std::cerr << "Error in input format" << std::endl;
	    return -1;
	}
    }

    bool use_alifold = true;

    // if the input has only one sequence and forced to do alifold
    if (!clp.force_alifold && mseq->num_of_rows()==1)
	use_alifold = false;
    if(mseq->has_annotation(MultipleAlignment::AnnoType::structure) && !clp.use_struct_constraints)
    {
	std::cerr << "Warning locarna_rnafold_pp: structure constraints will be ignored" << std::endl;
	mseq->set_annotation(MultipleAlignment::AnnoType::structure,SequenceAnnotation() );

    }
    
    PFoldParams pfoldparams(clp.no_lonely_pairs, clp.opt_stacking);

    RnaEnsemble rna_ensemble(*mseq,pfoldparams, clp.opt_in_loop, use_alifold);

    size_t len = mseq->length();

    int num_bps = 0;
    // ----------------------------------------
    // print base pair probabilities
    std::cout << "//**** print base pair probabilities****" << std::endl;
    for( size_t i=1; i <= len; i++ ) {
	for( size_t j=i+3+1; j <= len; j++ ) {

	    double p = rna_ensemble.arc_prob(i,j);
	    if (p > clp.min_prob) { // apply filter
		std::cout <<  i << "\t" << j<<"\t" << p << std::endl;
		num_bps++;
	    }
	}
    }
    std::cout << "//****num_bps: " << num_bps << "********************************" << std::endl;

    if(mseq->num_of_rows() == 1)
    {
	const char *sequence=mseq->seqentry(0).seq().str().c_str(); //wow such a sequence call!
	char *structure = (char *)space(mseq->length()+1);



	dangles=2;

	//for long sequence, compute special scaling factor
	pf_scale=1.0;
	fold_constrained=false;
	std::cout << "locarna_rna_fold:\tseq:" << sequence << std::endl << "sturcture: " << structure << std::endl;

	//double energy =
	pf_fold(sequence,structure);

	FLT_OR_DBL *probs = export_bppm();
	plist *pl;
	assign_plist_from_pr(&pl, probs, mseq->length(), clp.min_prob);


	// read the pair_info structures from array pl
	// and write the base pair probabilities to out
	//
	for(plist *pair=pl; pair->i!=0 && pair->j!=0; pair++) {
	    if ( pair->p  >  clp.min_prob) {

		double pnew = rna_ensemble.arc_prob(pair->i,pair->j);
		std::cout << pair->i << "\t"
			<< pair->j << "\t"
			<< pair->p << "\t" << sqrt(pair->p) << "\t" << pnew << "\t" << sqrt(pnew) <<std::endl;
	    }
	}
	free(pl);

	free_arrays();
	free_pf_arrays();

	free(structure);
    }
    // write pp file

    //set the appropriate ostream from input_file or std::cout
    std::ofstream of;
    std::streambuf * buff;
    if (clp.output_file.length() == 0)
    {
	buff = std::cout.rdbuf();
    }
    else
    {
	of.open(clp.output_file.c_str());
	buff = of.rdbuf();
    }
    std::ostream out_stream(buff);

    RnaData* rna_data_ptr;
    if (clp.opt_in_loop)
    {
	ExtRnaData ext_rna_data(rna_ensemble, clp.min_prob, clp.prob_basepair_in_loop_threshold, clp.prob_unpaired_in_loop_threshold);
	ext_rna_data.write_pp(out_stream, clp.min_prob, clp.prob_basepair_in_loop_threshold, clp.prob_unpaired_in_loop_threshold);//toask: why should give the theresholds again
    }
    else
    {
	RnaData rna_data(rna_ensemble, clp.min_prob);
	rna_data.write_pp(out_stream,clp.min_prob);
    }

    return 0;
}

#else
/** 
 * \brief Main in case the programs are not linked to the Vienna RNA
 */
int
main(int argc, char **argv) {

    bool process_success=process_options(argc,argv,my_options);

       if ( clp.test ) {
	std::cout << "0";
	return -1;
    }
    std::cerr << "ERROR: locarna_rnafold_pp requires linking against Vienna librna.\n";
    std::cerr << "This program was compiled without configure option --enable-librna."<<std::endl;
    std::cerr << "Please reconfigure and recompile to activate this program. \n";
    return -1;
}
#endif // HAVE_LIBRNA
