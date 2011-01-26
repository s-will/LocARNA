// Determine the deviation of an alignment to a reference alignment
// Copyright Sebastian Will


#include <stdlib.h>
#include <iostream>
#include <string>
#include "multiple_alignment.hh"

const std::string 
VERSION_STRING = (std::string)PACKAGE_STRING; 


void usage() {
    std::cout 
	<< "locarna_deviation - compute the deviation of an alignment to a reference alignemnt"<<std::endl<<std::endl
	<< "Usage: deviation <aln-file> <ref-aln-file>"<<std::endl<<std::endl
	<< "Options:"<<std::endl<<std::endl
	<< " <aln-file>        alignment file in clustalw format" << std::endl<< std::endl
	<< " <aln-ref-file>    reference alignment file in clustalw format" << std::endl<< std::endl
	<< std::endl
	<< "The sequences of <aln-file> have to be contained"<<std::endl
	<< "in the alignment of <aln-ref-file>." << std::endl;
    
}

int
main(int argc, char **argv) {
    
    if (argc==2) {
	if ((std::string)argv[1]=="--version" || (std::string)argv[1]=="-v") {
	    std::cout << "locarna_deviation ("<<VERSION_STRING<<")"<<std::endl;
	}
	if ((std::string)argv[1]=="--help" || (std::string)argv[1]=="-h") {
	    usage();
	}
	exit(0);
    }
    
    if (argc!=3) {
	usage();
	exit(-1);
    }
    
    MultipleAlignment ma((std::string)argv[1]);
    MultipleAlignment refma((std::string)argv[2]);
    
    std::cout << refma.deviation(ma) << std::endl;
    
    exit(0);
}
