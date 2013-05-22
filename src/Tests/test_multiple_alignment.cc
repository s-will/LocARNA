#include <iostream>
#include <LocARNA/multiple_alignment.hh>
#include <LocARNA/sequence.hh>

using namespace LocARNA;

/** @file some unit tests for MultipleAlignment and Sequence
    
    Tests writing from fasta and clustalw format, basic functionality
    and simple checks
*/

int
main(int argc, char **argv) {
    
    // create simple alignment
    MultipleAlignment ma("seqA","seqB","A--CGT-U&CC-CG-CU");
    
    //! test whether ma is proper
    if (!ma.is_proper()) {
	return 1;
    }
    
    // create simple alignment from file
    MultipleAlignment *ma2=0L;
    try {
	ma2 = new MultipleAlignment("Tests/archaea.aln");
	if (!ma2->is_proper()) throw(failure("Wrong format"));
	if (ma2->empty()) throw(failure("Wrong format"));

	//! test whether ma2 is proper
	if (!ma2->is_proper()) {
	    return 2;
	}

	//! test whether ma2 has correct size
	if (!ma2->row_number() == 6) {
	    return 3;
	}
	
	//! test whether ma2 contains "fdhA" and "fruA" 
	if (!ma2->contains("fdhA") || !ma2->contains("fruA")) {
	    return 4;
	}

    } catch (failure &f) {
	return 5;
    }
    
    Sequence seq = *ma2;
    delete ma2;
    
    std::string name_str = "hdrA";
    std::string seq_str  = "GG--CACCACUCGAAGGCUA-------------AG-CCAAAGUGGUG--CU";
    seq.append(Sequence::SeqEntry(name_str,seq_str));
    
    //! test whether seq is proper
    if (!seq.is_proper()) {
    	return 6;
    }

    //! test whether seq has correct size
    if (!seq.row_number() == 7) {
    	return 7;
    }
    
    //! test whether seq contains "fdhA" and "hdrA" 
    if (!seq.contains("fdhA") || !seq.contains(name_str)) {
    	    return 8;
    }
    
    size_t index = seq.index(name_str);
    if (seq.seqentry(index).seq().to_string() != seq_str) {
	return 9;
    }
    
    bool ok=false;
    MultipleAlignment *ma3=0L;
    try {
    	ma3 = new MultipleAlignment("Tests/archaea-wrong.fa",MultipleAlignment::FASTA);
	if (!ma3->is_proper()) throw(failure("Wrong format"));
	if (ma3->empty()) throw(failure("Wrong format"));
    } catch(failure &f) {
	ok=true;
    }
    if (!ok) {
    	return 10;
    }

    ok=false;
    try {
	MultipleAlignment ma4("Tests/archaea.fa",MultipleAlignment::CLUSTAL);
	if (!ma4.is_proper()) throw(failure("Wrong format"));
	if (ma4.empty()) throw(failure("Wrong format"));
	ma4.write_debug(std::cerr);
    } catch(failure &f) {
	ok=true;
    }
    if (!ok) {
	return 11;
    }

    MultipleAlignment *ma5;
    try {
	ma5 = new MultipleAlignment("Tests/archaea.fa",MultipleAlignment::FASTA);
	if (!ma5->is_proper()) throw(failure("Wrong format"));
	if (ma5->empty()) throw(failure("Wrong format"));
    } catch(failure &f) {
	return 12;
    }

    seq=*ma5;
    delete ma5;
    
    //! test whether seq is proper
    if (!seq.is_proper()) {
    	return 106;
    }

    //! test whether seq has correct size
    if (!seq.row_number() == 7) {
    	return 107;
    }

    return 0;
}
