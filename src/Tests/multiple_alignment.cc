#include <cassert>
#include <iostream>
#include <sstream>
#include <LocARNA/multiple_alignment.hh>
#include <LocARNA/alignment.hh>
#include <LocARNA/sequence.hh>

using namespace LocARNA;

#undef NDEBUG

/** @file some unit tests for MultipleAlignment and Sequence
    
    Tests writing from fasta and clustalw format, basic functionality
    and simple checks
*/

int
main(int argc, char **argv) {
    
    // create simple alignment
    MultipleAlignment ma("seqA","seqB",
			 "A-CGT-U",
			 "CCCG-CU");
    
    //! test whether ma is proper
    assert(ma.is_proper());
    
    // create simple alignment from file
    MultipleAlignment *ma2=0L;
    
    ma2 = new MultipleAlignment("Tests/archaea.aln");
    if (!ma2->is_proper()) throw(failure("Wrong format"));
    if (ma2->empty()) throw(failure("Wrong format"));
    
    //! test whether ma2 is proper
    assert(ma2->is_proper());
    
    //! test whether ma2 has correct size
    assert(ma2->num_of_rows() == 6);
    
    //! test whether ma2 contains "fdhA" and "fruA" 
    assert (ma2->contains("fdhA") && ma2->contains("fruA"));
    
    Sequence seq = (*ma2).as_sequence();
    delete ma2;
    
    std::string name_str = "hdrA";
    std::string seq_str  = "GG--CACCACUCGAAGGCUA-------------AG-CCAAAGUGGUG--CU";
    seq.append(Sequence::SeqEntry(name_str,seq_str));
    
    //! test whether seq is proper
    assert(seq.is_proper());

    //! test whether seq has correct size
    assert(7 == seq.num_of_rows());
    
    //! test whether seq contains "fdhA" and "hdrA" 
    assert(seq.contains("fdhA") && seq.contains(name_str));
    
    size_t index = seq.index(name_str);
    assert(seq.seqentry(index).seq().str() == seq_str);
    
    bool ok=false;
    MultipleAlignment *ma3=0L;

    try {
    	ma3 = new MultipleAlignment("Tests/archaea-wrong.fa",MultipleAlignment::FormatType::FASTA);
	if (!ma3->is_proper()) throw(failure("Wrong format"));
	if (ma3->empty()) throw(failure("Wrong format"));
    } catch(failure &f) {
	ok=true;
    }
    assert(ok);

    ok=false;
    try {
	MultipleAlignment ma4("Tests/archaea.fa",MultipleAlignment::FormatType::CLUSTAL);
	if (!ma4.is_proper()) throw(failure("Wrong format"));
	if (ma4.empty()) throw(failure("Wrong format"));
	//ma4.write_debug(std::cerr);
    } catch(failure &f) {
	ok=true;
    }
    assert(ok);

    MultipleAlignment *ma5;
    
    ma5 = new MultipleAlignment("Tests/archaea.fa",MultipleAlignment::FormatType::FASTA);
    if (!ma5->is_proper()) throw(failure("Wrong format"));
    if (ma5->empty()) throw(failure("Wrong format"));
    
    seq=(*ma5).as_sequence();
    delete ma5;
    
    //! test whether seq is proper
    assert(seq.is_proper());

    //! test whether seq has correct size
    assert(seq.num_of_rows() == 3);
    
    {
	// write and read again test
	
	std::ostringstream out; 
	seq.write(out);

	std::istringstream in(out.str()); 
	MultipleAlignment seq2(in,MultipleAlignment::FormatType::CLUSTAL);
	
	assert(seq.num_of_rows() == seq2.num_of_rows());
	assert(seq.length() == seq2.length());
    }

    return 0;
}
