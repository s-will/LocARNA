#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "multiple_alignment.hh"

#include "sequence.hh"
#include "alignment.hh"

#include <limits>

#include <stdlib.h>

using namespace std;

namespace LocARNA {

MultipleAlignment::MultipleAlignment(std::istream &in) {
    if (!in.good()) {
	throw(std::ifstream::failure("Cannot read input stream."));
    }

    read_aln_clustalw(in);

    create_name2idx_map();
}



MultipleAlignment::MultipleAlignment(const std::string &filename) {
    try {
	ifstream in(filename.c_str());
    
	if (!in.good()) {
	    throw(std::ifstream::failure("Cannot read file"+filename+"."));
	}
	
	read_aln_clustalw(in);
	in.close();
    } catch (std::ifstream::failure e) {
	throw failure("Cannot construct multiple alignment: "+(std::string)e.what());
    }

    create_name2idx_map();
}



MultipleAlignment::MultipleAlignment(const Sequence &sequence) {
    size_type rows = sequence.get_rows();
    
    for (size_type i=0; i<rows; i++) {
	std::string name = sequence.names()[i];
	
	std::string seq="";
	
	for (size_type j=1; j<=sequence.length(); ++j) {
	    seq += sequence[j][i];
	}
	
	alig.push_back(SeqEntry(name,seq));
    }
    
    create_name2idx_map();
}
    

MultipleAlignment::MultipleAlignment(const std::string &nameA, const std::string &nameB, const std::string &alistrings) {
    size_t separator_pos=alistrings.find("&");
    
    if (separator_pos==string::npos) {
	throw failure("No separator in string of alignment strings."); 
    }
    
    const string aliA = alistrings.substr(0,separator_pos);
    const string aliB = alistrings.substr(separator_pos+1);
    
    if (aliA.length() != aliB.length()) {
	throw failure("Alignment strings of unequal length."); 
    }
    
    alig.push_back(SeqEntry(nameA,aliA));
    alig.push_back(SeqEntry(nameB,aliB));

    create_name2idx_map();
}

MultipleAlignment::MultipleAlignment(const Alignment &alignment) {
    const std::vector<int> &a=alignment.get_a();
    const std::vector<int> &b=alignment.get_b();

    const Sequence &seqA=alignment.get_seqA();
    const Sequence &seqB=alignment.get_seqB();
    
    std::vector<std::string> aliA(seqA.get_rows(),"");
    std::vector<std::string> aliB(seqB.get_rows(),"");
    
    std::vector<int>::size_type alisize = a.size();
    
    for (size_type i=0; i<alisize; i++) {
	if ( a[i]==-1 ) {
	    for (size_type k=0; k<seqA.get_rows(); k++) {
		aliA[k] += '-';
	    }
	} else {
	    for (size_type k=0; k<seqA.get_rows(); k++) {
		aliA[k] += seqA[a[i]][k];
	    }
	}
	if ( b[i]==-1 ) {
	    for (size_type k=0; k<seqB.get_rows(); k++) {
		aliB[k] += '-';
	    }
	} else {
	    for (size_type k=0; k<seqB.get_rows(); k++) {
		aliB[k] += seqB[b[i]][k];
	    }
	}
    }
    
    const std::vector<std::string> namesA = seqA.names();
    const std::vector<std::string> namesB = seqB.names();
    
    for (size_type k=0; k<seqA.get_rows();k++) {
	alig.push_back(SeqEntry(namesA[k],aliA[k]));
    }
    for (size_type k=0; k<seqB.get_rows();k++) {
	alig.push_back(SeqEntry(namesB[k],aliB[k]));
    }
}


void
MultipleAlignment::create_name2idx_map() {
    for (std::vector<SeqEntry>::size_type i=0;
	 i<alig.size();
	 ++i) {
	name2idx[alig[i].name()]=i;
    }
}

 
void
MultipleAlignment::read_aln_clustalw(std::istream &in) {
        
    std::string name;
    std::string seqstr;
    
    std::map<std::string,std::string> seq_map;
    std::vector<std::string> names;
    
    std::string line;

    getline(in,line);
    
    // accept and ignore a header line starting with CLUSTAL
    if (line.substr(0,7) == "CLUSTAL")  {
	getline(in,line);
    }
    
    do {	
	if (line.length()>0 && line[0]!=' ') {
	    std::istringstream in(line);
	    in >> name >> seqstr;
	    
	    if (seq_map.find(name)==seq_map.end()) {
		seq_map[name] = "";
		names.push_back(name);
	    }
	    seq_map[name] += seqstr;
	}
    } while (getline(in,line) && line.substr(0,6) != "CLUSTAL"); // stop when reading a "CLUSTAL" header line again, this allows reading multiple clustal entries from one file
    
    // store the name/sequence pairs in the vector alig
    for (std::vector<std::string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
	alig.push_back( SeqEntry(*it,seq_map[*it]) );
    }
}

bool
MultipleAlignment::SeqEntry::is_gap_symbol(char c) {
    return gap_symbols.find(c)!=string::npos;
}

MultipleAlignment::size_type
MultipleAlignment::SeqEntry::length_wogaps() const {
    size_type len = 0;
    for ( size_type col = 1 ; col <= seq_.length() ; col++) {
	if ( ! is_gap_symbol( seq_[col] ) ) { 
	    len++;
	}
    }
    return len;
}


MultipleAlignment::size_type
MultipleAlignment::SeqEntry::pos_to_col(size_type pos) const {

    if (pos==0) return 0; // special case: position 0 maps to column 0

    // iterate over the columns col in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0; // the current position in the sequence (without gaps)
    
    for ( size_t col = 1 ; col <= seq_.length() ; col++) {
	    if ( ! is_gap_symbol( seq_[col] ) ) { 
		curr_pos++;
		if (pos==curr_pos) return col;
	    }
    }
    return seq_.length()+1;
}

MultipleAlignment::SeqEntry::pos_pair_t
MultipleAlignment::SeqEntry::col_to_pos(size_type col) const {
    // std::cout << "col_to_pos : " << col << "/" << seq().length() << std::endl;
    assert(col<=seq().length()+1);
    // iterate over the positions in the sequence until you've read pos number of non-gap characters
    size_t curr_pos = 0;
 
    if (col==0) { // special case
	return pos_pair_t(0,0);
    }

    if (col==seq().length()+1) { // special case
	for (size_t i = 1 ; i <= seq().length() ; i++) {
	    if (! is_gap_symbol( seq_[i] )) { 
		curr_pos++;
	    }
	}
	return pos_pair_t(curr_pos+1,curr_pos+1);
    }

    
    for (size_t i = 1 ; i <= col ; i++) {
	if (! is_gap_symbol( seq_[i] )) { 
	    curr_pos++;
	}
    }
    
    // if column col contains a gap, then return (curr_pos, curr_pos + 1)
    // if column col contains a non-gap, then return (curr_pos, curr_pos)
    if (is_gap_symbol( seq_[col] ) ) {
    	return pos_pair_t(curr_pos, curr_pos + 1);
    } else {
    	return pos_pair_t(curr_pos, curr_pos);
    }
}

bool
MultipleAlignment::is_proper() const {
	// Iterate through all sequences (with gaps) and check for equal lengths
	bool proper = true;
	std::vector<SeqEntry>::const_iterator it = alig.begin();
	size_t length = it->seq().length();
	++it;
	while (proper && (it != alig.end())) {
	    proper = proper && (length == it->seq().length());
	}
	return proper;
}

bool
MultipleAlignment::contains(string name) const {
    for (std::vector<SeqEntry>::const_iterator it = begin(); it != end(); ++it) {
	if (name == it->name()) { return true; }
    }
    return false;
}

		   
MultipleAlignment::size_type
MultipleAlignment::deviation2(const string1 &a1,
			      const string1 &a2,
			      const string1 &ref1,
			      const string1 &ref2
			      ) {
    
    //determine the cut deviation between to pariwise alignments
    // compute the maximum over the minimum distance of cut c1 in a to a cut c2 in b 
    
    size_type d=0;
    
    size_type i1=0;
    size_type j1=0;
    for (size_type c1=0; c1<=a1.length(); c1++) {
	
	size_type dprime=numeric_limits<size_type>::max();
	
	if (c1>0) {
	    if (!SeqEntry::is_gap_symbol(a1[c1])) i1++;
	    if (!SeqEntry::is_gap_symbol(a2[c1])) j1++;
	}

	size_type i2=0;
	size_type j2=0;
	for (size_type c2=0; c2<=ref1.length(); c2++) {
	    
	    if (c2>0) {
		if (!SeqEntry::is_gap_symbol(ref1[c2])) i2++;
		if (!SeqEntry::is_gap_symbol(ref2[c2])) j2++;
	    }

	    size_type dc =
		abs((long int)i1-(long int)i2)
		+abs((long int)j1-(long int)j2);

	    dprime = std::min(dprime,dc);
	}

	d=std::max(d,dprime);
    }
    return d;
}


MultipleAlignment::size_type
MultipleAlignment::deviation(const MultipleAlignment &ma) const {
    
    //determine maximal deviation2 of pairwise subalignments

    size_type d=0;
    
    //traverse all pairs of sequences in *this and ma
    for (size_type x=0; x<ma.alig.size(); x++) {
	const std::string &namex = ma.seqentry(x).name();
	for (size_type y=x+1; y<ma.alig.size(); y++) {
	    const std::string &namey = ma.seqentry(y).name();
	    
	    size_type d2 =
		deviation2(ma.seqentry(x).seq(),ma.seqentry(y).seq(),
			   seqentry(namex).seq(),seqentry(namey).seq()
			   );

	    d = std::max(d, d2);
	}
    }
    return d;
}

		   
void
MultipleAlignment::print_debug(ostream &out) const {
    
    for (size_type i=0; i<alig.size(); ++i) {
	out << alig[i].name() << " \t" << alig[i].seq().to_string()<<std::endl;
    } 
    
}


} // end namespace 
