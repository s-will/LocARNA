#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "alphabet.hh"
#include "sequence.hh"
#include "rna_data.hh"
#include "basepairs.hh"
#include "alignment.hh"
#include "multiple_alignment.hh"
#include "sequence_annotation.hh"

#include <limits>

#include <math.h>
#include <stdlib.h>

using namespace std;

namespace LocARNA {
    
    // allowed tags are part of the extended aln and pp file format definitions
    std::vector<std::string>
    initial_annotation_tags() {
	std::vector<std::string> as;
	as.resize(3);
	
	as[MultipleAlignment::AnnoType::structure]="S";
	as[MultipleAlignment::AnnoType::fixed_structure]="FS";
	as[MultipleAlignment::AnnoType::anchors]="A";
	
	return as;
    }
    const std::vector<std::string> MultipleAlignment::annotation_tags = initial_annotation_tags();
    
    MultipleAlignment::MultipleAlignment() 
	: alig_(),
	  annotations_(),
	  name2idx_() {
    }
    
    MultipleAlignment::MultipleAlignment(std::istream &in, FormatType::type format)
    	: alig_(),
	  annotations_(),
	  name2idx_() {
	
	if (!in.good()) {
	    throw(failure("Cannot read input stream."));
	}

	if (format==FormatType::FASTA) {
	    read_aln_fasta(in);
	} else if (format==FormatType::CLUSTAL) {
	    read_aln_clustalw(in);
	} else {
	    throw failure("Unknown format.");
	}
	
	create_name2idx_map();
    }
    
    MultipleAlignment::MultipleAlignment(const std::string &filename, FormatType::type format)
	: alig_(),
	  annotations_(),
	  name2idx_() {
	
	try {
	    ifstream in(filename.c_str());
    
	    if (!in.is_open()) {
		throw(std::ifstream::failure("Cannot open file "+filename+" for reading."));
	    }
	
	    if (format==FormatType::FASTA) {
		read_aln_fasta(in);
	    } else if (format == FormatType::CLUSTAL) {
		read_aln_clustalw(in);
	    } else {
		throw failure("Unknown format.");
	    }

	
	    in.close();
	} catch (std::ifstream::failure &e) {
	    throw failure("Cannot construct multiple alignment: "+(std::string)e.what());
	}
		
	create_name2idx_map();
    }

    MultipleAlignment::MultipleAlignment(const std::string &name, 
					 const std::string &sequence)
	: alig_(),
	  annotations_(),
	  name2idx_() {
	
	alig_.push_back(SeqEntry(name,sequence));
    	
    	create_name2idx_map();
    }

    MultipleAlignment::MultipleAlignment(const std::string &nameA,
					 const std::string &nameB, 
					 const std::string &aliA,
					 const std::string &aliB)
	: alig_(),
	  annotations_(),
	  name2idx_() {
	
	if (aliA.length() != aliB.length()) {
    	    throw failure("Alignment strings of unequal length."); 
    	}
	
    	alig_.push_back(SeqEntry(nameA,aliA));
    	alig_.push_back(SeqEntry(nameB,aliB));
	
    	create_name2idx_map();
    }

    MultipleAlignment::MultipleAlignment(const Alignment &alignment, bool only_local)
	:alig_(),
	 annotations_(),
	 name2idx_() {
	
	SequenceAnnotation
	    anchors(alignment.alignment_edges(only_local),
		 alignment.seqA().annotation(AnnoType::anchors),
		 alignment.seqB().annotation(AnnoType::anchors));
	
	if (! anchors.duplicate_names() ) {
	    set_annotation(AnnoType::anchors, anchors);
	}
			
	init(alignment.alignment_edges(only_local),
	     alignment.seqA(),
	     alignment.seqB());
    }
    
    MultipleAlignment::MultipleAlignment(const AlignmentEdges &edges,
					 const Sequence &seqA,
					 const Sequence &seqB)
	:alig_(),
	 annotations_(),
	 name2idx_() {

	SequenceAnnotation 
	    anchors(edges,
		 seqA.annotation(AnnoType::anchors),
		 seqB.annotation(AnnoType::anchors));
	
	if (! anchors.duplicate_names() ) {
	    set_annotation(AnnoType::anchors,
			   anchors);
	}
	
	init(edges,seqA,seqB);
    }
    
    void
    MultipleAlignment::init(const AlignmentEdges &edges,
			    const Sequence &seqA,
			    const Sequence &seqB) {
	
	std::vector<std::string> aliA(seqA.num_of_rows(),"");
	std::vector<std::string> aliB(seqB.num_of_rows(),"");
    
	std::vector<int>::size_type alisize = edges.size();
    
	for (size_type i=0; i<alisize; i++) {
	    if ( edges.first[i].is_gap() ) {
		for (size_type k=0; k<seqA.num_of_rows(); k++) {
		    aliA[k] += gap_symbol(edges.first[i].gap());
		}
	    } else {
		for (size_type k=0; k<seqA.num_of_rows(); k++) {
		    aliA[k] += seqA.column(edges.first[i])[k];
		}
	    }
	    if ( edges.second[i].is_gap() ) {
		for (size_type k=0; k<seqB.num_of_rows(); k++) {
		    aliB[k] += gap_symbol(edges.second[i].gap());
		}
	    } else {
		for (size_type k=0; k<seqB.num_of_rows(); k++) {
		    aliB[k] += seqB.column(edges.second[i])[k];
		}
	    }
	}
    
	for (size_type k=0; k<seqA.num_of_rows();k++) {
	    alig_.push_back(SeqEntry(seqA.seqentry(k).name(),aliA[k]));
	}
	for (size_type k=0; k<seqB.num_of_rows();k++) {
	    alig_.push_back(SeqEntry(seqB.seqentry(k).name(),aliB[k]));
	}

    	create_name2idx_map();
    }

    
    MultipleAlignment::~MultipleAlignment() {
    }

    const Sequence & MultipleAlignment::as_sequence() const {
	return static_cast<const Sequence &>(*this);
    }

    void
    MultipleAlignment::create_name2idx_map() {
	for (std::vector<SeqEntry>::size_type i=0;
	     i<alig_.size();
	     ++i) {
	    name2idx_[alig_[i].name()]=i;
	}
    }

    void
    MultipleAlignment::read_aln_clustalw(std::istream &in) {

	std::string name;
	std::string seqstr;
    
	std::map<std::string,std::string> seq_map;
	std::vector<std::string> names;
    
	std::string line;

	std::vector<std::string> anchors;

	std::string structure_string="";
	std::string fixed_structure_string="";
	
	alig_.clear();
	
	get_nonempty_line(in,line);
	
	// accept and ignore a header line starting with CLUSTAL
	if (has_prefix(line,"CLUSTAL"))  {
	    get_nonempty_line(in,line);
	}
	
	const std::string anchors_tag   = "#"+annotation_tags[AnnoType::anchors];
	const std::string structure_tag = "#"+annotation_tags[AnnoType::structure];
	const std::string fixed_structure_tag = "#"+annotation_tags[AnnoType::fixed_structure];
	
	do {
	    if (line[0]=='#') {
		if (has_prefix(line,"#END")) { // recognize END for use in pp files
		    // section end
		    break;
		}
		else if (has_prefix(line,anchors_tag)) {
		    // anchor constraint
		    std::istringstream in(line.substr(2));
		    int idx;
		    std::string astr;
		    in >> idx >> astr;
		    
		    if (idx < 1) {
			throw syntax_error_failure("Invalid index in anchor specification.");
		    }
		    
		    if ((unsigned int)idx > anchors.size() + 1) {
			throw syntax_error_failure("Non-contiguous anchor specification.");
		    }
		    		    
		    if ((unsigned int)idx > anchors.size()) {
			anchors.resize(idx);
		    }
		    
		    anchors[idx-1] += astr;
		} else if (has_prefix(line,structure_tag) 
			   || has_prefix(line,fixed_structure_tag)) {
		    std::istringstream in(line);
		    std::string tag;
		    std::string str;
		    in >> tag >> str;
		    
		    if ( in.fail() ) {
			throw syntax_error_failure("Invalid tag with structure prefix "
						   +structure_tag+".");
		    }
		    
		    if (tag==structure_tag) {
			structure_string += str;
		    } else if (tag==fixed_structure_tag) {
		    	fixed_structure_string += str;
		    } else {
			throw syntax_error_failure("Unknown tag with structure prefix "
						   +structure_tag+".");
		    }
		}
	    } else {
	    	std::istringstream in(line);
		in >> name >> seqstr;
		if ( in.fail() || name.empty() || seqstr.empty()) {
		    throw syntax_error_failure("Unexpected line while reading clustal format");
		}
		
		if (seq_map.find(name)==seq_map.end()) {
		    seq_map[name] = "";
		    names.push_back(name);
		}
		seq_map[name] += seqstr;
	    }
	    // stop when reading a "CLUSTAL" header line again, this
	    // allows reading multiple clustal entries from one file
	} while (get_nonempty_line(in,line) 
		 && !has_prefix(line,"CLUSTAL"));
    
	// store the name/sequence pairs in the vector alig
	for (std::vector<std::string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
	    alig_.push_back( SeqEntry(*it,seq_map[*it]) );
	}

	// check anchor constraints
	if (anchors.size()>0) {
	    size_t len = anchors[0].size();
	    if (len != length()) {
		throw syntax_error_failure("Anchor string(s) of wrong length.");
	    }
	    for (size_t i=1; i<anchors.size(); i++) {
		if ( anchors[i].size() != len ) {
		    throw syntax_error_failure("Anchor strings of unequal length.");
		}
	    }
	    
	    // initialize sequence_anchors_ with anchors
	    set_annotation( AnnoType::anchors, SequenceAnnotation(anchors) );
	}
	
	// check and set structure string
	if (!structure_string.empty()) {
	    if (structure_string.length() != length()) {
		throw syntax_error_failure("Structure annotation of wrong length.");
	    }
	    try {
		set_annotation( AnnoType::structure,
				SequenceAnnotation(structure_string) );
	    } catch(failure &f) {
		throw syntax_error_failure((std::string)"Wrong structure annotation."
					   +f.what());
	    }
	}

	// check and set fixed structure string
	if (!fixed_structure_string.empty()) {
	    if (fixed_structure_string.length() != length()) {
		throw syntax_error_failure("Fixed structure annotation of wrong length.");
	    }
	    try {
		set_annotation( AnnoType::fixed_structure,
				SequenceAnnotation(fixed_structure_string) );
	    } catch(failure &f) {
		throw syntax_error_failure((std::string)"Wrong fixed structure annotation."
					   +f.what());
	    }
	}
	
    }
    
    void
    MultipleAlignment::read_aln_fasta(std::istream &in) {
	std::string name;
	std::string description;
	
	std::string line;
    
	alig_.clear();
	
	getline(in,line);
	
	while(in) {
	    
	    if (line.length()>0 && line[0]=='>') {
		std::istringstream sin(line);
		sin.get(); // this eats '>'
		sin >> name; // gets the first non-whitespace substring after '>' of the line
		
		if (sin.fail() || name.empty()) {
		    throw syntax_error_failure("Cannot read sequence header after '>'");
		}
		
		// todo: this does not eat off blanks at begining of description yet
		std::stringbuf sb;
		sin.get(sb);
		description=sb.str();

		std::string seqstr="";
    		getline(in,line);
		while((in) && (line.size()==0 || line[0]!='>')) {
		    // remove whitespace and add to seqstr
		    std::istringstream sin(line);
		    std::string seqstr1;
		    while(sin >> seqstr1) {
			seqstr += seqstr1;
		    }
		    getline(in,line);
		}
		
		alig_.push_back( SeqEntry(name,description,seqstr) );
	    } else {
		throw syntax_error_failure("Unexpected line while reading fasta");
	    }
	}
    }
    

    size_type
    MultipleAlignment::SeqEntry::length_wogaps() const {
	size_type len = 0;
	for ( size_type col = 1 ; col <= seq_.length() ; col++) {
	    if ( ! is_gap_symbol( seq_[col] ) ) { 
		len++;
	    }
	}
	return len;
    }
    
    
    size_type
    MultipleAlignment::SeqEntry::pos_to_col(pos_type pos) const {

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
    MultipleAlignment::SeqEntry::col_to_pos(pos_type col) const {
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
    
    void MultipleAlignment::normalize_rna_symbols() {
	for (std::vector<SeqEntry>::iterator it = alig_.begin();
	     alig_.end() != it; ++it) {
	    std::string seq = it->seq().str();
	    normalize_rna_sequence(seq);
	    it->set_seq( seq );
	}
    }

    const SequenceAnnotation &
    MultipleAlignment::annotation(const AnnoType::type &annotype) const {

	assert(0<=annotype && annotype<num_of_annotypes());

    	annotation_map_t::const_iterator it = annotations_.find(annotype);

    	if ( it != annotations_.end() ) {
    	    return it->second;
    	} else {
	    return SequenceAnnotation::empty_instance();
    	}
    }
    bool
    MultipleAlignment::is_proper() const {
	if (empty()) return true; // empty alignment is proper
	
	// if not empty check that all alignment strings have the same length	
	size_t length = alig_.front().seq().length();
	
	for (std::vector<SeqEntry>::const_iterator it = alig_.begin();
	     alig_.end() != it; ++it) {
	    if (it->seq().length() != length) return false;
	}
	return true;
    }

    bool
    MultipleAlignment::contains(string name) const {
	for (std::vector<SeqEntry>::const_iterator it = begin(); it != end(); ++it) {
	    if (name == it->name()) { return true; }
	}
	return false;
    }

    size_type
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
		if (!is_gap_symbol(a1[c1])) i1++;
		if (!is_gap_symbol(a2[c1])) j1++;
	    }

	    size_type i2=0;
	    size_type j2=0;
	    for (size_type c2=0; c2<=ref1.length(); c2++) {
	    
		if (c2>0) {
		    if (!is_gap_symbol(ref1[c2])) i2++;
		    if (!is_gap_symbol(ref2[c2])) j2++;
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


    size_type
    MultipleAlignment::deviation(const MultipleAlignment &ma) const {
    
	//determine maximal deviation2 of pairwise subalignments
    
	size_type d=0;
    
	//traverse all pairs of sequences in *this and ma
	for (size_type x=0; x<ma.alig_.size(); x++) {
	    const std::string &namex = ma.seqentry(x).name();
	    for (size_type y=x+1; y<ma.alig_.size(); y++) {
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
    
    double
    MultipleAlignment::sps(const MultipleAlignment &ma, bool compalign) const {
	double sps=0.0;
	
	
	//traverse all pairs of sequences in *this and ma
	for (size_type x=0; x<ma.alig_.size(); x++) {
	    const std::string &namex = ma.seqentry(x).name();
	    for (size_type y=x+1; y<ma.alig_.size(); y++) {
		const std::string &namey = ma.seqentry(y).name();
		
		size_t len1 = seqentry(namex).length_wogaps(); 
		size_t len2 = seqentry(namey).length_wogaps(); 
		
		// compute pairwise score
				
		double match_score =
		    pairwise_match_score(ma.seqentry(x),ma.seqentry(y),
					 seqentry(namex),seqentry(namey),
					 compalign
					 );
		
		if (compalign) {
		    match_score +=
			pairwise_match_score(ma.seqentry(y),ma.seqentry(x),
					     seqentry(namey),seqentry(namex),
					     compalign
					     );
		
		    match_score /= len1+len2;
		} else {
		    match_score *= 2;
		    match_score /= 
			count_matches(ma.seqentry(x),ma.seqentry(y))
			+ count_matches(seqentry(namey),seqentry(namex));
		}
			
		sps += match_score;
	    }
	}
	
	size_t K=ma.num_of_rows();
	
	return sps*2.0/K/(K-1);
    }

    double
    MultipleAlignment::avg_deviation_score(const MultipleAlignment &ma) const {
	double score=0.0;
	
	
	//traverse all pairs of sequences in *this and ma
	for (size_type x=0; x<ma.alig_.size(); x++) {
	    const std::string &namex = ma.seqentry(x).name();
	    for (size_type y=x+1; y<ma.alig_.size(); y++) {
		const std::string &namey = ma.seqentry(y).name();
		
		size_t len1 = seqentry(namex).length_wogaps(); 
		size_t len2 = seqentry(namey).length_wogaps(); 
		
		// compute pairwise score
		double pw_score =
		    pairwise_deviation_score(ma.seqentry(x),ma.seqentry(y),
					     seqentry(namex),seqentry(namey));
		
		pw_score +=
		    pairwise_deviation_score(ma.seqentry(y),ma.seqentry(x),
					     seqentry(namey),seqentry(namex));
		
		pw_score /= len1+len2;
			
		score += pw_score;
	    }
	}
	
	size_t K=ma.num_of_rows();
	
	return score*2.0/K/(K-1);
    }

    
    std::vector<int> 
    MultipleAlignment::match_vector(const string1 &s,
				    const string1 &t) {
	size_t len=s.length();
	assert(len == t.length());
	
	size_t posT=0;
	
	std::vector<int> res;
	res.push_back(-1);

	for (size_t col=1; col<=len; col++) {
	    int entry=-1;
	    if (!is_gap_symbol(t[col])) {
		posT++;
		entry=posT;
	    }
	    if (!is_gap_symbol(s[col])) {
		res.push_back(entry);
	    }
	}

	return res;
    }

    std::vector<int> 
    MultipleAlignment::match_vector2(const string1 &s,
				    const string1 &t) {
	size_t len=s.length();
	assert(len == t.length());
	
	size_t posT=0;
	
	std::vector<int> res;
	res.push_back(-1);

	for (size_t col=1; col<=len; col++) {
	    if (!is_gap_symbol(t[col])) {
		posT++;
	    }
	    if (!is_gap_symbol(s[col])) {
		res.push_back(posT);
	    }
	}

	return res;
    }


    double
    MultipleAlignment::cmfinder_realignment_score(const MultipleAlignment &ma) const {
	size_t matches=0; // count matches in ma
	size_t exclusive_matches=0; // count matches in ma that are not found in reference
	
	//traverse all pairs of sequences in *this and ma
	for (size_type x=0; x<ma.alig_.size(); x++) {
	    const std::string &namex = ma.seqentry(x).name();
	    for (size_type y=x+1; y<ma.alig_.size(); y++) {
		const std::string &namey = ma.seqentry(y).name();
				
		// compute pairwise score
				
		matches += count_matches(ma.seqentry(x),ma.seqentry(y));

		exclusive_matches += count_exclusive_matches(ma.seqentry(x),ma.seqentry(y),
							     seqentry(namex),seqentry(namey)
							     );
		
	    }
	}
	
	// std::cout << exclusive_matches << " " << matches << " " << std::endl;
	
	return exclusive_matches / (double)matches;
    } 
    
    size_t
    MultipleAlignment::count_matches(const SeqEntry &a1,
				     const SeqEntry &a2) {

	size_t matches=0;

	size_t len_a=a1.seq().length();
	
	assert(len_a == a2.seq().length());

	for (size_t col=1; col<=len_a; col++) {
	    if ((!is_gap_symbol(a1.seq()[col]))
		&& (!is_gap_symbol(a2.seq()[col]))
		) {
		matches++;
	    }
	}
	
	return matches;
    }
    
    size_t
    MultipleAlignment::count_exclusive_matches(const SeqEntry &a1,
					       const SeqEntry &a2,
					       const SeqEntry &ref1,
					       const SeqEntry &ref2
					       ) {
	size_t matches=0;
	
	size_t len_a=a1.seq().length();
	size_t len_ref=ref1.seq().length();
	
	assert(len_a == a2.seq().length());
	assert(len_ref == ref2.seq().length());
	
	const std::vector<int> matchvecA = match_vector(a1.seq(),a2.seq());
	const std::vector<int> matchvecRef = match_vector(ref1.seq(),ref2.seq());
	
	// count exclusive matches alignments a
	size_t i=1;
	for (size_t col=1; col<=len_a; col++) {
	    if (!is_gap_symbol(a1.seq()[col])) {
		
		if (
		    (matchvecA[i]!=-1)
		    &&
		    (matchvecA[i]!=matchvecRef[i])) {
		    matches+=1;
		}
		i++;
	    }
	}
	
	return matches;
    }

    
    double
    MultipleAlignment::pairwise_match_score(const SeqEntry &a1,
					    const SeqEntry &a2,
					    const SeqEntry &ref1,
					    const SeqEntry &ref2,
					    bool score_common_gaps
					    ) {
	size_t s=0;
	
	size_t len_a=a1.seq().length();
	size_t len_ref=ref1.seq().length();
	
	assert(len_a == a2.seq().length());
	assert(len_ref == ref2.seq().length());
	
	const std::vector<int> matchvecA = match_vector(a1.seq(),a2.seq());
	const std::vector<int> matchvecRef = match_vector(ref1.seq(),ref2.seq());
	
	// count common matches in the two alignments
	size_t i=1;
	for (size_t col=1; col<=len_a; col++) {
	    if (!is_gap_symbol(a1.seq()[col])) {
		
		// add one to score, if position i of the first
		// sequence is matched to the same position by both
		// alignments
		if ((score_common_gaps || (matchvecA[i]!=-1))
		    && 
		    (matchvecA[i]==matchvecRef[i])) {
		    s+=1;
		}
		i++;
	    }
	}
	
	return (double)s;
    }
    
    double
    MultipleAlignment::pairwise_deviation_score(const SeqEntry &a1,
						const SeqEntry &a2,
						const SeqEntry &ref1,
						const SeqEntry &ref2
						) {
	double s=0.0;
	
	size_t len_a=a1.seq().length();
	size_t len_ref=ref1.seq().length();
	
	assert(len_a == a2.seq().length());
	assert(len_ref == ref2.seq().length());
	
	const std::vector<int> matchvecA = match_vector2(a1.seq(),a2.seq());
	const std::vector<int> matchvecRef = match_vector2(ref1.seq(),ref2.seq());
	
	// count common matches in the two alignments
	size_t i=1;
	for (size_t col=1; col<=len_a; col++) {
	    if (!is_gap_symbol(a1.seq()[col])) {
		
		double j_A = matchvecA[i] + ((matchvecA[i]==matchvecA[i-1])?0.5:0); 
		double j_Ref = matchvecRef[i] + ((matchvecRef[i]==matchvecRef[i-1])?0.5:0); 
		
		s += fabs(j_A-j_Ref);
		
		i++;
	    }
	}
	
	return s;
    }
    
    void
    MultipleAlignment::write_debug(ostream &out) const {
	for (size_type i=0; i<alig_.size(); ++i) {
	    out << alig_[i].name() << " \t" << alig_[i].seq().str()<<std::endl;
	} 
    
    }

    std::string
    MultipleAlignment::consensus_sequence() const {
	std::string cs="";
	//iterate over columns and built up consensus sequence 
	for (size_type i=1; i<=length(); ++i) {
	    
	    map<char,size_t> tab;
	    
	    // iterate over sequences and count character
	    for (std::vector<SeqEntry>::const_iterator it=alig_.begin(); alig_.end()!=it; ++it) {
		size_t c=it->seq()[i];
		if (tab.end()==tab.find(c)) tab[c]=0;
		tab[c]++;
	    }
	    
	    size_type cur_max=0;
	    char max_char='N';
	    for(map<char,size_t>::const_iterator it=tab.begin(); tab.end()!=it; ++it) {
		if (it->second > cur_max) {
		    cur_max=it->second;
		    max_char=it->first;
		}
	    }
	    
	    cs.push_back(max_char);
	}
	
	return cs;
    }

    std::ostream &
    MultipleAlignment::write_name_sequence_line(std::ostream &out,
						const std::string &name,
						const std::string &sequence) const{
	std::streamsize ow=out.width(26);
	out << std::left << name << " ";
	out.width(ow);
	out << sequence << std::endl;
	
	return out;
    }

    std::ostream &
    MultipleAlignment::write(std::ostream &out,
			     size_type start, 
			     size_type end) const
    {
	assert(1<=start);
	assert(end+1>=start);
	    	
	std::string structure_string = 
	    annotation(AnnoType::structure).single_string();
	
	std::string fixed_structure_string = 
	    annotation(AnnoType::fixed_structure).single_string();
	
	for (size_type i=0; i<alig_.size(); i++) {
	    const std::string seq = alig_[i].seq().str();
	    assert(end <= seq.length()); 
	    
	    write_name_sequence_line(out,alig_[i].name(),
				     seq.substr(start-1,end+1-start));
	}

	// write sequence anchors
	if (has_annotation(AnnoType::anchors)) {
	    	
	    for(size_t i=0; i<annotation(AnnoType::anchors).name_length(); i++) {
		const std::string anchors_tag   = "#"+annotation_tags[AnnoType::anchors];
		std::ostringstream name;
		name << anchors_tag << (i+1);
		
		std::string anchor_string = annotation(AnnoType::anchors).annotation_string(i);
		
		write_name_sequence_line(out,
					 name.str(),
					 anchor_string.substr(start-1,end-start+1));
	    }
	}
	
	if (has_annotation(AnnoType::structure)) {
	    const std::string structure_tag   = "#"+annotation_tags[AnnoType::structure];
	    write_name_sequence_line(out,
				     structure_tag,
				     structure_string.substr(start-1,end-start+1));
	}

	if (has_annotation(AnnoType::fixed_structure)) {
	    const std::string structure_tag   = "#"+annotation_tags[AnnoType::fixed_structure];
	    write_name_sequence_line(out,
				     structure_tag,
				     fixed_structure_string.substr(start-1,end-start+1));
	}
	
	return out;
    }

    std::ostream &
    MultipleAlignment::write(std::ostream &out) const {
	return
	    write(out,1,length());
    }

    std::ostream &
    MultipleAlignment::write(std::ostream &out,size_t width) const {
	size_t start=1;
	do {
	    size_t end = std::min(length(),start+width-1);
	    write(out,start,end);
	    start=end+1;
	} while (start <= length() && out << std::endl);
	
	return out;
    }
    
    void 
    MultipleAlignment::reverse() {
	for (std::vector<SeqEntry>::iterator it=alig_.begin(); alig_.end()!=it; ++it) {
	    it->reverse();
	}
    }
    
    bool 
    MultipleAlignment::checkAlphabet(const Alphabet<char> &alphabet) const {
	for (const_iterator it=begin(); end()!=it; ++it) {
	    for (size_type i=1; i<=it->seq().length(); i++) {
		if (!alphabet.in( it->seq()[i] )) {
		    return false;
		}
	    }
	}
	
    	return true;
    }

    void MultipleAlignment::append(const SeqEntry &seqentry) {
	name2idx_[seqentry.name()]=alig_.size(); // keep map name->index up to date 	
	alig_.push_back(seqentry);
    }

    void MultipleAlignment::prepend(const SeqEntry &seqentry) {
	alig_.insert(alig_.begin(),seqentry);
	
	name2idx_.clear();
	create_name2idx_map();
    }
    
    void
    MultipleAlignment::operator += (const AliColumn &c) {
	for (size_type i=0; i<alig_.size(); ++i) {
	    alig_[i].push_back(c[i]);
	}
    }

    void
    MultipleAlignment::operator += (char c) {    
	for (size_type i=0; i<alig_.size(); ++i) {
	    alig_[i].push_back(c);
	}
    }

    std::ostream &
    operator << (std::ostream &out, const MultipleAlignment &ma) {
	ma.write(out);
	return out;
    }
    

} // end namespace 
