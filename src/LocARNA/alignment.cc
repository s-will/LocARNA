#include <math.h>

#include "alignment.hh"
#include "alignment_impl.hh"

#include "sequence.hh"
#include "basepairs.hh"
#include "rna_data.hh"
#include "rna_structure.hh"
#include "anchor_constraints.hh"
#include "multiple_alignment.hh"
#include "string1.hh"
#include "plusvector.hh"

#ifdef HAVE_LIBRNA
extern "C" {
#  include <string.h>
#  include <ViennaRNA/data_structures.h>
#  include <ViennaRNA/fold_vars.h>
#  include <ViennaRNA/alifold.h>
}
#endif

namespace LocARNA {

    Alignment::Alignment(const Sequence &seqA, 
			 const Sequence &seqB)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	clear();
    }
    

    Alignment::Alignment(const Sequence &seqA, 
			 const Sequence &seqB, 
			 const std::string &alistrA,
			 const std::string &alistrB)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	
	assert(alistrA.length()==alistrB.length());
	
	size_t alilen = alistrA.length();
	
#     ifndef NDEBUG
	// assert that there is no gap only column in input
	bool no_gaponly=true;
	for(size_t k=0; k<alilen && no_gaponly; ++k) {
	    if (is_gap_symbol(alistrA[k]) && is_gap_symbol(alistrB[k])) {
		no_gaponly=false;
	    }
	}
	assert(no_gaponly);

	// assert correct number of gaps
	size_t gapsA=0;
	size_t gapsB=0;
	for(size_t k=0; k<alilen; ++k) {
	    if (is_gap_symbol(alistrA[k])) {
		gapsA++;
	    }
	    if (is_gap_symbol(alistrB[k])) {
		gapsB++;
	    }
	}
	assert(alistrA.length() == seqA.length() + gapsA);
	assert(alistrB.length() == seqB.length() + gapsB);
#     endif

	
	size_t posA=0;
	size_t posB=0;
	
	for(size_t k=0; k<alilen; ++k) {
	    if (is_gap_symbol(alistrA[k])) {
		// in case of locality gaps, no alignment edge is appended;
		// skip directly to the next position in the alignment strings
		if (alistrA[k]=='~') continue;
		
		// add gap entry
		pimpl_->a_.push_back(-1);
	    }
	    else {
		posA++;
		pimpl_->a_.push_back(posA);
	    }

	    if (is_gap_symbol(alistrB[k])) {
		// in case of locality gaps, no alignment edge is appended;
		// skip directly to the next position in the alignment strings
		if (alistrB[k]=='~') continue;
		pimpl_->b_.push_back(-1);
	    } else {
		posB++;
		pimpl_->b_.push_back(posB);
	    }
	}
    }

    Alignment::Alignment(const Alignment &alignment)
	: pimpl_(new AlignmentImpl(*alignment.pimpl_))
    {
    }
    
    Alignment::~Alignment() { delete pimpl_; }

    Alignment &
    Alignment::operator =(const Alignment &alignment)
    {
	pimpl_ = new AlignmentImpl(*alignment.pimpl_);
	return *this;
    }
    
    
    void
    Alignment::set_structures(const RnaStructure &structureA,const RnaStructure &structureB) {
	pimpl_->strA_ = structureA.to_string();
	pimpl_->strB_ = structureB.to_string();
	assert(pimpl_->strA_.length() == pimpl_->seqA_.length());
	assert(pimpl_->strB_.length() == pimpl_->seqB_.length());
    }

    void
    Alignment::set_consensus_structure(const RnaStructure &structure) {
	set_structures(structure,structure);
    }

    void Alignment::clear() {
	pimpl_->strA_.resize(pimpl_->seqA_.length()+1);
	pimpl_->strB_.resize(pimpl_->seqB_.length()+1);
	fill(pimpl_->strA_.begin(),pimpl_->strA_.end(),'.');
	fill(pimpl_->strB_.begin(),pimpl_->strB_.end(),'.');
	
	pimpl_->a_.clear();
	pimpl_->b_.clear();
    }

    void 
    Alignment::append(int i, int j) {
	pimpl_->a_.push_back(i);
	pimpl_->b_.push_back(j);
    }

    void
    Alignment::add_basepairA(int i, int j) {
	pimpl_->strA_[i]='(';
	pimpl_->strA_[j]=')';
    }
    
    void
    Alignment::add_basepairB(int i, int j) {
	pimpl_->strB_[i]='(';
	pimpl_->strB_[j]=')';
    }
    
    
    const Alignment::edge_vector_t 
    Alignment::alignment_edges(bool only_local) const {
	edge_vector_t edges;
	
	const std::vector<int> &a=pimpl_->a_;
	const std::vector<int> &b=pimpl_->b_;

	size_t alisize=a.size();

	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B
	
	for (size_type i=0; i<alisize; i++) {
	    for (int j=lastA; j<a[i]; j++) {
		if (!only_local) edges.push_back(edge_t(j,-2));
		lastA++;
	    }
	    for (int j=lastB; j<b[i]; j++) {
		if (!only_local) edges.push_back(edge_t(-2,j));
		lastB++;
	    }
	    int x=-1;
	    int y=-1;
	    if ( a[i] >= 0 ) {
		x = a[i];
		lastA++;
	    }
	    if ( b[i] >= 0 ) {
		y = b[i];
		lastB++;
	    }
	    edges.push_back(edge_t(x,y));
	}
	
	if (!only_local) { 
	    for (size_type j=lastA; j<=pimpl_->seqA_.length(); j++) {
		edges.push_back(edge_t(j,-2));
		lastA++;
	    }
	    for (size_type j=lastB; j<=pimpl_->seqB_.length(); j++) {
		edges.push_back(edge_t(-2,j));
		lastB++;
	    }
	}
	return edges;
    }
    
    
    std::string
    AlignmentImpl::consensus_constraint_string(const AnchorConstraints &seqConstraints, 
					       const Alignment::edge_vector_t &edges) {
	std::string cons_constraint="";
	if (!seqConstraints.empty()) {
	    
	    std::vector<std::string> seqCStrings(seqConstraints.name_size());

	    std::string noname="";
	    for (size_type j=0; j<seqConstraints.name_size(); ++j) noname+=".";
	    
	    for (size_type i=0; i<edges.size(); i++) {
		const std::string &nameA = (edges[i].first>0) ?seqConstraints.get_name_a(edges[i].first) :"";
		const std::string &nameB = (edges[i].second>0)?seqConstraints.get_name_b(edges[i].second):"";
		
		std::string name=noname;

		if (nameA!="") name =  nameA;
		else if (nameB!="") name = nameB;

		for (size_type j=0; j<seqConstraints.name_size(); ++j) {
		    seqCStrings[j] += name[j];
		}
	    }

	    for (size_type j=0; j<seqConstraints.name_size(); ++j) {
		if (j!=0) cons_constraint += "#";
		cons_constraint += seqCStrings[j];
	    }
	}
	return cons_constraint;
    }


    size_type
    Alignment::local_startA() const {return pimpl_->a_[0];}
    
    size_type
    Alignment::local_endA() const {return pimpl_->a_[pimpl_->a_.size()-1];}
    
    size_type
    Alignment::local_startB()  const {return pimpl_->b_[0];}
	
    size_type
    Alignment::local_endB() const {return pimpl_->b_[pimpl_->b_.size()-1];}

    const Sequence &
    Alignment::seqA() const {return pimpl_->seqA_;} 

    const Sequence &
    Alignment::seqB() const {return pimpl_->seqB_;} 

    // const std::vector<int> &
    // Alignment::get_a() const {return pimpl_->a_;} 

    // const std::vector<int> &
    // Alignment::get_b() const {return pimpl_->b_;} 


    void AlignmentImpl::write_debug(std::ostream &out) const {

	for (size_type i=0; i<a_.size(); i++) {
	    out << a_[i] << " ";
	}
	out<<std::endl;
	for (size_type i=0; i<b_.size(); i++) {
	    out << b_[i] << " ";
	}
	out<<std::endl;

    }


    std::string
    Alignment::dot_bracket_structureA(bool only_local) const {
	edge_vector_t edges = alignment_edges(only_local);
	return pimpl_->dot_bracket_structure(pimpl_->strA_,edges,1);
    }

    std::string
    Alignment::dot_bracket_structureB(bool only_local) const {
	edge_vector_t edges = alignment_edges(only_local);
	return pimpl_->dot_bracket_structure(pimpl_->strB_,edges,2);
    }

    std::string
    AlignmentImpl::dot_bracket_structure(const std::string &str,
					 const Alignment::edge_vector_t &x,
					 size_t ab) {
	std::string s;
	for (size_t i=0; i<x.size(); i++) {
	    int xi=(ab==1)?x[i].first:x[i].second;
	    if (xi>0) {
		s.push_back(str[xi]);
	    } else if (xi==-1) {
		s.push_back('-');
	    } else if (xi==-2) {
		s.push_back('~');
	    }
	}
	return s;
    }
    
}
