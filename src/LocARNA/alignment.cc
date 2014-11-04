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

extern "C" {
#  include <string.h>
#  include <ViennaRNA/data_structures.h>
#  include <ViennaRNA/fold_vars.h>
#  include <ViennaRNA/alifold.h>
}

namespace LocARNA {

    Alignment::Alignment(const Sequence &seqA, 
			 const Sequence &seqB)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	clear();
    }
    

    Alignment::Alignment(const Sequence &seqA, 
    			 const Sequence &seqB, 
			 const edges_t &edges)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	
	// append all non-locality edges
	for (size_t k=0; k<edges.size(); k++) {
	    const edge_end_t &x = edges.first[k];
	    const edge_end_t &y = edges.second[k];
	    
	    if (x.is_gap() && y.is_gap()) {
		throw failure("Invalid alignment edges");
	    }
	    if (x.is_pos() && (x<1 || x>seqA.length())) {
		throw failure("Alignment edge out of range (first sequence).");
	    }
	    if (y.is_pos() && (y<1 || y>seqB.length())) {
		throw failure("Alignment edge out of range (second sequence).");
	    }
	    
	    if (!((x.is_gap() && x.gap()==Gap::locality)
		  ||
		  (y.is_gap() && y.gap()==Gap::locality))) {
		append(x,y);
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
    Alignment::append(edge_end_t i, edge_end_t j) {
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
    
    
	void 
	Alignment::add_deleted_basepairA(int i, int j) {
	pimpl_->strA_[i]='(';
	pimpl_->strA_[j]=')';
	}

	void 
	Alignment::add_deleted_basepairB(int i, int j) {
	pimpl_->strB_[i]='(';
	pimpl_->strB_[j]=')';
	}

    const Alignment::edges_t 
    Alignment::alignment_edges(bool only_local) const {

	edge_ends_t endsA;
	edge_ends_t endsB;
	
	const edge_ends_t &a=pimpl_->a_;
	const edge_ends_t &b=pimpl_->b_;

	size_t alisize=a.size();

	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B
	
	for (size_type i=0; i<alisize; i++) {
	    if (a[i].is_pos()) {
		for (size_t j=lastA; j<a[i]; j++) {
		    if (!only_local) {
			endsA.push_back(j);
			endsB.push_back(Gap::locality);
		    }
		    lastA++;
		}
	    }
	    if (b[i].is_pos()) {
		for (size_t j=lastB; j<b[i]; j++) {
		    if (!only_local) {
			endsA.push_back(Gap::locality);
			endsB.push_back(j);
		    }
		    lastB++;
		}
	    }
	    if ( a[i].is_pos() ) {
		lastA++;
	    }
	    if ( b[i].is_pos() ) {
		lastB++;
	    }
	    endsA.push_back(a[i]);
	    endsB.push_back(b[i]);
	}
	
	if (!only_local) { 
	    for (size_type j=lastA; j<=pimpl_->seqA_.length(); j++) {
		endsA.push_back(j);
		endsB.push_back(Gap::locality);
		lastA++;
	    }
	    for (size_type j=lastB; j<=pimpl_->seqB_.length(); j++) {
		endsA.push_back(Gap::locality);
		endsB.push_back(j);
		lastB++;
	    }
	}

	return edges_t(endsA,endsB);
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


    void
    AlignmentImpl::write_debug(std::ostream &out, const Alignment::edge_ends_t &ends) {
	for (size_type i=0; i<ends.size(); i++) {
	    if (ends[i].is_pos()) {
		out << ends[i] << " ";
	    } else {
		out << "g" << ends[i].gap().idx() << " ";
	    }
	}
	out<<std::endl;
    }
    
    void
    AlignmentImpl::write_debug(std::ostream &out) const {
	write_debug(out,a_);
	write_debug(out,b_);
    }


    std::string
    Alignment::dot_bracket_structureA(bool only_local) const {
	edges_t edges = alignment_edges(only_local);
	return pimpl_->dot_bracket_structure(pimpl_->strA_,edges.first);
    }

    std::string
    Alignment::dot_bracket_structureB(bool only_local) const {
	edges_t edges = alignment_edges(only_local);
	return pimpl_->dot_bracket_structure(pimpl_->strB_,edges.second);
    }

    std::string
    AlignmentImpl::dot_bracket_structure(const std::string &str,
					 const Alignment::edge_ends_t &x) {
	std::string s;
	for (size_t i=0; i<x.size(); i++) {
	    if (x[i].is_pos()) {
		s.push_back(str[x[i]]);
	    } else if (x[i].is_gap()) {
		s.push_back(gap_symbol(x[i].gap()));
	    }
	}
	return s;
    }

    Alignment::edge_ends_t 
    Alignment::alistr_to_edge_ends(const std::string alistr)
    {
	edge_ends_t xs;
	size_t i=1;
	
	for (size_t k=0; k<alistr.size(); k++) {
	    edge_end_t x;
	    if (is_gap_symbol(alistr[k])) {
		x = gap_code(alistr[k]);
	    } else {
		x = i;
		i++;
	    }
	    xs.push_back(x);
	}
	return xs;
    }

}
