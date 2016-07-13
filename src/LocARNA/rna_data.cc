#include <stdio.h>
#include <ctype.h> // import isspace

#include <math.h> // import log

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "aux.hh"
#include "pfold_params.hh"
#include "alignment.hh"
#include "rna_ensemble.hh"
#include "sequence_annotation.hh"
#include "rna_data_impl.hh"
#include "ext_rna_data_impl.hh"
#include "rna_structure.hh"

#include "LocARNA/global_stopwatch.hh"

extern "C" {
#   include <ViennaRNA/energy_const.h> // import TURN
#   include <ViennaRNA/MEA.h>
}


namespace LocARNA {

    RnaData::RnaData(const RnaEnsemble &rna_ensemble,
		     double p_bpcut,
		     double max_bps_length_ratio,
		     const PFoldParams &pfoldparams)
	: pimpl_(new RnaDataImpl(this,
				 p_bpcut)) {
	init_from_rna_ensemble(rna_ensemble,pfoldparams);

	if (max_bps_length_ratio > 0) {
	    pimpl_->drop_worst_bps(max_bps_length_ratio*pimpl_->sequence_.length());
	}
    }


    RnaData::RnaData(const std::string &filename,
		     double p_bpcut,
		     double max_bps_length_ratio,
		     const PFoldParams &pfoldparams
		     )
	: pimpl_(new RnaDataImpl(this,
				 p_bpcut
				 )) {
	bool complete=
	    read_autodetect(filename,
			    pfoldparams.stacking());
    	
	if (!complete) {
	    // recompute all probabilities
	    RnaEnsemble 
		rna_ensemble(pimpl_->sequence_,
			     pfoldparams,false,true); // use given parameters, no in loop, use alifold
	    
	    // initialize from RnaEnsemble; note: method is virtual
	    init_from_rna_ensemble(rna_ensemble,
				   pfoldparams);
	}
	
	if (max_bps_length_ratio > 0) {
	    pimpl_->drop_worst_bps(max_bps_length_ratio*pimpl_->sequence_.length());
	}

    }
    
    // do almost nothing
    RnaData::RnaData(double p_bpcut)
	: pimpl_(new RnaDataImpl(this,
				 p_bpcut)) {
    }
    
    // "consensus" constructor
    RnaData::RnaData(const RnaData &rna_dataA,
		     const RnaData &rna_dataB,
		     const Alignment &alignment,
		     double p_expA,
		     double p_expB,
		     bool only_local
		     )
	: pimpl_(new RnaDataImpl(this,
				 rna_dataA,
				 rna_dataB,
				 alignment.alignment_edges(only_local),
				 p_expA,
				 p_expB)) {
    }
    
    RnaData::~RnaData() {
	delete pimpl_;
    }

    bool
    RnaData::has_stacking() const {
	return pimpl_->has_stacking_;
    }

    void
    RnaData::set_anchors(const SequenceAnnotation &anchors) {
	pimpl_->sequence_.set_annotation(MultipleAlignment::AnnoType::anchors,anchors);
    }

    // "consensus" constructor
    RnaDataImpl::RnaDataImpl(RnaData *self,
			     const RnaData &rna_dataA,
			     const RnaData &rna_dataB,
			     const Alignment::edges_t &edges,
			     double p_expA,
			     double p_expB) 
	:self_(self),
	 sequence_(edges,rna_dataA.sequence(),rna_dataB.sequence()),
	 p_bpcut_(),
	 arc_probs_(0.0),
	 arc_2_probs_(0.0),
	 has_stacking_(false)
    {
	init_as_consensus_dot_plot(edges,
				   rna_dataA,
				   rna_dataB,
				   p_expA,
				   p_expB,
				   rna_dataA.has_stacking() && rna_dataB.has_stacking());
    }


    // do almost nothing
    RnaDataImpl::RnaDataImpl(RnaData *self,
			     double p_bpcut)
	:self_(self),
	 sequence_(),
	 p_bpcut_(p_bpcut),
	 arc_probs_(0.0),
	 arc_2_probs_(0.0),
	 has_stacking_(false)
    {
    }
    

    ExtRnaData::ExtRnaData(const std::string &filename,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut,
			   double max_bps_length_ratio,
			   double max_uil_length_ratio,
			   double max_bpil_length_ratio,
			   const PFoldParams &pfoldparams)
	: 
	RnaData(p_bpcut),
	ext_pimpl_(new ExtRnaDataImpl(this,
                                      p_bpilcut,
                                      p_uilcut)) {

	bool complete=
	    read_autodetect(filename,pfoldparams.stacking());
    	
	if (!complete) {
	    // recompute all probabilities
	    RnaEnsemble 
		rna_ensemble(sequence(),
			     pfoldparams,true,true); // use given parameters, in-loop, use alifold
	    
	    // initialize
	    init_from_rna_ensemble(rna_ensemble,pfoldparams);
	}
	
	if (max_bps_length_ratio > 0) {
	    ext_pimpl_->drop_worst_bps(max_bps_length_ratio*length());
	}
	if (max_uil_length_ratio > 0) {
	    ext_pimpl_->drop_worst_uil(max_uil_length_ratio*length());
	}
	if (max_bpil_length_ratio > 0) {
	    ext_pimpl_->drop_worst_bpil_precise(max_bpil_length_ratio);
	}

    }

    ExtRnaDataImpl::ExtRnaDataImpl(ExtRnaData *self,
				   double p_bpilcut,
				   double p_uilcut)
	:self_(self),
	 p_bpilcut_(p_bpilcut),
	 p_uilcut_(p_uilcut),
	 arc_in_loop_probs_(arc_prob_matrix_t(0.0)),
	 unpaired_in_loop_probs_(arc_prob_vector_t(0.0)),
	 has_in_loop_probs_(false)
    {
    }

    ExtRnaData::ExtRnaData(const RnaEnsemble &rna_ensemble,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut,
			   double max_bps_length_ratio,
			   double max_uil_length_ratio,
			   double max_bpil_length_ratio,
			   const PFoldParams &pfoldparams)
	:
	RnaData(rna_ensemble,
		p_bpcut,
		max_bps_length_ratio,
		pfoldparams),
	ext_pimpl_(new ExtRnaDataImpl(this,
				  p_bpilcut,
				  p_uilcut)) {
	init_from_rna_ensemble(rna_ensemble,pfoldparams);

	if (max_uil_length_ratio > 0) {
	    	ext_pimpl_->drop_worst_uil(max_uil_length_ratio*length());
	}
	if (max_bpil_length_ratio > 0) {
	    	ext_pimpl_->drop_worst_bpil(max_bpil_length_ratio*length());
	}
	
    }

    ExtRnaData::~ExtRnaData() {
	delete ext_pimpl_;
    }

    bool
    RnaData::read_autodetect(const std::string &filename,
			     bool stacking) {
	bool failed=true;  //flag for signalling a failed attempt to
			   //read a certain file format
	
	bool sequence_only=false; 
	// does the file format contain only sequence information or
	// should we assume that the base pair probabilities are given

	pimpl_->has_stacking_=stacking;

	// try dot plot ps format
	if (failed) {
	    sequence_only=false;
	    failed=false;
	    try {
		read_ps(filename);
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
	    } catch (wrong_format_failure &f) {
		failed=true;
	    }
	}

	// try pp 2.0
	if (failed) {
	    sequence_only=false;
	    failed=false;
	    try {
		// std::cerr << "Try reading pp "<<filename<<" ..."<<std::endl;
		read_pp(filename);
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
	    } catch (wrong_format_failure &f) {
		failed=true;
	    }
	    // if (failed) std::cerr << "  ... did not succeed."<<std::endl;
	    // else std::cerr << "  ... success."<<std::endl;
	}
	

	// try fasta format
	if (failed) {
	    sequence_only = true;
	    failed=false;
	    try {
		MultipleAlignment ma(filename, MultipleAlignment::FormatType::FASTA);
		pimpl_->sequence_ = ma;
		// even if reading does not fail, we still want to
		// make sure that the result is reasonable. Otherwise,
		// we assume that the file is in a different format.
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
	    } catch (failure &f) {
		failed=true;
	    }
	}
	
	// try old pp
	if (failed) {
	    sequence_only=false;
	    failed=false;
	    try {
		// std::cerr << "Try reading old pp "<<filename<<" ..."<<std::endl;
		read_old_pp(filename);
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
	    } catch (wrong_format_failure &f) {
		failed=true;
	    }
	    // if (failed) std::cerr << "  ... did not succeed."<<std::endl;
	    // else std::cerr << "  ... success."<<std::endl;
	}
        
        // try stockholm format
        if (failed) {
            sequence_only=true;
	    failed=false;
	    try {
		// std::cerr << "Try reading stockholm "<<filename<<" ..."<<std::endl;
		MultipleAlignment ma(filename, MultipleAlignment::FormatType::STOCKHOLM);
                
		pimpl_->sequence_ = ma;
	    	// even if reading does not fail, we still want to
		// make sure that the result is reasonable. Otherwise,
		// we assume that the file is in a different format.
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}

                // handle potential fixed structure in MA
		if (!failed) {
                    typedef MultipleAlignment::AnnoType TA;

		    if (ma.has_annotation(TA::fixed_structure)) {
			init_from_fixed_structure(ma.annotation(TA::fixed_structure),
						  stacking);
			sequence_only=false;
		    }
		}

            } catch (wrong_format_failure &f) {
                failed=true;
            }catch (syntax_error_failure &f) {
		throw failure((std::string)"RnaData: Cannot read input data from stockholm file.\n\t"+f.what());
	    } catch (failure &f) {
		failed=true;
	    }
        }

	// try clustal format
	if (failed) {
	    sequence_only=true;
	    failed=false;
	    try {
		// std::cerr << "Try reading clustal "<<filename<<" ..."<<std::endl;
		MultipleAlignment ma(filename, MultipleAlignment::FormatType::CLUSTAL);

		pimpl_->sequence_ = ma;
	    	// even if reading does not fail, we still want to
		// make sure that the result is reasonable. Otherwise,
		// we assume that the file is in a different format.
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
		
		// handle potential fixed structure in MA
		if (!failed) {
                    typedef MultipleAlignment::AnnoType TA;

		    if (ma.has_annotation(TA::fixed_structure)) {
			init_from_fixed_structure(ma.annotation(TA::fixed_structure),
						  stacking);
			sequence_only=false;
		    }
		}
		
	    } catch (syntax_error_failure &f) {
		throw failure((std::string)"RnaData: Cannot read input data from clustal file.\n\t"+f.what());
	    }
	    catch (failure &f) {
		failed=true;
	    }
	    // if (failed) std::cerr << "  ... did not succeed."<<std::endl;
	    // else std::cerr << "  ... success."<<std::endl;
	}
	
	if (failed) {
	    throw failure("RnaData: Cannot read input data from file.");
	}
	
	pimpl_->sequence_.normalize_rna_symbols();

	// now, we have the sequence but not necessarily all required probabilities!
	// return whether (re)computation of probabilities is required
	return !sequence_only && inloopprobs_ok();
    }

    void
    RnaData::init_from_fixed_structure(const SequenceAnnotation &structure,
				       bool stacking) {
	pimpl_->init_from_fixed_structure(structure,stacking);
    }
    
    void
    RnaDataImpl::init_from_fixed_structure(const SequenceAnnotation &structure,
					   bool stacking) {
	assert(structure.length() == sequence_.length());
	RnaStructure rna_structure(structure.single_string());
	
	p_bpcut_=1.0;
	
	for (RnaStructure::const_iterator it=rna_structure.begin();
	     rna_structure.end() != it; ++it) {
	    arc_probs_(it->first,it->second)=1.0;
	
	    if (stacking) {
		if (rna_structure.contains(RnaStructure::bp_t(it->first+1,it->second-1))) {
		    arc_2_probs_(it->first,it->second)=1.0;
		}
	    }
	}
	has_stacking_=stacking;
    }

    void
    ExtRnaData::init_from_fixed_structure(const SequenceAnnotation &structure,
					  bool stacking) {
	RnaData::init_from_fixed_structure(structure,stacking);
	ext_pimpl_->init_from_fixed_structure(structure);
    }

    void
    ExtRnaDataImpl::init_fixed_unpaired_in_loop(size_t i,
						size_t j,
						const RnaStructure &rna_structure) {
	for (size_t k=i+1; k<j; ++k) {
	    bool contained=true;
	    for (RnaStructure::const_iterator it2=rna_structure.begin();
		 contained && rna_structure.end() != it2; ++it2) {
		if (i < it2->first
		    && it2->first <= k
		    && k <= it2->second
		    && it2->second < j) {
		    
		    contained=false;
		    // k=it2->second-1; // micro-optimization
		}
	    }
	    if (contained) {
		std::cout << "UIL "<<i<<" "<<j<<": "<<k<<std::endl;
		unpaired_in_loop_probs_.ref(i,j)[k] = 1.0;
	    }
	    
	}
    }

    void
    ExtRnaDataImpl::init_fixed_basepairs_in_loop(size_t i,
						 size_t j,
						 const RnaStructure &rna_structure) {
	for (RnaStructure::const_iterator it2=rna_structure.begin();
	     rna_structure.end() != it2; ++it2) {
	    
	    bool contained=true;
	    
	    if (it2->first <= i
		|| j <= it2->second) {
		contained=false;
	    }
	    for (RnaStructure::const_iterator it3=rna_structure.begin();
		 contained && rna_structure.end() != it3; ++it3) {
		if (i < it3->first
		    && it3->first < it2->first
		    && it2->second < it3->second
		    && it3->second < j) {
		    
		    contained=false;
		}
	    }
	    if (contained) {
		std::cout << "BPIL "<<i<<" "<<j<<" "
			  << ": "<<it2->first<<" "<<it2->second<<" "
			  << std::endl;
		arc_in_loop_probs_.ref(i,j)(it2->first,it2->second) = 1.0;
	    }
	    
	}
    }
    

    void
    ExtRnaDataImpl::init_from_fixed_structure(const SequenceAnnotation &structure) {
	RnaStructure rna_structure(structure.single_string());
	
	// initialize in loop probabilities
	//
	// for each base pair, we enumerate the bases and base pairs
	// that are contained in its loop.
	//
	// NOTE: we implement a very inefficient approach!  note
	// that this seems acceptable (for a start), since -as generated here-
	// rna_structure contains only a linear number of base pairs
	// (which limits run-time to cubic).
	//
	
	for (RnaStructure::const_iterator it=rna_structure.begin();
	     rna_structure.end() != it; ++it) {
	    
	    init_fixed_unpaired_in_loop(it->first,it->second,rna_structure);

	    init_fixed_basepairs_in_loop(it->first,it->second,rna_structure);
	}

	// external loop

	init_fixed_unpaired_in_loop(0,rna_structure.length()+1,rna_structure);
	
	init_fixed_basepairs_in_loop(0,rna_structure.length()+1,rna_structure);
	
	
    }

    void
    RnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
				    const PFoldParams &pfoldparams) {
	pimpl_->init_from_rna_ensemble(rna_ensemble,pfoldparams);
    }

    void
    ExtRnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
				       const PFoldParams &pfoldparams) {
	RnaData::init_from_rna_ensemble(rna_ensemble,pfoldparams);
	ext_pimpl_->init_from_ext_rna_ensemble(rna_ensemble);
    }

    void
    RnaDataImpl::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
					const PFoldParams &pfoldparams) {
	assert(rna_ensemble.has_base_pair_probs());
	
	// ----------------------------------------
	// init sequence	
	sequence_ = rna_ensemble.multiple_alignment();
	size_t len = sequence_.length();

	// ----------------------------------------
	// init base pair probabilities
	arc_probs_.clear();
	for( size_t i=1; i <= len; i++ ) {
	    for( size_t j=i+TURN+1; j <= len; j++ ) {
		
		double p = rna_ensemble.arc_prob(i,j);
		
		// commented out: drop lonely base pairs
		// double p_out = 
		//     (i>1 && j<len)
		//     ? rna_ensemble.arc_prob(i-1,j+1)
		//     : 0.0;
		// double p_in = rna_ensemble.arc_prob(i+1,j-1);
		
		if ( p > p_bpcut_ ) { // apply filter
		    
		    // if (!pfoldparams.noLP() 
		    // 	|| p_in>p_bpcut_ 
		    // 	|| p_out>p_bpcut_ ) {
		    arc_probs_(i,j)=p;
		    // }  else {
		    //std::cout << "Drop "<< i<<" " <<j << " "<<p<<" "<<p_in<<" " <<p_out<<std::endl;
		    //}
		}
	    }
	}
	
	// ----------------------------------------
	// init stacking probabilities
	arc_2_probs_.clear();
	has_stacking_ = pfoldparams.stacking();	
	if (has_stacking_) {
	    for( size_t i=1; i <= len; i++ ) {
		for( size_t j=i+TURN+3; j <= len; j++ ) {
		    double p2 = rna_ensemble.arc_2_prob(i,j);
		    if (p2 > p_bpcut_) { // apply filter to joint probability !
			arc_2_probs_(i,j)=p2;
		    }
		}
	    }
	}

	// all set
	return;
    }

    void
    ExtRnaDataImpl::init_from_ext_rna_ensemble(const RnaEnsemble &rna_ensemble) {
	// initialize in loop probabilities
	// (usually, this is called after RnaDataImpl::init_from_rna_ensemble)
	assert(rna_ensemble.has_in_loop_probs());

	size_t len = self_->length();
	
	// ----------------------------------------
	// init base pair probabilities
	arc_in_loop_probs_.clear();
	
	// ------------------------------
	// construct helper data structure for efficiency:
	// map left ends to right ends of all arcs in arc_probs_
	std::vector<std::vector<size_t> > right_ends;
	right_ends.resize(len+1);
	for(arc_prob_matrix_t::const_iterator it = self_->arc_probs_begin();
	    self_->arc_probs_end()!=it; ++it) {
	    pos_type i = it->first.first;
	    pos_type j = it->first.second;
	    right_ends[i].push_back(j);
	}
	for(std::vector<std::vector<size_t> >::iterator it = right_ends.begin();
	    right_ends.end()!=it; ++it) {
	    sort(it->begin(),it->end());
	}
	// end constructing helper data structure

	// in loop
	for(arc_prob_matrix_t::const_iterator it = self_->arc_probs_begin();
	    self_->arc_probs_end()!=it; ++it) {
	    pos_type i = it->first.first;
	    pos_type j = it->first.second;
	    arc_prob_matrix_t m_ij(0.0);
	    
	    for(size_t ip=i+1; ip < j; ip++ ) {
		// for( size_t jp=ip+TURN+1; jp < j; jp++ ) {
		for(std::vector<size_t>::const_iterator jpit = right_ends[ip].begin();
		    right_ends[ip].end()!=jpit && *jpit<j; ++jpit) {
		    size_t jp = *jpit;

		    double p = rna_ensemble.arc_in_loop_prob(ip,jp,i,j);
		    
	    	    if ( p > p_bpilcut_ ) {
	    		m_ij(ip,jp)=p;
	    	    }
	    	}
	    }

	    
	    // set only if not empty; use set instead of assignment,
	    // to avoid the comparison of complex SparseMatrix objects
	    if (!m_ij.empty()) {
		arc_in_loop_probs_.set(i,j,m_ij);
	    }
	}

	// external
	arc_prob_matrix_t m_ext(0.0);
	for( size_t ip=1; ip < len; ip++ ) {
	    for(std::vector<size_t>::const_iterator jpit = right_ends[ip].begin();
		right_ends[ip].end()!=jpit; ++jpit) {
		size_t jp = *jpit;

		double p = rna_ensemble.arc_external_prob(ip,jp);
		
		if ( p > p_bpilcut_ ) {
		    m_ext(ip,jp) = p;
		}
	    }
	}

	// set only if not empty; see above
	if (!m_ext.empty()) {
	    arc_in_loop_probs_.set(0,self_->length()+1,m_ext);
	}
	
	// ----------------------------------------
	// init unpaired probabilities
	unpaired_in_loop_probs_.clear();
	
	// in loop
	for(arc_prob_matrix_t::const_iterator it = self_->arc_probs_begin();
	    self_->arc_probs_end()!=it; ++it) {
	    pos_type i = it->first.first;
	    pos_type j = it->first.second;
	    arc_prob_vector_t v_ij(0.0);
	    
	    for( size_t k=i+1; k < j; k++ ) {
		double p = rna_ensemble.unpaired_in_loop_prob(k,i,j);
		if ( p > p_uilcut_ ) {
		    v_ij[k] = p;
		}
	    }
	    
	    // set only if not empty; see above
	    if (!v_ij.empty()) {
		unpaired_in_loop_probs_.set(i,j,v_ij);
	    }
	}

	// external
	arc_prob_vector_t v_ext(0.0);
	for( size_t k=1; k <= len; k++ ) {
	    double p = rna_ensemble.unpaired_external_prob(k);
	    
	    if ( p > p_uilcut_ ) {
		v_ext[k] = p;
	    }
	}

	// set only if not empty; see above
	if (!v_ext.empty()) {
	    unpaired_in_loop_probs_.set(0,self_->length()+1,v_ext);
	}
	
	
	// set flag
	has_in_loop_probs_=true;

	// all set
	return;
    } // end method init_from_ext_rna_ensemble

    bool
    ExtRnaData::inloopprobs_ok() const {
	return ext_pimpl_->has_in_loop_probs_;
    }

    const Sequence &
    RnaData::sequence() const {
	return pimpl_->sequence_.as_sequence();
    }

    const MultipleAlignment &
    RnaData::multiple_alignment() const {
	return pimpl_->sequence_;
    }

    size_type
    RnaData::length() const {
	return pimpl_->sequence_.length();
    }
    
    double
    RnaData::arc_cutoff_prob() const {
	return pimpl_->p_bpcut_;
    }
    
    double 
    RnaData::arc_prob(pos_type i, pos_type j) const {
	return pimpl_->arc_probs_(i,j);
    }
    
    RnaData::arc_probs_const_iterator
    RnaData::arc_probs_begin() const {
	return pimpl_->arc_probs_.begin();
    }

    RnaData::arc_probs_const_iterator
    RnaData::arc_probs_end() const {
	return pimpl_->arc_probs_.end();
    }


    double
    RnaData::joint_arc_prob(pos_type i, pos_type j) const {
	return pimpl_->arc_2_probs_(i,j);
    }
    
    double 
    RnaData::stacked_arc_prob(pos_type i, pos_type j) const {
	assert(pimpl_->arc_probs_(i+1,j-1)!=0);
	
	return
	    pimpl_->arc_2_probs_(i,j)
	    /
	    pimpl_->arc_probs_(i+1,j-1);
    }

    double 
    RnaData::prob_paired_upstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=i+1; j<=length(); j++) {
	    prob_paired += pimpl_->arc_probs_(i,j); 
	}
	
	return prob_paired;
    }
    
    double
    RnaData::prob_paired_downstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=1; j<i; j++) {
	    prob_paired += pimpl_->arc_probs_(j,i); 
	}
	
	return prob_paired;
    }
    
    double
    RnaData::prob_unpaired(size_type i) const {
	return 
	    1.0
	    - prob_paired_upstream(i)
	    - prob_paired_downstream(i);
    }
    
    double
    ExtRnaData::arc_in_loop_cutoff_prob() const {
	return ext_pimpl_->p_bpilcut_;
    }
	
    double 
    ExtRnaData::arc_in_loop_prob(pos_type i, pos_type j,pos_type p, pos_type q) const {
	ExtRnaDataImpl::arc_prob_matrix_t m_pq = ext_pimpl_->arc_in_loop_probs_(p,q);
	return m_pq(i,j);
    }
    
    double 
    ExtRnaData::arc_external_prob(pos_type i, pos_type j) const {
	ExtRnaDataImpl::arc_prob_matrix_t m_ext = ext_pimpl_->arc_in_loop_probs_(0,length()+1);
	return m_ext(i,j);
    }
    
    double
    ExtRnaData::unpaired_in_loop_cutoff_prob() const {
	return ext_pimpl_->p_uilcut_;
    }
    
    double 
    ExtRnaData::unpaired_in_loop_prob(pos_type k,pos_type p, pos_type q) const {
	ExtRnaDataImpl::arc_prob_vector_t v_pq = ext_pimpl_->unpaired_in_loop_probs_(p,q);
	return v_pq[k];
    }
    
    double 
    ExtRnaData::unpaired_external_prob(pos_type k) const {
	ExtRnaDataImpl::arc_prob_vector_t v_ext = ext_pimpl_->unpaired_in_loop_probs_(0,length()+1);
	return v_ext[k];
    }

    void RnaData::read_ps(const std::string &filename) {
	
	std::ifstream in(filename.c_str());
	std::string line;
		
	getline(in,line);
	if (line!="%!PS-Adobe-3.0 EPSF-3.0") {
	    throw wrong_format_failure();
	}
	
	bool contains_stacking=false;
	
	std::string seqname = "seq"; // default sequence name
	
	while (getline(in,line) && !has_prefix(line,"/sequence")) {
	    if (pimpl_->has_stacking_ && has_prefix(line,"% Probabilities for stacked pairs")) {
		contains_stacking=true;
	    } else if (has_prefix(line,"%delete next line to get rid of title")) {
		getline(in,line);
		std::istringstream in2(line);
		std::string s;
		while(in2 >> s) {
		    if (s.length()>=2 && s[0]=='(' && s[s.length()-1]==')') {
			seqname = s.substr(1,s.length()-2);
			break;
		    }
		}
		if (seqname.empty()) {
		    throw syntax_error_failure("improper title specification");
		}
	    }
	}
	
	if (!contains_stacking) {
	    pimpl_->has_stacking_=false;
	}

	if (!has_prefix(line,"/sequence")) {
	    throw syntax_error_failure("no sequence tag");
	}
	
	std::string seqstr="";
	{
	    std::string s;
	    while (in >> s && s!=")") {
		s = s.substr(0,s.size()-1); // chop of last character
		// cout << s <<endl;
		seqstr+=s;
	    }
	}
	
	//! sequence characters should be upper case, and 
	//! Ts translated to Us
	normalize_rna_sequence(seqstr);
        
	pimpl_->sequence_.append(MultipleAlignment::SeqEntry(seqname,seqstr));
	
	while (getline(in,line)) {
	    if (line.length()>4) {
		std::string type=line.substr(line.length()-4);
		if (type == "ubox"
		    ||
		    type == "lbox"
		    ) {
		
		    std::istringstream ss(line);
		    unsigned int i,j;
		    double p;
		    ss >> i >> j >> p;
		
		    p*=p;
		
		    //std::cout << i << " " << j << std::endl;
		
		    if (! (1<=i && i<j && j<=pimpl_->sequence_.length())) {
			std::cerr << "WARNING: Input dotplot "<<filename<<" contains invalid line " << line << " (indices out of range)" << std::endl;
		    } else {
			if (type=="ubox") {
			    pimpl_->arc_probs_(i,j)=p;
			}
			else if (pimpl_->has_stacking_ && type=="lbox") { // read a stacking probability
			    pimpl_->arc_2_probs_(i,j)=p; // we store the joint probability of (i,j) and (i+1,j-1)
			}
		    }
		}
	    }
	}
    } // end read_ps

    void RnaData::read_old_pp(const std::string &filename) {
	
	std::ifstream in(filename.c_str());
	
	std::string name;
	std::string seqstr;
	
	
	// ----------------------------------------
	// read sequence/alignment
    
	std::map<std::string,std::string> seq_map;
    
	std::string line;
	std::string sequence_anchor_string="";
	
	bool contains_stacking=false;
	
	while (get_nonempty_line(in,line) && line!="#") { // iterate through lines; stop at the first line that equals "#"
	    std::istringstream in(line);
	    in >> name >> seqstr;
	    
	    if ( in.fail() ) {
		throw wrong_format_failure();
	    }
	    
	    
	    if (name != "SCORE:") { // ignore the (usually first) line that begins with SCORE:
		if (name == "#C") {
		    sequence_anchor_string += seqstr;
		} else {
		    normalize_rna_sequence(seqstr);
		    seq_map[name] += seqstr;
		}
	    }
	}
	
	if ( line!="#" ) {
	    throw wrong_format_failure();
	}
	
	if (sequence_anchor_string!="") {
	    pimpl_->sequence_.set_annotation( MultipleAlignment::AnnoType::anchors, 
					      SequenceAnnotation( sequence_anchor_string ));
	}
	
	for (std::map<std::string,std::string>::iterator it=seq_map.begin(); it!=seq_map.end(); ++it) {
	    // std::cout << "SEQ: " << it->first << " " << it->second << std::endl;
	    pimpl_->sequence_.append(MultipleAlignment::SeqEntry(it->first,it->second));
	}
	
	
	// ----------------------------------------
	// read base pairs
    
	int i,j;
	double p;

	// std::cout << "LEN: " << len<<std::endl;
	
	while( get_nonempty_line(in,line) ) {
	    std::istringstream in(line);
	    
	    in>>i>>j>>p;
      
	    if ( in.fail() ) {
		throw syntax_error_failure("Invalid line \"" +line+ "\" does not specify base pair probability.");
	    }
	    
	    if (i>=j) {
		throw syntax_error_failure("Error in PP input line \""+line+"\" (i>=j).\n");
	    }
      
	    pimpl_->arc_probs_(i,j)=p;
	    
	    double p2;
	    
	    if (pimpl_->has_stacking_ && (in >> p2)) {
		pimpl_->arc_2_probs_(i,j)=p2; // p2 is joint prob of (i,j) and (i+1,j-1)
		contains_stacking=true;
	    }
	}
	
	if (!contains_stacking) {
	    pimpl_->has_stacking_=false;
	}
	
    } // end read_old_pp
    
    

    void RnaData::read_pp(const std::string &filename) {
	std::ifstream in(filename.c_str());
	read_pp(in);
    }
    
    std::istream &
    RnaData::read_pp(std::istream &in) {
	
	std::string line;
	
	// check header
	
	getline(in,line);
	if (!has_prefix(line,"#PP 2")) {
	    throw  wrong_format_failure();
	}

	pimpl_->read_pp_sequence(in);

	get_nonempty_line(in,line);
	if( line == "#SECTION BASEPAIRS" ) {
	    pimpl_->read_pp_arc_probabilities(in);
	} else {
	    throw syntax_error_failure("Expected base pair section header.");
	}
    
	return in;
    }

    std::istream &
    RnaDataImpl::read_pp_sequence(std::istream &in) {
	
	sequence_ = MultipleAlignment(in,MultipleAlignment::FormatType::PP);
	sequence_.normalize_rna_symbols();
	
	return in;
    }

    std::istream &
    ExtRnaData::read_pp(std::istream &in) {
	RnaData::read_pp(in);
	
	std::string line;
	get_nonempty_line(in,line);
	
	if( line == "#SECTION INLOOP" ) {
	    ext_pimpl_->read_pp_in_loop_probabilities(in);
	    ext_pimpl_->has_in_loop_probs_=true;
	} else {
	    ext_pimpl_->has_in_loop_probs_=false;
	}
	return in;
    }

    
    std::istream &
    RnaDataImpl::read_pp_arc_probabilities(std::istream &in) {
	// ----------------------------------------
	// read base pairs
    
	bool contains_stacking=false;

	// std::cout << "LEN: " << len<<std::endl;
	std::string line;
	while( get_nonempty_line(in,line) ) {
	    if (line[0]=='#') {
		// keyword line
		
		if (has_prefix(line,"#END")) {
		    // section end
		    break;
		} else if (has_prefix(line,"#BPCUT")) {
		    std::istringstream in(line);
		    std::string dummy;
		    double p;
		    in >> dummy >> p;
		    if ( in.fail() ) {
			throw syntax_error_failure("Cannot parse line \""+line+"\" in base pairs section.");
		    }
		    p_bpcut_ = std::max(p,p_bpcut_);
		}  else if (has_prefix(line,"#STACK")) {
		    contains_stacking = true;
		}
	    } else {
		std::istringstream in(line);
		size_t i,j;
		double p;
		in>>i>>j>>p;
      
		if ( in.fail() ) {
		    throw syntax_error_failure("Cannot parse line \""+line+"\" in base pairs section.");
		}
		
		if (!(1<=i && i<j && j<=sequence_.length())) {
		    throw syntax_error_failure("Invalid indices in PP input line \""+line+"\".");
		}
		
		if (p>p_bpcut_) {
		    arc_probs_(i,j)=p;
		
		    double p2;
		
		    if (has_stacking_ && (in >> p2)) {
			if (p2>p_bpcut_) {
			    arc_2_probs_(i,j)=p2; // p2 is joint prob of (i,j) and (i+1,j-1)
			}
		    }
		}
	    }
	}

	if (!contains_stacking && arc_2_probs_.size()>0) {
	    throw syntax_error_failure("Stacking probabilties found but stack keyword is missing.");
	}
	
	return in;
    }

    std::string
    read_pp_in_loop_block(const std::string &firstline,std::istream &in) {
	size_t pos = firstline.find(":");
	assert(pos != std::string::npos);
	
	std::string block = firstline.substr(pos+1);
	
	if (block.size()==0) { 
	    return block;
	}
	
	std::string line;
	while ( block[block.size()-1]=='\\' && getline(in,line) ) {
	    block = block.substr(0,block.size()-1);
	    block += line;
	}
	
	return block;
    }

    std::istream &
    ExtRnaDataImpl::read_pp_in_loop_probabilities(std::istream &in) {
    	std::string line;
	while( get_nonempty_line(in,line) ) {
	    if (has_prefix(line,"#END")) {
		// section ends
		break;
	    } else if (line[0]=='#') {
		read_pp_in_loop_probabilities_kwline(line);
	    } else {
		read_pp_in_loop_probabilities_line(line);
	    }
	    
	}
	return in;
    }

    void
    ExtRnaDataImpl::read_pp_in_loop_probabilities_kwline(const std::string &line) {
	if (has_prefix(line,"#BPILCUT")) {
	    std::istringstream in(line);
	    std::string dummy;
	    double p;
	    in >> dummy >> p;
	    if ( in.fail() ) {
		throw syntax_error_failure("Cannot parse line \""
					   +line+"\" in in-loop section.");
	    }
	    p_bpilcut_ = std::max(p,p_bpilcut_);
	} else if (has_prefix(line,"#UILCUT")) {
	    std::istringstream in(line);
	    std::string dummy;
	    double p;
	    in >> dummy >> p;
	    if ( in.fail() ) {
		throw syntax_error_failure("Cannot parse line \""
					   +line+"\" in in-loop section.");
	    }
	    p_uilcut_ = std::max(p,p_uilcut_);
	}
    }
    
    void
    ExtRnaDataImpl::read_pp_in_loop_probabilities_line(const std::string &line) {
	size_t i;
	size_t j;

	std::istringstream in(line);
	std::string sep;
	in>>i>>j>>sep;
		    
	if (sep!=":") {
	    throw syntax_error_failure("Invalid line \""
				       +line+"\" in in-loop section.");
	}
	if (!((1<=i && i<j && j<=self_->length()) // regular loop
	      || (i==0 && j==self_->length()+1) // external loop
	      )) {
	    throw syntax_error_failure("Index error in PP input line \""
				       +line+"\" (i>=j).");
	}
	
		
	std::string block_string = read_pp_in_loop_block(line,in);
		
	std::vector<std::string> blocks;
	split_at_separator(block_string,';',blocks);
		
	if (blocks.size() != 2) {
	    std::cerr << "Faulty block: "<<block_string<<std::endl;
	    throw syntax_error_failure("Invalid in loop probabilitity "
				       "specification at line \""
				       +line+"\"");
	}
		
	{
	    std::istringstream in1(blocks[0]);
	    size_t ip;
	    size_t jp;
	    double p;
	    while(in1 >> ip >> jp >> p) {
		if (!(i<ip && ip<jp && jp<j)) {
		    throw 
			syntax_error_failure("Index error in in-loop "
					     "specification.");
		}
		arc_in_loop_probs_.ref(i,j).set(ip,jp,p);
	    }
	}
		
	{
	    std::istringstream in2(blocks[1]);
	    size_t kp;
	    double p;
	    while(in2 >> kp >> p) {
		if (!(i<kp && kp<j)) {
		    throw 
			syntax_error_failure("Index error in in-loop "
					     "specification.");
		}
		unpaired_in_loop_probs_.ref(i,j)[kp] = p;
	    }
	}
    }
    
    std::ostream &
    RnaData::write_pp(std::ostream &out,
		      double p_outbpcut) const {
	
	out << "#PP 2.0"
	    << std::endl  
	    << std::endl;
	
	pimpl_->write_pp_sequence(out);
	
	pimpl_->write_pp_arc_probabilities(out,p_outbpcut,pimpl_->has_stacking_);
	
	return out;
    }
    
    std::ostream &
    ExtRnaData::write_pp(std::ostream &out,
    			 double p_outbpcut,
			 double p_outbpilcut,
    			 double p_outuilcut
    			 ) const {
	
	RnaData::write_pp(out,p_outbpcut);
	
	ext_pimpl_->write_pp_in_loop_probabilities(out,
					       p_outbpcut,
					       p_outbpilcut,
					       p_outuilcut);
	
    	return out;
    }
    

    std::ostream &
    RnaDataImpl::write_pp_sequence(std::ostream &out) const {
	out << sequence_;
	
	out << std::endl
	    << "#END" << std::endl;
	
	return out;
    }


    /** 
     * @brief output format for probabilities in pp files
     * use limited precision; use scientific notation if it is shorter
     */
    std::string
    format_prob(double prob) {
	std::ostringstream outd;
	outd.precision(3);
	outd << prob;
	
	if (outd.str().length()<=6) {
	    return outd.str();
	}

	std::ostringstream outs;
	outs.setf( std::ios::scientific, std:: ios::floatfield );
	outs.precision(2);
	outs << prob;
	
	std::string s=outs.str();
	size_t pos = s.find("e-0");
	if (pos!=std::string::npos) {
	    s.replace(pos,3,"e-");
	}
	
	return s;
    }
    
    /**
     * @brief Write arc probabilities
     * 
     * Writes arc and stacking probabilities to stream; filters by
     * probability threshold p_outbpcut
     */
    std::ostream &
    RnaDataImpl::write_pp_arc_probabilities(std::ostream &out,
					    double p_outbpcut,
					    bool stacking
					    ) const {
	
	out << std::endl
	    << "#SECTION BASEPAIRS" << std::endl
	    << std::endl
	    << "#BPCUT "<<format_prob(std::max(p_bpcut_,p_outbpcut)) << std::endl;
	
	if (stacking) {
	    out << "#STACK"<<std::endl;
	}
	out <<std::endl;
	
	// assume that for each entry in arc_2_probs_ there is a corresponding entry in arc_probs_
#     ifndef NDEBUG
	for (arc_prob_matrix_t::const_iterator it = arc_2_probs_.begin();
	     arc_2_probs_.end() != it;
	     ++it) {
	    assert(arc_probs_(it->first.first,it->first.second)!=0.0);
	}
#     endif
	
	for (arc_prob_matrix_t::const_iterator it = arc_probs_.begin();
	     arc_probs_.end() != it;
	     ++it) {
	    size_t i=it->first.first;
	    size_t j=it->first.second;
	    if (it->second > p_outbpcut) {
		out << i << " " << j << " " << format_prob(it->second);
		if (stacking && has_stacking_ && arc_2_probs_(i,j)>p_bpcut_) {
		    out << " " << format_prob(arc_2_probs_(i,j));
		}
		out << std::endl;
	    }
	}

	out << std::endl
	    << "#END" << std::endl;
	
	return out;
    }

    std::ostream &
    ExtRnaDataImpl::write_pp_in_loop_probabilities(std::ostream &out,
						   double p_outbpcut,
						   double p_outbpilcut,
						   double p_outuilcut
						   ) const {
	out << std::endl
	    << "#SECTION INLOOP" << std::endl
	    << std::endl
	    << "#BPILCUT " << format_prob(std::max(p_bpilcut_,p_outbpilcut)) << std::endl
	    << "#UILCUT  " << format_prob(std::max(p_uilcut_,p_outuilcut)) << std::endl
	    << std::endl;
	
	// write in-loop probabilities for all arcs with probability greater than p_outbpcut
	for (arc_prob_matrix_t::const_iterator it = self_->arc_probs_begin();
	     self_->arc_probs_end() != it;
	     ++it) {
	    if (it->second > p_outbpcut) {
		write_pp_in_loop_probability_line(out,
						  it->first.first,
						  it->first.second,
						  p_outbpilcut,
						  p_outuilcut);
	    }
	}

	// write in loop probs for external loop
	write_pp_in_loop_probability_line(out,
					  0,self_->length()+1,
					  p_outbpilcut,
					  p_outuilcut);
	
	out << std::endl
	    << "#END" << std::endl;
	
	return out;
    }
    
    std::ostream &
    ExtRnaDataImpl::write_pp_in_loop_probability_line(std::ostream &out,
						      size_t i, size_t j,
						      double p_bpilcut,
						      double p_uilcut) const {
	out << i << " " << j << " :";
	// if (arc_in_loop_probs_(i,j).size()>=5) {
	//     out << std::endl << "   ";
	// }
	
	write_pp_basepair_in_loop_probabilities(out, arc_in_loop_probs_(i,j),
						p_bpilcut);
	
	out << " ;"; // separate base pair and unpaired probabilities
	if (arc_in_loop_probs_(i,j).size()>=4 && unpaired_in_loop_probs_(i,j).size()>=4) {
	    out << "\\" << std::endl << "   ";
	}
	
	write_pp_unpaired_in_loop_probabilities(out, unpaired_in_loop_probs_(i,j),
						p_uilcut);
	out << std::endl;
	
	return out;
    }
    
    std::ostream &
    ExtRnaDataImpl::write_pp_basepair_in_loop_probabilities(std::ostream &out,
							    const arc_prob_matrix_t &probs,
							    double p_cut) const {
	for (arc_prob_matrix_t::const_iterator it=probs.begin(); probs.end()!=it; ++it) {
	    if (it->second > p_cut) {
		out << " " << it->first.first
		    << " " << it->first.second
		    << " " << format_prob(it->second);
	    }
	}
	return out;
    }

    std::ostream &
    ExtRnaDataImpl::write_pp_unpaired_in_loop_probabilities(std::ostream &out,
							    const arc_prob_vector_t &probs,
							    double p_cut) const {
	for (arc_prob_vector_t::const_iterator it=probs.begin(); probs.end()!=it; ++it) {
	    if (it->second > p_cut) {
		out << " " << it->first << " " << format_prob(it->second);
	    }
	}
	return out;
    }

    
    std::ostream &
    RnaData::write_size_info(std::ostream &out) const {
	out << "arcs: "<<pimpl_->arc_probs_.size();
	if (pimpl_->has_stacking_) {
	    out << "  stackings: "<<pimpl_->arc_2_probs_.size();
	}
	return out;
    }

    std::ostream &
    ExtRnaData::write_size_info(std::ostream &out) const {
	// count arcs in loop
	size_t num_arcs_in_loop=0;
	// count unpaired bases in loop
	size_t num_unpaired_in_loop=0;

	size_t len = length();
	for( size_t i=1; i <= len; i++ ) {
	    for( size_t j=i+1; j <= len; j++ ) {
		num_arcs_in_loop += 
		    arc_prob_matrix_t(ext_pimpl_->arc_in_loop_probs_(i,j)).size();
		num_unpaired_in_loop += 
		    ExtRnaDataImpl::arc_prob_vector_t(ext_pimpl_->unpaired_in_loop_probs_(i,j)).size();	
	    }
	}
	
	return
	    RnaData::write_size_info(out)
	    <<"  arcs in loops: " << num_arcs_in_loop << "  unpaireds in loops: " << num_unpaired_in_loop;
    }


    void
    RnaDataImpl::init_as_consensus_dot_plot(const Alignment::edges_t &edges,
					    const RnaData &rna_dataA,
					    const RnaData &rna_dataB,
					    double p_expA,
					    double p_expB,
					    bool stacking
					    ) {
	
	size_t rowsA = rna_dataA.sequence().num_of_rows();
	size_t rowsB = rna_dataB.sequence().num_of_rows();
	
	double p_minA = rna_dataA.arc_cutoff_prob();
	double p_minB = rna_dataB.arc_cutoff_prob();
	
	double p_minMean =
	    exp(
		(log(p_minA)*rowsA + log(p_minB)*rowsB)
		/ (rowsA+rowsB)
		);
	
	p_bpcut_ = p_minMean;
	
	for (size_type i=0; i<edges.size(); i++) {
	    for (size_type j=i+1; j<edges.size(); j++) {
		// here we compute consensus pair probabilities
		
		double pA =
		    (edges.first[i].is_gap() || edges.first[j].is_gap())
		    ? 0
		    : rna_dataA.arc_prob(edges.first[i], edges.first[j]);

		double pB =
		    (edges.second[i].is_gap() || edges.second[j].is_gap())
		    ? 0
		    : rna_dataB.arc_prob(edges.second[i], edges.second[j]);

		double p = consensus_probability(pA,pB,rowsA,rowsB,p_expA,p_expB);

		if (stacking) {
		    double st_pA =
			(edges.first[i].is_gap() || edges.first[j].is_gap())
			? 0
			: rna_dataA.joint_arc_prob(edges.first[i], edges.first[j]);
		    
		    double st_pB =
			(edges.second[i].is_gap() || edges.second[j].is_gap())
			? 0
			: rna_dataB.joint_arc_prob(edges.second[i], edges.second[j]);

		    double st_p = consensus_probability(st_pA,st_pB,rowsA,rowsB,p_expA,p_expB);

		    if (p > p_minMean || st_p > p_minMean) {
			arc_probs_(i+1,j+1) = p;
			arc_2_probs_(i+1,j+1) = st_p;
		    }

		} else {
		    if (p > p_minMean) {
			arc_probs_(i+1,j+1) = p;
		    }
		}
	    }
	}
    }
    
    double
    RnaDataImpl::consensus_probability(double pA, double pB,
				       size_t sizeA,size_t sizeB,
				       double p_expA, double p_expB) const {
	pA = std::max(std::min(p_expA,p_bpcut_*0.75), pA);
	pB = std::max(std::min(p_expB,p_bpcut_*0.75), pB);

	// weighted geometric mean
	double p = exp(
		       (log(pA)*sizeA + log(pB)*sizeB) /
		       //---------------------------------------------------
		       (sizeA + sizeB)
		       );
	
	return p;

	/*
	  would something like
	  if (pA<p_min*1.05) { pA = std::min(p_expA,p_min*0.75); }
	  work better???
	*/
    }
    
    void
    RnaDataImpl::drop_worst_bps(size_t keep) {
	typedef keyvec<arc_prob_matrix_t::key_t> kv_t;

	kv_t::vec_t vec;
	
	//std::copy(arc_probs_.begin(),arc_probs_.end(),vec.begin());
	for (arc_prob_matrix_t::const_iterator it=arc_probs_.begin();
	     arc_probs_.end() != it;
	     ++it ) {
	    vec.push_back(*it);
	}

	std::make_heap(vec.begin(),vec.end(),kv_t::comp);
	
	while(vec.size()>keep) {
	    const arc_prob_matrix_t::key_t &key = vec.front().first;
	    
	    arc_probs_(key.first,key.second) = 0.0;
	    arc_2_probs_(key.first,key.second) = 0.0;
	    std::pop_heap(vec.begin(),vec.end(),kv_t::comp);
	    vec.pop_back();
	}
    }

    void
    ExtRnaDataImpl::drop_worst_bps(size_t keep) {

        // access pimpl_ of parent RnaData object
	RnaDataImpl *rdimpl = static_cast<RnaData *>(self_)->pimpl_;
	rdimpl->drop_worst_bps(keep);
	
	// free unpaired in loop where arc prob is 0
	for (arc_prob_vector_matrix_t::const_iterator it = unpaired_in_loop_probs_.begin();
	     unpaired_in_loop_probs_.end() != it;
	     ++it) {
	    arc_prob_vector_matrix_t::key_t key = it->first;
	    if ( rdimpl->arc_probs_(key.first,key.second) == 0.0 ) {
		if (key.first==0) continue;
		unpaired_in_loop_probs_.reset(key.first,key.second);
	    }
	}
	
	// free base pairs in loop where arc prob is 0
	for (arc_prob_matrix_matrix_t::const_iterator it = arc_in_loop_probs_.begin();
	     arc_in_loop_probs_.end() != it;
	     ++it) {
	    arc_prob_matrix_matrix_t::key_t key = it->first; 
	    if ( rdimpl->arc_probs_(key.first,key.second) == 0.0 ) {
		if (key.first==0) continue;
		arc_in_loop_probs_.reset(key.first,key.second);
	    } else {
		for (arc_prob_matrix_t::const_iterator it2 = it->second.begin();
		     it->second.end() != it2;
		     ++it2) {
		    arc_prob_matrix_matrix_t::key_t key2 = it2->first; 
		    if ( rdimpl->arc_probs_(key2.first,key2.second) == 0.0 ) {
			arc_in_loop_probs_.ref(key.first,key.second)
			    .reset(key2.first,key2.second);
		    }	
		}
	    }
	}

    }

    void
    ExtRnaDataImpl::drop_worst_uil(size_t keep) {
	
	typedef std::pair< arc_prob_vector_matrix_t::key_t, arc_prob_vector_t::key_t > key_t;

	typedef RnaDataImpl::keyvec<key_t> kv_t;

	kv_t::vec_t vec;
	
	// push all uil probs with their key to vector vec
	for (arc_prob_vector_matrix_t::const_iterator it=unpaired_in_loop_probs_.begin();
	     unpaired_in_loop_probs_.end() != it;
	     ++it ) {
	    
	    for (arc_prob_vector_t::const_iterator it2=it->second.begin();
		 it->second.end() != it2;
		 ++it2 ) {
		
		vec.push_back(kv_t::kvpair_t(key_t(it->first,it2->first), it2->second));
	    }
	}

	std::make_heap(vec.begin(),vec.end(),kv_t::comp);
	
	while(vec.size()>keep) {
	    const key_t &key = vec.front().first;
	    
	    unpaired_in_loop_probs_.ref(key.first.first,key.first.second).reset(key.second);
	    
	    std::pop_heap(vec.begin(),vec.end(),kv_t::comp);
	    vec.pop_back();
	}

    }

    void
    ExtRnaDataImpl::drop_worst_bpil(size_t keep) {
	
	typedef std::pair< arc_prob_matrix_matrix_t::key_t, arc_prob_matrix_t::key_t > key_t;
	
	typedef RnaDataImpl::keyvec<key_t> kv_t;

	kv_t::vec_t vec;
	
	// push all uil probs with their key to vector vec
	for (arc_prob_matrix_matrix_t::const_iterator it=arc_in_loop_probs_.begin();
	     arc_in_loop_probs_.end() != it;
	     ++it ) {
	    
	    for (arc_prob_matrix_t::const_iterator it2=it->second.begin();
		 it->second.end() != it2;
		 ++it2 ) {
		
		vec.push_back(kv_t::kvpair_t(key_t(it->first,it2->first), it2->second));
	    }
	}

	std::make_heap(vec.begin(),vec.end(),kv_t::comp);
	
	while(vec.size()>keep) {
	    const key_t &key = vec.front().first;
	    
	    arc_in_loop_probs_.ref(key.first.first,key.first.second).reset(key.second.first,key.second.second);
	    
	    std::pop_heap(vec.begin(),vec.end(),kv_t::comp);
	    vec.pop_back();
	}

    }


    void
    ExtRnaDataImpl::drop_worst_bpil_precise(double ratio) {

	typedef std::pair< arc_prob_matrix_matrix_t::key_t, arc_prob_matrix_t::key_t > key_t;

	typedef RnaDataImpl::keyvec<key_t> kv_t;


	// push all uil probs with their key to vector vec
	for (arc_prob_matrix_matrix_t::const_iterator it=arc_in_loop_probs_.begin();
	     arc_in_loop_probs_.end() != it;
	     ++it ) {
		kv_t::vec_t vec;
		vec.clear();
	    for (arc_prob_matrix_t::const_iterator it2=it->second.begin();
		 it->second.end() != it2;
		 ++it2 ) {

		vec.push_back(kv_t::kvpair_t(key_t(it->first,it2->first), it2->second));
	    }
		double keep =  ratio * ((double)(it->first.second)- (double)(it->first.first) + 1) ;
		if (vec.size()> keep)
		{
			std::make_heap(vec.begin(),vec.end(),kv_t::comp);
			while(vec.size()> keep ) {
				const key_t &key = vec.front().first;
				arc_in_loop_probs_.ref(key.first.first,key.first.second).reset(key.second.first,key.second.second);

				std::pop_heap(vec.begin(),vec.end(),kv_t::comp);
				vec.pop_back();
			}
		}

	}


    }

    vrna_plist_t *
    RnaData::plist() const {
        std::vector<vrna_plist_t> plist;
        size_type len = length();
        for(size_t i=1; i<=len; ++i) {
            for(size_t j=1; j<=len; ++j) {
                double p=arc_prob(i,j);
                if (p>0) {
                    vrna_plist_t x;
                    x.i=i;
                    x.j=j;
                    x.p=p;
                    x.type=0;
                    plist.push_back(x);
                }
            }
        }
        
        // construct Vienna RNA / C - compatible array
        vrna_plist_t * c_plist = new vrna_plist_t[plist.size()+1];
        // and copy contents of the vecctor to this array
        copy(plist.begin(), plist.end(), c_plist);
        // mark end of list
        plist[plist.size()].i=0;
        plist[plist.size()].j=0;
        plist[plist.size()].p=0;
        plist[plist.size()].type=0;
        
        return c_plist;
    }

    std::string
    RnaData::mea_structure(double gamma) const {
        vrna_plist_t *pl = plist();
        
        char *c_structure = new char[length()+1];
        for (size_t i=0; i<length(); ++i) {c_structure[i]='.';}
        c_structure[length()]=0;
        
        // call RNAlib's mea function
        float mea = MEA(pl,c_structure,gamma);
        
        std::string structure(c_structure);
        
        delete [] c_structure;
        delete [] pl;
        
        return structure;
    }

} // end namespace LocARNA
