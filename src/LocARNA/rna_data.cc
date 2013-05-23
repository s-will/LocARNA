#include <string>
#include <fstream>
#include <sstream>
#include "pfold_params.hh"
#include "rna_data_impl.hh"
#include "rna_ensemble.hh"



namespace LocARNA {

    RnaData::RnaData(const RnaEnsemble &rna_ensemble,
		     double p_bpcut)
	: pimpl_(new RnaDataImpl(this,
				 rna_ensemble,
				 p_bpcut)) {
    }

    RnaData::RnaData(const std::string &filename,
		     double p_bpcut,
		     const PFoldParams &pfoldparams)
	: pimpl_(new RnaDataImpl(this,
				 filename,
				 p_bpcut,
				 pfoldparams)) {
    }

    
    // do almost nothing
    RnaData::RnaData(double p_bpcut)
	: pimpl_(new RnaDataImpl(this,
				 p_bpcut)) {
    }


    RnaData::~RnaData() {
	delete pimpl_;
    }


    RnaDataImpl::RnaDataImpl(RnaData *self,
			     const std::string &filename,
			     double p_bpcut,
			     const PFoldParams &pfoldparams)
	:self_(self),
	 sequence_(),
	 p_bpcut_(p_bpcut),
	 arc_probs_(0.0),
	 arc_2_probs_(0.0),
	 sequence_anchors_(),
	 stacking_probs_available_(false)
    {	
	// HACK: read_autodetect needs an initialized pimpl_ 
	// how could we do this more cleanly?
	self_->pimpl_=this;
	
	self_->read_autodetect(filename,pfoldparams);
    }

    RnaDataImpl::RnaDataImpl(RnaData *self,
			     const RnaEnsemble &rna_ensemble,
			     double p_bpcut)
	:self_(self),
	 sequence_(),
	 p_bpcut_(p_bpcut),
	 arc_probs_(0.0),
	 arc_2_probs_(0.0),
	 sequence_anchors_(),
	 stacking_probs_available_(false)
    {
	// HACK: we need an initialized pimpl_ 
	self_->pimpl_=this;

	self_->init_from_rna_ensemble(rna_ensemble);
    }

    // do almost nothing
    RnaDataImpl::RnaDataImpl(RnaData *self,
			     double p_bpcut)
	:self_(self),
	 sequence_(),
	 p_bpcut_(p_bpcut),
	 arc_probs_(0.0),
	 arc_2_probs_(0.0),
	 sequence_anchors_(),
	 stacking_probs_available_(false)
    {
    }
    
    ExtRnaData::ExtRnaData(const RnaEnsemble &rna_ensemble,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut)
	:
	RnaData(rna_ensemble,
		p_bpcut),				  
	pimpl_(new ExtRnaDataImpl(this,
				  rna_ensemble,
				  p_bpilcut,
				  p_uilcut)) {
    }
    
    ExtRnaData::ExtRnaData(const std::string &filename,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut,
			   const PFoldParams &pfoldparams)
	: 
	RnaData(p_bpcut),
	pimpl_(new ExtRnaDataImpl(this,
				  filename,
				  p_bpilcut,
				  p_uilcut,
				  pfoldparams)) {
    }

    ExtRnaData::~ExtRnaData() {
	delete pimpl_;
    }


    ExtRnaDataImpl::ExtRnaDataImpl(ExtRnaData *self,
				   const std::string &filename,
				   double p_bpilcut,
				   double p_uilcut,
				   const PFoldParams &pfoldparams)
	:self_(self),
	 p_bpilcut_(p_bpilcut),
	 p_uilcut_(p_uilcut),
	 arc_in_loop_probs_(arc_prob_matrix_t(0.0)),
	 unpaired_in_loop_probs_(arc_prob_vector_t(0.0)),
	 in_loop_probs_available_(false)
    {
	// HACK: read_autodetect needs an initialized pimpl_ 
	// how could we do this more cleanly?
	self_->pimpl_=this;
	
	self_->read_autodetect(filename,pfoldparams);
    }

    ExtRnaDataImpl::ExtRnaDataImpl(ExtRnaData *self,
				   const RnaEnsemble &rna_ensemble,
				   double p_bpilcut,
				   double p_uilcut)
	:self_(self),
	 p_bpilcut_(p_bpilcut),
	 p_uilcut_(p_uilcut),
	 arc_in_loop_probs_(arc_prob_matrix_t(0.0)),
	 unpaired_in_loop_probs_(arc_prob_vector_t(0.0)),
	 in_loop_probs_available_(false)
    {
	// HACK: we need an initialized pimpl_ 
	self_->pimpl_=this;

	self_->init_from_rna_ensemble(rna_ensemble);
    }

    /*
     * Implementation note: we use an annoyingly complex mechanism
     * involving virtual methods and calls to interface-object methods
     * from the implementation object's constructor to make this
     * method usable in RnaData and ExtRnaData while avoiding code
     * duplication. In hindsight, this seems rather pathetic. Could it
     * be replaced by a simpler mechanism?
     */
    void
    RnaData::read_autodetect(const std::string &filename,
			     const PFoldParams &pfoldparams
			     ) {
	bool failed=true;  //flag for signalling a failed attempt to
			   //read a certain file format
	
	bool recompute=false; // do we need to recompute probabilities

	if (failed) {
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
	
	if (failed) {
	    recompute = true;
	    failed=false;
	    try {
		MultipleAlignment ma(filename, MultipleAlignment::FASTA);
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
	
	if (failed) {
	    recompute = true;
	    failed=false;
	    try {
		MultipleAlignment ma(filename, MultipleAlignment::CLUSTAL);

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

	// if MultipleAlignment would support AUTO, we could replace the above by
	// if (failed) {
	//     pimpl_->sequence_ = MultipleAlignment(filename,MultipleAlignment::AUTO);
	//     recompute = true;
	// }

	
	if (failed) {
	    recompute=false;
	    failed=false;
	    try {
		read_pp(filename, pfoldparams.stacking());
		if (!pimpl_->sequence_.is_proper() || pimpl_->sequence_.empty() ) {
		    failed=true;
		}
	    } catch (wrong_format_failure &f) {
		failed=true;
	    }
	}
	
	if (failed) {
	    throw failure("Cannot read input data from file.");
	}
	
	// now, we have the sequence but not necessarily all required probabilities!
	if (!inloopprobs_ok()) {
	    recompute=true;
	}
	
	if (recompute) {
	    
	    // recompute all probabilities
	    RnaEnsemble rna_ensemble(pimpl_->sequence_,pfoldparams,!inloopprobs_ok(),true); // use given parameters, use alifold
	    
	    // initialize from (temporary) RnaEnsemble object; note that the method is virtual
	    init_from_rna_ensemble(rna_ensemble);
	}
	
	return;
    }

    void
    RnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble) {
	pimpl_->init_from_rna_ensemble(rna_ensemble);
    }

    void
    ExtRnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble) {
	RnaData::init_from_rna_ensemble(rna_ensemble);
	pimpl_->init_from_ext_rna_ensemble(rna_ensemble);
    }
    
    void
    RnaDataImpl::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble) {
	assert(rna_ensemble.has_base_pair_probs());
	
	// ----------------------------------------
	// init sequence	
	sequence_ = rna_ensemble.sequence();
	size_t len = sequence_.length();

	// ----------------------------------------
	// init base pair probabilities
	arc_probs_.clear();
	for( size_t i=1; i <= len; i++ ) {
	    for( size_t j=i+TURN+1; j <= len; j++ ) {
		
		double p = rna_ensemble.arc_prob(i,j);
		if (p > p_bpcut_) { // apply filter
		    arc_probs_(i,j)=p;
		}
	    }
	}
	
	// ----------------------------------------
	// init stacking probabilities
	arc_2_probs_.clear();
	stacking_probs_available_=rna_ensemble.has_stacking_probs();	
	if (stacking_probs_available_) {
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
	
	// in loop
	for(arc_prob_matrix_t::const_iterator it = self_->arc_probs_begin();
	    self_->arc_probs_end()!=it; ++it) {
	    pos_type i = it->first.first;
	    pos_type j = it->first.second;
	    arc_prob_matrix_t m_ij(0.0);
	    
	    for( size_t ip=i+1; ip < j; ip++ ) {
		for( size_t jp=ip+TURN+1; jp < j; jp++ ) {
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
	    for( size_t jp=ip+TURN+1; jp <= len; jp++ ) {
		double p = rna_ensemble.arc_external_prob(ip,jp);
		
		if ( p > p_bpilcut_ ) {
		    m_ext(ip,jp) = p;
		}
	    }
	}

	// set only if not empty; use set instead of assignment,
	// to avoid the comparison of complex SparseMatrix objects
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
	    
	    // set only if not empty; use set instead of assignment,
	    // to avoid the comparison of complex SparseMatrix objects
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

	// set only if not empty; use set instead of assignment,
	// to avoid the comparison of complex SparseMatrix objects
	if (!v_ext.empty()) {
	    unpaired_in_loop_probs_.set(0,self_->length()+1,v_ext);
	}
	
	
	// set flag
	in_loop_probs_available_=true;

	// all set
	return;
    }

    bool
    ExtRnaData::inloopprobs_ok() const {
	return pimpl_->in_loop_probs_available_;
    }

    const Sequence &
    RnaData::sequence() const {
	return pimpl_->sequence_;
    }

    size_type
    RnaData::length() const {
	return pimpl_->sequence_.length();
    }
    
    const std::string &
    RnaData::sequence_anchors() const {
	return pimpl_->sequence_anchors_;
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
	return pimpl_->p_bpilcut_;
    }
	
    double 
    ExtRnaData::arc_in_loop_prob(pos_type i, pos_type j,pos_type p, pos_type q) const {
	ExtRnaDataImpl::arc_prob_matrix_t m_pq = pimpl_->arc_in_loop_probs_(p,q);
	return m_pq(i,j);
    }
    
    double 
    ExtRnaData::arc_external_prob(pos_type i, pos_type j) const {
	ExtRnaDataImpl::arc_prob_matrix_t m_ext = pimpl_->arc_in_loop_probs_(0,length()+1);
	return m_ext(i,j);
    }
    
    double
    ExtRnaData::unpaired_in_loop_cutoff_prob() const {
	return pimpl_->p_uilcut_;
    }
    
    double 
    ExtRnaData::unpaired_in_loop_prob(pos_type k,pos_type p, pos_type q) const {
	ExtRnaDataImpl::arc_prob_vector_t v_pq = pimpl_->unpaired_in_loop_probs_(p,q);
	return v_pq[k];
    }
    
    double 
    ExtRnaData::unpaired_external_prob(pos_type k) const {
	ExtRnaDataImpl::arc_prob_vector_t v_ext = pimpl_->unpaired_in_loop_probs_(0,length()+1);
	return v_ext[k];
    }

    void RnaData::read_ps(const std::string &filename) {
	
	std::ifstream in(filename.c_str());
	std::string line;
		
	getline(in,line);
	if (line!="%!PS-Adobe-3.0 EPSF-3.0") {
	    throw wrong_format_failure();
	}

	std::string seqname = "seq";
	
	pimpl_->stacking_probs_available_=false;
	
	while (getline(in,line) && has_prefix(line,"/sequence")) {
	    if (has_prefix(line,"% Probabilities for stacked pairs")) {
		pimpl_->stacking_probs_available_=true;
	    } else if (has_prefix(line,"%delete next line to get rid of title")) {
		getline(in,line);
		std::istringstream in2(line);
		std::string s;
		while(in2 >> s && s.length()>=2 && s[0]!='(') {
		    seqname = s.substr(1,s.length()-2);
		}
	    }
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
        
	pimpl_->sequence_.append(Sequence::SeqEntry(seqname,seqstr));
	
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
			else if (pimpl_->stacking_probs_available_ && type=="lbox") { // read a stacking probability
			    pimpl_->arc_2_probs_(i,j)=p; // we store the joint probability of (i,j) and (i+1,j-1)
			}
		    }
		}
	    }
	}
    } // end read_ps


    void RnaData::read_pp(const std::string &filename,
			  double p_incut
			  ) {
	
	std::ifstream in(filename.c_str());
	
	std::string name;
	std::string seqstr;
	
	
	// ----------------------------------------
	// read sequence/alignment
    
	std::map<std::string,std::string> seq_map;
    
	std::string line;
    
	while (getline(in,line) && (line.length()==0 || line[0]!='#') ) { // iterate through lines; stop at the first line that starts with '#'
	    if (line.length()>0 && line[0]!=' ') { // ignore empty lines and lines that start with blank
		std::istringstream in(line);
		in >> name >> seqstr;
	    
		if ( in.fail() ) {
		    throw wrong_format_failure();
		}


		if (name != "SCORE:") { // ignore the (usually first) line that begins with SCORE:
		    if (name == "#C") {
			pimpl_->sequence_anchors_ += seqstr;
		    } else {
			normalize_rna_sequence(seqstr);
			seq_map[name] += seqstr;
		    }
		}
	    }
	}
	
	// check whether '#'-line contains further annotation
	{
	    char c;
	    double p;
	    std::string s;
	    
	    std::istringstream in(line);
	    
	    in >> c; // eat '#'
	    if ( c!='#' ) {
		throw wrong_format_failure();
	    }
	    
	    in >> p;	    
	    if (!in.fail()) {
		pimpl_->p_bpcut_ = p;
	    }

	    in >> s;	    
	    if (has_prefix(s,"stack")) {
		pimpl_->stacking_probs_available_=true;
	    }
	    
	}

	for (std::map<std::string,std::string>::iterator it=seq_map.begin(); it!=seq_map.end(); ++it) {
	    // std::cout << "SEQ: " << it->first << " " << it->second << std::endl;
	    pimpl_->sequence_.append(Sequence::SeqEntry(it->first,it->second));
	}
    
	// ----------------------------------------
	// read base pairs
    
	int i,j;
	double p;

	// std::cout << "LEN: " << len<<std::endl;
	
	while( getline(in,line) ) {
	    std::istringstream in(line);
      
	    in>>i>>j>>p;
      
	    if ( in.fail() ) continue; // skip lines that do not specify a base pair probability
      
	    if (i>=j) {
		std::ostringstream err;
		err << "Error in PP input line \""<<line<<"\" (i>=j).\n"<<std::endl;
		throw failure(err.str());
	    }
      
	    pimpl_->arc_probs_(i,j)=p;
      
	    double p2;
	    
	    if (in >> p2) {
		pimpl_->arc_2_probs_(i,j)=p2; // p2 is joint prob of (i,j) and (i+1,j-1)
		pimpl_->stacking_probs_available_ = true;
	    }
	}
    } // end read_pp

    std::ostream &
    RnaData::write_size_info(std::ostream &out) const {
	out << "arcs: "<<pimpl_->arc_probs_.size();
	if (pimpl_->stacking_probs_available_) {
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
		    arc_prob_matrix_t(pimpl_->arc_in_loop_probs_(i,j)).size();
		num_unpaired_in_loop += 
		    ExtRnaDataImpl::arc_prob_vector_t(pimpl_->unpaired_in_loop_probs_(i,j)).size();	
	    }
	}
	
	return
	    RnaData::write_size_info(out)
	    <<"  arcs in loops: " << num_arcs_in_loop << "  unpaireds in loops: " << num_unpaired_in_loop;
    }


    // -- writing stuff (old code from RnaEnsemble)
    // std::ostream &
    // RnaEnsemble::write_unpaired_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
	
    // 	// write lines for loops closed by base pairs
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type i=1; i<=get_length(); ++i) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type j=i+1; j<=get_length(); ++j) {
    // 		if (pimpl_->arc_probs_(i,j)>threshold1) {
    // 		    bool had_entries=false;
    // 		    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type k=i+1; k<=j-1; ++k) {
    // 			double p=prob_unpaired_in_loop(k,i,j);
    // 			if (p>threshold2) {
    // 			    if (!had_entries) {
    // 				out << i << " " << j; had_entries=true;
    // 			    }
    // 			    out << " " << k << " " << p;
    // 			}
    // 		    }
    // 		    if (had_entries) out << std::endl;
    // 		}
    // 	    }
    // 	}
	
    // 	// write lines for external loop
	
    // 	bool had_entries=false;
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type k=1; k<=get_length(); ++k) {
    // 	    double p=prob_unpaired_external(k);
    // 	    if (p>threshold2) {
    // 		if (!had_entries) {
    // 		    out << 0 << " " << (get_length()+1); had_entries=true;
    // 		}
    // 		out << " " << k << " " << p;
    // 	    }
    // 	}
    // 	if (had_entries) out << std::endl;
    // 	return out;
    // }
    
	
    // std::ostream &
    // RnaEnsemble::write_basepair_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
    // 	// write lines for loops closed by base pairs
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type i=1; i<=get_length(); ++i) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type j=i+1; j<=get_length(); ++j) {
    // 		if (pimpl_->arc_probs_(i,j)>threshold1) {
    // 		    bool had_entries=false;
    // 		    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type ip=i+1; ip<=j-1; ++ip) {
    // 			for(RnaEnsembleImpl::arc_prob_matrix_t::size_type jp=ip+1; jp<=j-1; ++jp) {
    // 			    double p=prob_basepair_in_loop(ip,jp,i,j);
    // 			    if (p>threshold2) {
    // 				if (!had_entries) {out << i << " " << j; had_entries=true;}
    // 				out << " " << ip << " " << jp << " " << p;
    // 			    }
    // 			}
    // 		    }
    // 		    if (had_entries) out << std::endl;
    // 		}
    // 	    }
    // 	}
	
    // 	// write lines for external loop
    // 	bool had_entries=false;
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type ip=1; ip<=get_length(); ++ip) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type jp=ip+1; jp<=get_length(); ++jp) {
    // 		double p=prob_basepair_external(ip,jp);
    // 		if (p>threshold2) {
    // 		    if (!had_entries) {out << 0 << " " << (get_length()+1); had_entries=true;}
    // 		    out << " " << ip << " " << jp << " " << p;
    // 		}
    // 	    }
    // 	}
    // 	if (had_entries) out << std::endl;
    // 	return out;
    // }
    
    
    // std::ostream &
    // RnaEnsemble::write_basepair_and_in_loop_probs(std::ostream &out,double threshold1,double threshold2,double threshold3, bool write_probs, bool diff_encoding) const {
	
    // 	size_t i=0;
    // 	size_t j=get_length()+1;

    // 	RnaEnsembleImpl::arc_prob_matrix_t::size_type last_i=i;
    // 	RnaEnsembleImpl::arc_prob_matrix_t::size_type last_j=j;
	

    // 	// write line for external loop
	
    // 	out << (diff_encoding?(int)i-(int)last_i:(int)i) << " " << (diff_encoding?(int)last_j-(int)j:(int)j)<< " 1 ;";
	
    // 	RnaEnsembleImpl::arc_prob_matrix_t::size_type last_k=i;
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type k=1; k<=get_length(); ++k) {
    // 	    double p=prob_unpaired_external(k);
    // 	    if (p>threshold2) {
    // 		out << " " << (diff_encoding?(k-last_k):k);
    // 		if (write_probs) out << " " << p;
    // 		last_k=k;
    // 	    }
    // 	}
    // 	out << ";";
	
    // 	RnaEnsembleImpl::arc_prob_matrix_t::size_type last_ip=i;
    // 	RnaEnsembleImpl::arc_prob_matrix_t::size_type last_jp=j;
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type ip=1; ip<=get_length(); ++ip) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type jp=get_length(); jp>ip; --jp) {
    // 		if (pimpl_->arc_probs_(ip,jp)>threshold1) {
    // 		    double p=prob_basepair_external(ip,jp);
    // 		    if (p>threshold3) {
    // 			out << " " << (diff_encoding?(int)ip-(int)last_ip:(int)ip) << " " << (diff_encoding?(int)last_jp-(int)jp:(int)jp);
    // 			if (write_probs) out << " " << p;
    // 			last_ip=ip;
    // 			last_jp=jp;
    // 		    }
    // 		}
    // 	    }
    // 	}
    // 	out << std::endl;

	
    // 	// write lines for internal loops
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type i=1; i<=get_length(); ++i) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type j=i+1; j<=get_length(); ++j) {
		
    // 		if (pimpl_->arc_probs_(i,j)>threshold1) {
		    
    // 		    out << (diff_encoding?(int)i-(int)last_i:(int)i) << " " << (diff_encoding?(int)last_j-(int)j:(int)j);
    // 		    last_i=i; 
    // 		    last_j=j;
		    
    // 		    // write base pair and stacking probability
    // 		    out << " " << pimpl_->arc_probs_(i,j);
    // 		    if (pimpl_->arc_2_probs_(i,j)>threshold1) {
    // 			out << " " << pimpl_->arc_2_probs_(i,j);
    // 		    }
    // 		    out << " ;";
		    
    // 		    // write unpaired in loop
    // 		    RnaEnsembleImpl::arc_prob_matrix_t::size_type last_k=i;
    // 		    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type k=i+1; k<=j-1; ++k) {
    // 			double p=prob_unpaired_in_loop(k,i,j);
    // 			if (p>threshold2) {
    // 			    out << " " << (diff_encoding?(k-last_k):k);
    // 			    if (write_probs) out << " " << p;
    // 			    last_k=k;
    // 			}
    // 		    }
		    
    // 		    out << " ;";
		    
    // 		    RnaEnsembleImpl::arc_prob_matrix_t::size_type last_ip=i;
    // 		    RnaEnsembleImpl::arc_prob_matrix_t::size_type last_jp=j;
		    
    // 		    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type ip=i+1; ip<=j-1; ++ip) {
    // 			for(RnaEnsembleImpl::arc_prob_matrix_t::size_type jp=j-1; jp>ip ; --jp) {
    // 			    double p=prob_basepair_in_loop(ip,jp,i,j);
    // 			    if (p>threshold3) {
				
    // 				out << " " << (diff_encoding?(int)ip-(int)last_ip:(int)ip) << " " << (diff_encoding?(int)last_jp-(int)jp:(int)jp);
    // 				if (write_probs) out << " " << p;
    // 				last_ip=ip;
    // 				last_jp=jp;
    // 			    }
    // 			}
    // 		    }
    // 		    out << std::endl;
    // 		}
    // 	    }
    // 	}

    // 	return out;	
    // }


    // std::ostream &
    // RnaEnsemble::write_basepair_probs(std::ostream &out,double threshold) const {
    // 	for(RnaEnsembleImpl::arc_prob_matrix_t::size_type i=1; i<=get_length(); ++i) {
    // 	    for(RnaEnsembleImpl::arc_prob_matrix_t::size_type j=i+1; j<=get_length(); ++j) {
		
    // 		if (pimpl_->arc_probs_(i,j)>threshold) {
    // 		    out << i << " " << j;
		    
    // 		    // write base pair and stacking probability
    // 		    out << " " << pimpl_->arc_probs_(i,j);
    // 		    if (pimpl_->arc_2_probs_(i,j)>threshold) {
    // 			out << " " << pimpl_->arc_2_probs_(i,j);
    // 		    }
    // 		}
    // 	    }
    // 	}
    // 	return out;
    // }

    // std::ostream &
    // RnaEnsemble::write_pp(std::ostream &out,
    // 			  int width,
    // 			  double thresh1,
    // 			  double thresh2,
    // 			  double thresh3) const
    // {
    
    // 	size_type length=get_length();
    
    // 	// write sequence
    // 	for(size_type k=1; k<=length; k+=width) {
    // 	    pimpl_->sequence_.write( out, k, std::min(length,k+width-1) );
    // 	    out <<std::endl;
    // 	}

    // 	// write constraints
    // 	if (pimpl_->seq_constraints_ != "") {
    // 	    out << "#C "<<pimpl_->seq_constraints_<<std::endl;
    // 	}
    
    // 	// write separator
    // 	out << std::endl;
    // 	out << "#" << std::endl;
    
    // 	// write probabilities
    
    // 	for (size_type i=1; i<=length; i++) {
    // 	    for (size_type j=i+1; j<=length; j++) {
    // 		double p=arc_prob(i,j);
    // 		if (p > thresh1) {
    // 		    out << i << " " << j << " " << p;
		
    // 		    // write joint probability if above threshold 1
    // 		    double p2=arc_2_prob(i,j);
    // 		    if ( p2 > thresh1 ) {
    // 			out << " " << p2;
    // 		    }
		    
    // 		    // write positions of bases with unpaired in loop probabilities above threshold 2
		
    // 		    for (size_type k=i+1; k<j; ++k) {
    // 		    }
		    
    // 		    // write positions of base pairs with in loop probabilities above threshold 3
		    
    // 		    out << std::endl;
    // 		}
    // 	    } // end for j
    // 	} // end for i
	
    // 	std::cerr << "Warning: rna_ensemble::write_pp not fully implemented!"<<std::endl;

    // 	return out;
    // }


} // end namespace LocARNA
