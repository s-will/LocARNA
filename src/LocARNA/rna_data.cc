#include "sparse_rna_data_impl.hh"

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
		     const PFoldParams pfoldparams)
	: pimpl_(new RnaDataImpl(this,
				 filename,
				 p_bpcut,
				 pfoldparams)) {
    }


    RnaData::~RnaData() {
	delete pimpl_;
    }


    RnaDataImpl::RnaDataImpl(const RnaData &self,
			     const std::string &filename,
			     double p_bpcut,
			     const PFoldParams pfoldparams)
	:self_(self),
	 p_bpcut_(p_bpcut)
    {
	self_->read_autodetect(filename,pfoldparams,false);
    }


    ExtRnaDataImpl::ExtRnaDataImpl(const RnaData &self,
				   const std::string &filename,
				   double p_bpcut,
				   const PFoldParams pfoldparams)
	:self_(self),
	 p_bpcut_(p_bpcut),
	 in_loop_probs_available_(false)
    {
	self_->read_autodetect(filename,pfoldparams,true);
    }
    	    
    
    ExtRnaData::ExtRnaData(const RnaEnsemble &rna_ensemble,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut)
	: pimpl_(new ExtRnaDataImpl(this,
				    rna_ensemble,
				    p_bpcut,
				    p_bpilcut,
				    p_uilcut)) {
    }
    
    ExtRnaData::ExtRnaData(const std::string &filename,
			   double p_bpcut,
			   double p_bpilcut,
			   double p_uilcut)
	: pimpl_(new ExtRnaDataImpl(this,
				    filename,
				    double p_bpcut,
				    double p_bpilcut,
				    double p_uilcut)) {
    }

    ExtRnaData::~ExtRnaData() {
	delete pimpl_;
    }


    
    void
    RnaData::read_autodetect(const std::string &filename,
			     double p_bpcut,
			     const PFoldParams pfoldparams,
			     bool inloopprobs
			     ) {
	
	bool failed=false; //flag for signalling a failed attempt to
			   //read a certain file format
	
	bool recompute=false; // do we need to recompute probabilities

	try {
	    read_ps(filename, pfoldparams->stacking_);
	} catch (wrong_format_failure &f) {
	    failed=true;
	}
	
	if (failed) {
	    try {
		read_pp(filename, pfoldparams->stacking_);
	    } catch (wrong_format_failure &f) {
		failed=true;
	    }
	}

	try {
	    pimpl_->sequence_ = MultipleAlignment(filename, MultipleAlignment::FASTA);
	} catch (failure &f) {
	    failed=true;
	}
	if (failed) {
	    failed=false;
	    try {
		pimpl_->sequence_ = MultipleAlignment(filename, MultipleAlignment::CLUSTAL);
	    } catch (failure &f) {
		failed=true;
	    }
	}
	
	// if MultipleAlignment would support AUTO, we could replace the above by
	// if (failed) {
	//     pimpl_->sequence_ = MultipleAlignment(filename,MultipleAlignment::AUTO);
	//     recompute = true;
	// }
	
	// now, we have the sequence but not necessarily all required probabilities!
	if (!inloopprobs_ok()) {
	    recompute=true;
	}
	
	if (recompute) {
	    // recompute all probabilities
	    RnaData rna_data(pimpl_->sequence_,pfoldparams,!inloopprobs_ok(),true); // use given parameters, use alifold
	    // initialize from (temporary) rnadata object
	    init_from_rna_data(rna_data,pfoldparams);
	}
	
	return;
    }


    RnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
				    const PFoldParams &pfoldparams) {
	throw failure("not implemented");
    }
    
    ExtRnaData::init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
				       const PFoldParams &pfoldparams) {
	throw failure("not implemented");
    }


    const Sequence &
    RnaData::sequence() const {
	return pimpl_->sequence_;
    }

    size_type
    RnaData::length() const {
	return pimpl_->sequence_->length();
    }
    
    const std::string &
    RnaData::sequence_anchors() const {
	return pimpl_->sequence_anchors_;
    }
    
    double
    RnaData::arc_cutoff_prob() const {
	return pimpl_->arc_cutoff_prob_;
    }
    
    double 
    RnaData::arc_prob(pos_type i, pos_type j) const {
	return pimpl_->arc_probs_(i,j);
    }
    
    double
    RnaData::joint_arc_prob(pos_type i, pos_type j) const {
	return pimpl_->arc_2_probs(i,j);
    }
    
    double 
    RnaData::stacked_arc_prob(pos_type i, pos_type j) const {
	assert(pimpl_->arc_probs_(i+1,j-1)!=0);
	
	return
	    pimpl_->arc_2_probs(i,j)
	    /
	    pimpl_->arc_probs_(i+1,j-1);
    }

    double 
    RnaData::prob_paired_upstream(size_type i) const {
	double prob_paired=0.0;
	
	for (size_type j=i+1; j<=get_length(); j++) {
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
    
    void RnaData::read_ps(const std::string &filename, 
			  double p_incut) {
	
	std::ifstream in(filename.c_str());
	std::string line;
		
	getline(in,s);
	if (s!="%!PS-Adobe-3.0 EPSF-3.0") {
	    throw wrong_format_failure();
	}

	std::string seqname = "seq";
	
	pimpl_->stacking_probs_available_=false;
	
	while (getline(in,line) && has_prefix(line,"/sequence")) {
	    if (has_prefix(line,"% Probabilities for stacked pairs")) {
		stacking_probs_available_=true;
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
		
		    if (! (1<=i && i<j && j<=sequence_.length())) {
			std::cerr << "WARNING: Input dotplot "<<filename<<" contains invalid line " << line << " (indices out of range)" << std::endl;
		    } else {
			if (type=="ubox") {
			    pimpl_->arc_probs_.set(i,j,p);
			}
			else if (pimpl_->stacking_probs_available_ && type=="lbox") { // read a stacking probability
			    pimpl_->arc_2_probs_.set(i,j,p); // we store the joint probability of (i,j) and (i+1,j-1)
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
			seq_constraints_ += seqstr;
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
		p_cut_ = p;
	    }

	    in >> s;	    
	    if (has_prefix(s,"stack")) {
		stacking_probs_available_=true;
	    }
	    
	}

	for (std::map<std::string,std::string>::iterator it=seq_map.begin(); it!=seq_map.end(); ++it) {
	    // std::cout << "SEQ: " << it->first << " " << it->second << std::endl;
	    sequence_.append(Sequence::SeqEntry(it->first,it->second));
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
      
	    pimpl_->arc_probs_.set(i,j,p);
      
	    double p2;
	    
	    if (in >> p2) {
		pimpl_->arc_2_probs.set(i,j,p2); // p2 is joint prob of (i,j) and (i+1,j-1)
		pimpl_->stacking_probs_available_ = true;
	    }
	}
    } // end read_pp


    void
    RnaDataImpl::set_arc_probs_from_McCaskill_bppm(double threshold, bool stacking) {
	clear_arc_probs();
	
	for( size_t i=1; i <= sequence_.length(); i++ ) {
	    for( size_t j=i+1; j <= sequence_.length(); j++ ) {
		
		double p= McCmat_->get_bppm(i,j);
		
		if (p >= threshold) { // apply filter
		    set_arc_prob(i,j,p);
		}
	    }
	}
	if(stacking){
	    plist* pl= stackProb(threshold);
	    plist* p= pl;
	    while(pl->i!=0){
		set_arc_2_prob(pl->i, pl->j, pl->p);
		pl++;
	    }
	    free(p);
	    // 	  pl= stackProb(threshold);
	    // 	  free(pl);
	    pl= NULL;
	    p= NULL;
	}
	
	pair_probs_available_=true;
	stacking_probs_available_=stacking;
    }

} // end namespace LocARNA
