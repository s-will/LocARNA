#include <fstream>
#include <sstream>
#include <map>
#include <limits>

#include "aux.hh"
#include "rna_data.hh"

#include "alphabet.hh"

#include "multiple_alignment.hh"

#ifdef HAVE_LIBRNA
extern "C" {
//#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/loop_energies.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/pair_mat.h>
#include <ViennaRNA/alifold.h>

    FLT_OR_DBL *alipf_export_bppm(void);
}

#include "mcc_matrices.hh"

#endif // HAVE_LIBRNA


namespace LocARNA {


    // ------------------------------------------------------------
    // implementation of class RnaData
    //
    RnaData::RnaData(const std::string &file,
		     bool readPairProbs,
		     bool readStackingProbs,
		     bool readInLoopProbs)
	:
	sequence_(),
	pair_probs_available_(false),
	stacking_probs_available_(false),
	in_loop_probs_available_(false),
	arc_probs_(0),
	arc_2_probs_(0),
	seq_constraints_(""),
	McCmat_(0L),
	used_alifold_(false),
	min_free_energy_(std::numeric_limits<double>::infinity()),
	min_free_energy_structure_("")
    {
	init_from_file(file, readPairProbs, readStackingProbs, readInLoopProbs);
    }
        
    RnaData::RnaData(const Sequence &sequence)
	: sequence_(sequence),
	  pair_probs_available_(false),
	  stacking_probs_available_(false),
	  in_loop_probs_available_(false),
	  arc_probs_(0), // init sparse matrix with default 0
	  arc_2_probs_(0), // init sparse matrix with default 0
	  seq_constraints_(""), // empty sequence constraints string
	  McCmat_(0L), // 0 pointer
	  used_alifold_(false),
	  min_free_energy_(std::numeric_limits<double>::infinity()),
	  min_free_energy_structure_("")
    {
    }
    
    RnaData::~RnaData() {
#ifdef HAVE_LIBRNA
	if (McCmat_) {
	    delete McCmat_;
	}
#endif
    }

#ifdef HAVE_LIBRNA

    void
    RnaData::compute_ensemble_probs(const PFoldParams &params,bool inLoopProbs, bool use_alifold) {
	
	stopwatch.start("bpp");
	
	assert(use_alifold || sequence_.row_number()==1);

	used_alifold_=use_alifold;

	// run McCaskill and get access to results
	// in McCaskill_matrices
	if (!use_alifold) {
	    compute_McCaskill_matrices(params,inLoopProbs);
	} else {
	    make_pair_matrix();
	    compute_McCaskill_alifold_matrices(params,inLoopProbs);
	}
	
	// initialize the object from base pair probabilities
	// Use the same proability threshold as in RNAfold -p !
	set_arc_probs_from_McCaskill_bppm(10e-6,params.stacking);
	
	// since we either have local copies of all McCaskill pf arrays
	// or don't need them anymore,
	// we can free the ones of the Vienna lib
	if (!used_alifold_) {
	    free_pf_arrays();
	} else {
	    free_alipf_arrays();
	}

	pair_probs_available_=true;
	in_loop_probs_available_=inLoopProbs;

	stopwatch.stop("bpp");
    }
    
    void
    RnaData::compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs) {
	assert(sequence_.row_number()==1);

	// use MultipleAlignment to get pointer to c-string of the
	// first (and only) sequence in object sequence.
	//
	size_t length = sequence_.length();
	
	char c_sequence[length+1];
	std::string seqstring = MultipleAlignment(sequence_).seqentry(0).seq().to_string();
	strcpy(c_sequence,seqstring.c_str());
	
	char c_structure[length+1];
	
	// ----------------------------------------
	// set folding parameters
	if (params.noLP) {noLonelyPairs=1;}
	
	
	// ----------------------------------------
	// call fold for setting the pf_scale
	min_free_energy_ = fold(c_sequence,c_structure);
	min_free_energy_structure_ = c_structure;
	// std::cout << c_structure << std::endl;
	free_arrays();
	
	// set pf_scale
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	pf_scale = exp(-min_free_energy_/kT/length);

	// ----------------------------------------
	// call pf_fold
	pf_fold(c_sequence,c_structure);
	
	// ----------------------------------------
	// get McC data structures and copy
	// 
	// since the space referenced by pointers in McCmat will be
	// overwritten by the next call to pf_fold, we have to copy
	// the data structures if we want to keep them.
	//
	McCmat_ = new McC_matrices_t(c_sequence,inLoopProbs); // makes local copy (if needed)
	
	// precompute further tables (expMLbase, scale, qm2) for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	if (inLoopProbs) {
	    expMLbase_.resize(length+1);
	    scale_.resize(length+1);
	    
	    // ----------------------------------------
	    // from scale_pf_params
	    //
	    kT = McCmat_->pf_params->kT;   /* kT in cal/mol  */
	    
	    /* scaling factors (to avoid overflows) */
	    if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
		pf_scale = exp(-(-185+(McCmat_->pf_params->temperature-37.)*7.27)/kT);
		if (pf_scale<1) pf_scale=1;
	    }
	    
	    scale_[0] = 1.;
	    scale_[1] = 1./pf_scale;
	    expMLbase_[0] = 1;
	    expMLbase_[1] = McCmat_->pf_params->expMLbase * scale_[1];
	    for (size_t i=2; i<=sequence_.length(); i++) {
		scale_[i] = scale_[i/2]*scale_[i-(i/2)]; // scale_[i] = pow(scale_[1],(double)i)
		expMLbase_[i] = pow(McCmat_->pf_params->expMLbase, (double)i) * scale_[i];
	    }
	
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2();
	}
    }


    void
    RnaData::compute_McCaskill_alifold_matrices(const PFoldParams &params, bool inLoopProbs) {
	
	size_t length = sequence_.length();
	size_t n_seq = sequence_.row_number();
	
	// ----------------------------------------
	// write sequences to array of C-strings
	MultipleAlignment ma(sequence_);
	char **sequences = new char*[n_seq+1];
	for (size_t i=0; i<n_seq; i++) {
	    sequences[i]=new char[length+1];
	    std::string seqstring = ma.seqentry(i).seq().to_string();
	    strncpy(sequences[i],seqstring.c_str(),length+1);
	}
	sequences[n_seq]=NULL; //sequences has to be NULL terminated for alifold() etc

	const char **c_sequences=const_cast<const char **>(sequences);

	// reserve space for structure
	char *c_structure = new char[length+1];
	
	// ----------------------------------------
	// set folding parameters
	if (params.noLP) {noLonelyPairs=1;}
	
	// ----------------------------------------
	// call fold for setting the pf_scale
	min_free_energy_ = alifold(c_sequences,c_structure);
	min_free_energy_structure_ = c_structure;
	// std::cout << c_structure << std::endl;
	free_alifold_arrays();
	
	// set pf_scale
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	pf_scale = exp(-min_free_energy_/kT/length);
	
	
	// ----------------------------------------
	// call pf_fold
	alipf_fold(c_sequences,c_structure,NULL);
	
	// ----------------------------------------
	// get McC data structures and copy
	// 
	// since the space referenced by pointers in McCmat will be
	// overwritten by the next call to pf_fold, we have to copy
	// the data structures if we want to keep them.
	//
	McCmat_=new McC_ali_matrices_t(n_seq,length,inLoopProbs); // makes local copy (if needed)
	
	// precompute further tables (expMLbase, scale, qm2) for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	if (inLoopProbs) {
	    expMLbase_.resize(length+1);
	    scale_.resize(length+1);
	    
	    // ----------------------------------------
	    // from scale_pf_params
	    //
	    double scaling_factor=McCmat_->pf_params->pf_scale;
	    kT = McCmat_->pf_params->kT / n_seq;   /* kT in cal/mol  */
	    
	    
	    /* scaling factors (to avoid overflows) */
	    if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
		scaling_factor = exp(-(-185+(McCmat_->pf_params->temperature-37.)*7.27)/kT);
		if (scaling_factor<1) scaling_factor=1;
		McCmat_->pf_params->pf_scale=scaling_factor;
	    }
	    scale_[0] = 1.;
	    scale_[1] = 1./scaling_factor;

	    expMLbase_[0] = 1;
	    expMLbase_[1] = McCmat_->pf_params->expMLbase/scaling_factor;
	    for (size_t i=2; i<=sequence_.length(); i++) {
		scale_[i] = scale_[i/2]*scale_[i-(i/2)];
		expMLbase_[i] = pow(McCmat_->pf_params->expMLbase, (double)i) * scale_[i];
	    }
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2_ali();
	}

	//free c_sequences and c_structure
	for (size_t i=0; i<n_seq; i++) {
	    delete [] c_sequences[i];
	}
	delete [] c_sequences;
	delete [] c_structure;
    }
    
#endif // HAVE_LIBRNA

    void 
    RnaData::init_from_file(const std::string &filename,
			    bool readPairProbs,
			    bool readStackingProbs,
			    bool readInLoopProbs			       
			  ) {
	assert(!readStackingProbs || readPairProbs);
	assert(!readInLoopProbs || readPairProbs);

	std::ifstream in(filename.c_str());
	if (! in.is_open()) {
	    std::ostringstream err;
	    err << "Cannot open "<<filename<<" for reading.";
	    throw failure(err.str());
	}

	// set to true if probs become available
	in_loop_probs_available_=false;
	pair_probs_available_=false;
	stacking_probs_available_=false;

	std::string s;
	// read first line and decide about file-format
	in >> s;
	in.close();
	if (s == "%!PS-Adobe-3.0") {
	    // try reading as dot.ps file (as generated by RNAfold)
	    read_ps(filename,readPairProbs,readStackingProbs);
	} else if (s.substr(0,7) == "CLUSTAL" || s[0]=='>') {
	    //read to multiple alignment object
	    MultipleAlignment ma(filename,
				 (s[0]=='>')
				 ?MultipleAlignment::FASTA
				 :MultipleAlignment::CLUSTAL);
	    // convert to sequence
	    sequence_ = Sequence(ma);
	} else {
	    // try reading as PP-file (proprietary format, which is easy to read and contains pair probs)
	    read_pp(filename,readPairProbs,readStackingProbs,readInLoopProbs);
	}
	
	// DUMP for debugging
	//std::cout << arc_probs_ << std::endl;
	//std::cout << arc_2_probs_ << std::endl;
    }

    void RnaData::read_ps(const std::string &filename,
			  bool readPairProbs,
			  bool readStackingProbs
			  ) {
	assert(!readStackingProbs || readPairProbs);
		
	in_loop_probs_available_=false;
	pair_probs_available_=readPairProbs;
	
	bool contains_stacking_probs=false;
	
	std::ifstream in(filename.c_str());
    	
	std::string s;
	while (in >> s && s!="/sequence") {
	    if (s=="stacked") contains_stacking_probs=true;
	}

	stacking_probs_available_ = readStackingProbs && contains_stacking_probs;
	
	in >> s; in >> s;

	std::string seqstr="";
	while (in >> s && s!=")") {
	    s = s.substr(0,s.size()-1); // chop of last character
	    // cout << s <<endl;
	    seqstr+=s;
	}
	
	std::string seqname = seqname_from_filename(filename);
	
	//! sequence characters should be upper case, and 
	//! Ts translated to Us
	normalize_rna_sequence(seqstr);
        
	sequence_.append(Sequence::SeqEntry(seqname,seqstr));
            
	std::string line;
	
	// return when reading of pair probs is not wanted
	if (!readPairProbs) {return;}
	
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
			    set_arc_prob(i,j,p);
			}
			else if (readStackingProbs && contains_stacking_probs && type=="lbox") { // read a stacking probability
			    set_arc_2_prob(i,j,p); // we store the joint probability of (i,j) and (i+1,j-1)
			}
		    }
		}
	    }
	}
    }


    void RnaData::read_pp(const std::string &filename,
			  bool readPairProbs,
			  bool readStackingProbs,
			  bool readInLoopProbs
			  ) {
	
	assert(!readStackingProbs || readPairProbs);
	assert(!readInLoopProbs || readPairProbs);
	
	in_loop_probs_available_=false; // this is not (yet) supported!
	pair_probs_available_=readPairProbs;
	
	std::ifstream in(filename.c_str());
	
	std::string name;
	std::string seqstr;
	
    
	// ----------------------------------------
	// read sequence/alignment
    
	std::map<std::string,std::string> seq_map;
    
	std::string line;
    
	while (getline(in,line) && line!="#" ) {
	    if (line.length()>0 && line[0]!=' ') {
		std::istringstream in(line);
		in >> name >> seqstr;
	    
		normalize_rna_sequence(seqstr);
	    
		if (name != "SCORE:") { // ignore the (usually first) line that begins with SCORE:
		    if (name == "#C") {
			seq_constraints_ += seqstr;
		    } else {
			seq_map[name] += seqstr;
		    }
		}
	    }
	}
    
	for (std::map<std::string,std::string>::iterator it=seq_map.begin(); it!=seq_map.end(); ++it) {
	    // std::cout << "SEQ: " << it->first << " " << it->second << std::endl;
	    sequence_.append(Sequence::SeqEntry(it->first,it->second));
	}
    
	// return when reading of pair probs is not wanted
	if (!readPairProbs) {return;}
	
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
      
	    set_arc_prob(i,j,p);
      
	    if (readStackingProbs) {
		double p2;
	  
		if (in >> p2) {
		    set_arc_2_prob(i,j,p2); // p2 is joint prob of (i,j) and (i+1,j-1)
		    stacking_probs_available_ = true;
		}
	    }
	}
    }
    
    void RnaData::clear_arc_probs() {
	arc_probs_.clear();
	arc_2_probs_.clear();
	pair_probs_available_=false;
	stacking_probs_available_=false;
    }


#ifdef HAVE_LIBRNA
    void
    RnaData::set_arc_probs_from_McCaskill_bppm(double threshold, bool stacking) {
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

    void
    RnaData::compute_Qm2(){
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);

	size_type len = sequence_.length();
	

	std::vector<FLT_OR_DBL> qqm(len+2,0);
	std::vector<FLT_OR_DBL> qqm1(len+2,0);
	
	//qm1.resize((len+1)*(len+2)/2);
	qm2_.resize((len+1)*(len+2)/2);
	
	// initialize qqm1
	for (size_type i=1; i<=len; i++) {
	    qqm[i]=0;
	    qqm1[i]=0;
	}
	
	for(size_type j=TURN+2; j<=len; j++) {
	    // --------------------
	    // one column of Qm1, which will be needed in the calculation of Qm2  
	    for(size_type i=j-TURN-1; i>=1; i--) {
		char type=MCm->get_ptype(i,j);
		qqm[i]= qqm1[i]*expMLbase_[1];
		if(type) {
		    qqm[i] +=
			MCm->get_qb(i,j)
			* exp_E_MLstem(type,
				       (i>1)   ? MCm->S1[i-1] : -1, 
				       (j<len) ? MCm->S1[j+1] : -1,  
				       MCm->pf_params);
		}
		
		//qm1[McCmat->iidx(i,j)]=qqm[i];

		assert(qqm[i] <= MCm->get_qm(i,j));
		assert((!frag_len_geq(i,j-1,TURN+2)) || qqm1[i] <= MCm->get_qm(i,j-1));
	    }
	    	    
	    // --------------------
	    // calculates column "j" of the Qm2 matrix
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*(TURN+2)+1; i>=1; i--) {
		    qm2_[MCm->iidx(i,j)] = 0;
		    for(size_type k = i + TURN+1; (k+1)+TURN+1 <= j; k++) {
			qm2_[MCm->iidx(i,j)] +=
			    MCm->get_qm(i,k)*qqm[k+1];
		    }
		    assert(qm2_[MCm->iidx(i,j)] <= MCm->get_qm(i,j));
		}
	    }
	    	    
	    // --------------------
	    // swap qqm and qqm1 (in constant time)
	    qqm1.swap(qqm);
	}
    }


    void
    RnaData::compute_Qm2_ali(){
	assert(used_alifold_);
	assert(McCmat_);
	McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t *>(this->McCmat_);

	size_type len   = sequence_.length();
	size_type n_seq = sequence_.row_number();
	
	std::vector<FLT_OR_DBL> qqm(len+2,0);
	std::vector<FLT_OR_DBL> qqm1(len+2,0);
	std::vector<int> type(n_seq);
	
	qm2_.resize((len+1)*(len+2)/2);
	
	// initialize qqm1
	for (size_type i=1; i<=len; i++)
	    qqm1[i]=0;
	
	for(size_type j= TURN+2; j<=len; j++) {

	    // --------------------
	    // first, calculate one row of matrix Qm1, which is needed
	    // in the subsequent calculation of Qm2
	    //
	    for(size_type i=j-TURN-1; i>=1; i--) {
		
		// get base pair types for i,j of all sequences
		for (size_t s=0; s<n_seq; ++s) {
		    type[s] = pair[MCm->S[s][i]][MCm->S[s][j]];
		    if (type[s]==0) type[s]=7;
		}
		
		qqm[i]= qqm1[i]*expMLbase_[1];
		
		FLT_OR_DBL qbt1=1.0; // collects contribution "inner basepair of multiloop"
		for (size_t s=0; s<n_seq; s++) {
		    qbt1 *= exp_E_MLstem(type[s], 
					 i>1 ? MCm->S5[s][i] : -1,
					 j<len ? MCm->S3[s][j] : -1,
					 MCm->pf_params);
		}
		qqm[i] += MCm->get_qb(i,j) * qbt1;

	    }
	    
	    // --------------------
	    // calculate a row of the matrix Qm2
	    //
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*TURN-3; i>=1; i--) {
		    qm2_[MCm->iidx(i+1,j-1)] = 0;
		    for(size_type k = i+TURN+2; k< j-TURN-2; k++) {
			qm2_[MCm->iidx(i+1,j-1)] +=
			    MCm->get_qm(i+1,k)*qqm1[k+1];
		    }
		}
	    }
	    
	    // --------------------
	    // swap row qqm and qqm1 (in constant time)
	    qqm1.swap(qqm);
	}
    }

    int
    RnaData::ptype_of_admissible_basepair(size_type i,size_type j) const {
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);

    	int type = MCm->get_ptype(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if ((type==0)
	    || (((type==3)||(type==4))&&no_closingGU)
	    || (MCm->get_qb(i,j)==0.0)
	    || (get_arc_prob(i,j)==0.0))
	    {
		return 0;
	    }
	
	return type;
    }

    double RnaData::prob_unpaired_in_loop_ali(size_type k,size_type i,size_type j) const {
    	assert(frag_len_geq(i,j,TURN+2));
	assert(i<k);
	assert(k<j);
	
	if (!in_loop_probs_available_) return 1.0;
	
	McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t*>(this->McCmat_);
	
        size_t n_seq = sequence_.row_number();
	
	// immediately return 0.0 if i and j do not pair
	if (get_arc_prob(i,j)==0.0 || MCm->get_qb(i,j)==0.0) {return 0.0;}
	
	// get base pair types for i,j of all sequences
	std::vector<int> type(n_seq);
	
	for (size_t s=0; s<n_seq; ++s) {
	    type[s] = pair[MCm->S[s][i]][MCm->S[s][j]];
	    if (type[s]==0) type[s]=7;
	}

	// ------------------------------------------------------------
	// hairpin contribution
        //
	
	FLT_OR_DBL H=1.0;
	
	for (size_t s=0; s<n_seq; s++) {
	    size_t u = MCm->a2s[s][j-1]-MCm->a2s[s][i];
	    if (MCm->a2s[s][i]<1) continue;
	    char loopseq[10];
	    if (u<7){
		strncpy(loopseq, MCm->Ss[s]+MCm->a2s[s][i]-1, 10);
	    }
	    H *= exp_E_Hairpin(u, type[s],
			       MCm->S3[s][i], MCm->S5[s][j], 
			       loopseq, 
			       MCm->pf_params);
        }
        H *= scale_[j-i+1];
    
	// ------------------------------------------------------------
	// interior loop contributions
	//
	
	FLT_OR_DBL I = 0.0;
	
	// case 1: i<k<i´<j´<j
	for (size_t ip=k+1; ip <= std::min(i+MAXLOOP+1,j-TURN-2); ip++) {
	    for (size_t jp = std::max(ip+TURN+1 + MAXLOOP ,j-1 + ip-i-1 ) - MAXLOOP; jp<j; jp++) {
		
		FLT_OR_DBL qloop=1.0;
		
		if (MCm->get_qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = MCm->a2s[s][ip-1] - MCm->a2s[s][i];
		    size_t u2 = MCm->a2s[s][j-1] - MCm->a2s[s][jp];
		    
		    int type_2 = pair[MCm->S[s][jp]][MCm->S[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    MCm->S3[s][i],
					    MCm->S5[s][j],
					    MCm->S5[s][ip],
					    MCm->S3[s][jp],
					    MCm->pf_params
					    );
		}
		
		I += MCm->get_qb(ip,jp) * scale_[ip-i+j-jp] * qloop;
	    }
	}
 
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip <= std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    for (size_t jp = std::max(ip+TURN+1 + MAXLOOP,j-1+ ip-i-1) - MAXLOOP; jp<k; jp++) {
		
		FLT_OR_DBL qloop=1.0;
		
		if (MCm->get_qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = MCm->a2s[s][ip-1] - MCm->a2s[s][i];
		    size_t u2 = MCm->a2s[s][j-1] - MCm->a2s[s][jp];
		    
		    int type_2 = pair[MCm->S[s][jp]][MCm->S[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    MCm->S3[s][i],
					    MCm->S5[s][j],
					    MCm->S5[s][ip],
					    MCm->S3[s][jp],
					    MCm->pf_params
					    );
		}
		
		I += MCm->get_qb(ip,jp) * scale_[ip-i+j-jp] * qloop;
	    }
	}
    
    
	// ------------------------------------------------------------
	// multiloop contributions
        //

	FLT_OR_DBL M = 0.0;
	
	// no base pair <= k:   i....k-----qm2-------j
	// valid entries of qm2_ have space for 2 inner base pairs,
	// i.e. at least length of "(...)(...)" (for TURN=3)
	if ( frag_len_geq(k+1, j-1, 2*(TURN+2)) ) {
	    M += qm2_[MCm->iidx(k+1,j-1)] * expMLbase_[k-i];
	}
	
	// no base pair >= k
	if ( frag_len_geq(i+1,k-1,2*(TURN+2)) ) {
	    M += qm2_[MCm->iidx(i+1,k-1)] * expMLbase_[j-k];
	}
	
	// base pairs <k and >k
	if ( frag_len_geq(i+1,k-1,TURN+2) && frag_len_geq(k+1,j-1,TURN+2) ) {
	    M += MCm->get_qm(i+1,k-1) * expMLbase_[1] *  MCm->get_qm(k+1,j-1);
	}
	
	// multiply with contribution for closing of multiloop
	
	for (size_t s=0; s<n_seq; s++) {
	    int tt = rtype[type[s]];
	    
	    M *= MCm->pf_params->expMLclosing 
		* exp_E_MLstem(tt,MCm->S5[s][j],MCm->S3[s][i], MCm->pf_params);
	}
	M *= scale_[2];
	
	FLT_OR_DBL Qtotal=H+I+M;
	
	double kTn   = MCm->pf_params->kT/10.;   /* kT in cal/mol  */
	
	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	Qtotal *= exp(MCm->get_pscore(i,j)/kTn);
	
	FLT_OR_DBL p_k_cond_ij = Qtotal/MCm->get_qb(i,j); 
	
	FLT_OR_DBL res = p_k_cond_ij * MCm->get_bppm(i,j);
	
	return res;
    }
    
    double RnaData::prob_unpaired_in_loop(size_type k,size_type i,size_type j) const {
	assert(i+TURN+1 <= j);
	assert(i<k);
	assert(k<j);
	
	if (!in_loop_probs_available_) return 1.0;
	
	
	if (used_alifold_) {
	    return prob_unpaired_in_loop_ali(k, i, j);
	}
	
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);

	const char *c_sequence = MCm->sequence;
	
	
	int type = ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}

	// ------------------------------------------------------------
	// Hairpin loop energy contribution
	
	size_t u=j-i-1;
	FLT_OR_DBL H = exp_E_Hairpin(u, type, MCm->S1[i+1], MCm->S1[j-1],
				     c_sequence+i-1, MCm->pf_params) * scale_[u+2];
	
	// ------------------------------------------------------------
	// Interior loop energy contribution
	FLT_OR_DBL I = 0.0;

	// case 1: i<k<i´<j´<j
	for (size_t ip=k+1; ip <= std::min(i+MAXLOOP+1,j-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp = std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<j; jp++) {
		int type2 = MCm->get_ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += MCm->get_qb(ip,jp) 
			* (scale_[u1+j-jp+1] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 MCm->S1[i+1],MCm->S1[j-1],
					 MCm->S1[ip-1],MCm->S1[jp+1], MCm->pf_params));
		}
	    }
	}
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip<=std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp=std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<k; jp++) {
		int type2 = MCm->get_ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += MCm->get_qb(ip,jp)
			* (scale_[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 MCm->S1[i+1],MCm->S1[j-1],
					 MCm->S1[ip-1],MCm->S1[jp+1], MCm->pf_params));
		}
	    }
	}
	
	// ------------------------------------------------------------
	// Multiple loop energy contribution
	FLT_OR_DBL M = 0.0;

	FLT_OR_DBL M1=0.0;
	FLT_OR_DBL M2=0.0;
	FLT_OR_DBL M3=0.0;

	// bases <=k unpaired
	if ( frag_len_geq(k+1, j-1, 2*(TURN+2)) ) {
	    M1 = expMLbase_[frag_len(i+1,k)] * qm2_[MCm->iidx(k+1,j-1)];
	}
	
	// bases >=k unpaired
	if ( frag_len_geq(i+1, k-1, 2*(TURN+2)) ) {
	    M2 = qm2_[MCm->iidx(i+1,k-1)] * expMLbase_[frag_len(k,j-1)];
	}
	
	// innner base pairs left and right of k
	if ( frag_len_geq(i+1,k-1,TURN+2) && frag_len_geq(k+1,j-1,TURN+2) ) {
	    M3 = MCm->get_qm(i+1,k-1) * expMLbase_[1] *  MCm->get_qm(k+1,j-1);
	}
	
	M=M1+M2+M3;
	
	// multiply with contribution for closing of multiloop
	M *= MCm->pf_params->expMLclosing 
	    * exp_E_MLstem(rtype[type],MCm->S1[j-1],MCm->S1[i+1], MCm->pf_params)
	    * scale_[2];
	
	FLT_OR_DBL Qtotal = H+I+M;

	FLT_OR_DBL p_k_cond_ij = Qtotal/MCm->get_qb(i,j); 

	FLT_OR_DBL res = p_k_cond_ij * MCm->get_bppm(i,j);
	
	return res;
    }

    double RnaData::prob_unpaired_external(size_type k) const {
	assert(1<=k);
	assert(k<=sequence_.length());
	
	if (!in_loop_probs_available_) return 1.0;
	
	return (McCmat_->q1k[k-1] * scale_[1] * McCmat_->qln[k+1]) / McCmat_->qln[1];
    }


    double
    RnaData::prob_basepair_in_loop_ali(size_type ip,
				       size_type jp,
				       size_type i,
				       size_type j) const {

	if (!in_loop_probs_available_) return 1.0;
	
	McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t *>(this->McCmat_);
	
	size_t n_seq = sequence_.row_number();
	
	// note: the following tests cover the case that the distances of i,j and ip,jp are too small

	// immediately return 0.0 if i and j do not pair
	if (get_arc_prob(i,j)==0.0 || MCm->get_qb(i,j)==0.0) {return 0.0;}
	
	// immediately return 0.0 when ip and jp cannot pair
	if (get_arc_prob(ip,jp)==0.0 || MCm->get_qb(ip,jp)==0.0) {return 0.0;}
	
	assert(frag_len_geq(i,j,TURN+4));
	assert(frag_len_geq(ip,jp,TURN+2));
	assert(i<ip);
	assert(jp<j);
	
	// ------------------------------------------------------------
	// get base pair types
	//
	std::vector<int> type(n_seq);
	std::vector<int> type2(n_seq);
	
	for (size_t s=0; s<n_seq; ++s) {
	    type[s] = pair[MCm->S[s][i]][MCm->S[s][j]];
	    if (type[s]==0) type[s]=7;

	    type2[s] = pair[MCm->S[s][ip]][MCm->S[s][jp]];
	    if (type2[s]==0) type2[s]=7;
	}


	// note: I and M are computed without factor get_qb(ip,jp),
	// which is multiplied only in the end.
	
	// ------------------------------------------------------------
	// Interior loop energy contribution
	//
	FLT_OR_DBL I=0.0;
	
	if ((frag_len(i,ip)+frag_len(jp,j))<=MAXLOOP) {
	    I = 1.0;
	    for (size_t s=0; s<n_seq; s++) {
		
		size_t u1 = MCm->a2s[s][ip-1] - MCm->a2s[s][i];
		size_t u2 = MCm->a2s[s][j-1] - MCm->a2s[s][jp];
		
		I *= exp_E_IntLoop(u1,u2,
				   type[s], rtype[type2[s]],
				   MCm->S3[s][i],
				   MCm->S5[s][j],
				   MCm->S5[s][ip],
				   MCm->S3[s][jp],
				   MCm->pf_params);
	    }
	    I *= scale_[ip-i+j-jp];
	}
	
	// ------------------------------------------------------------
	// Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase_[frag_len(i+1,ip-1)] * MCm->get_qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += MCm->get_qm(i+1,ip-1) * expMLbase_[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += MCm->get_qm(i+1,ip-1) * MCm->get_qm(jp+1,j-1);
	}
	
	for (size_t s=0; s<n_seq; s++) {
	    // multiply with factor for inner base pair
	    M *= exp_E_MLstem(type2[s],
			      MCm->S5[s][ip],
			      MCm->S3[s][jp],
			      MCm->pf_params);
	    
	    // multiply with factors for closing base pair
	    M *= MCm->pf_params->expMLclosing
		* exp_E_MLstem(rtype[type[s]],
			       MCm->S5[s][j],
			       MCm->S3[s][i],
			       MCm->pf_params);
	}
	
	M *= scale_[2]; // scale for closing base pair
	
	// ------------------------------------------------------------
	FLT_OR_DBL Qtotal= (I+M);
	
	Qtotal *= MCm->get_qb(ip,jp);

	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	double kTn   = MCm->pf_params->kT/10.;   /* kT in cal/mol  */
	Qtotal *= exp(MCm->get_pscore(i,j)/kTn);

	
	return
	    (Qtotal/MCm->get_qb(i,j))
	    *MCm->get_bppm(i,j);
    }

    double
    RnaData::prob_basepair_in_loop(size_type ip,
				   size_type jp,
				   size_type i,
				   size_type j) const {
	
	if (!in_loop_probs_available_) return 1.0;

	if (used_alifold_) {
	    return prob_basepair_in_loop_ali(ip, jp, i, j);
	}
	
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);
	

	// note: I and M are computed without factor get_qb(ip,jp),
	// which is multiplied only in the end.
	
	int type=ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}
	
	int type2=ptype_of_admissible_basepair(ip,jp);
	
	// immediately return 0.0 when ip and jp cannot pair
	if (type2==0) {return 0.0;}
	
	assert(frag_len_geq(i,j,TURN+4));
	assert(frag_len_geq(ip,jp,TURN+2));
	assert(i<ip);
	assert(jp<j);

	//calculating the Interior loop energy contribution
	//
	FLT_OR_DBL I=0.0;
	
	int u1 =(int)(ip-i-1);
	int u2 =(int)(j-jp-1);
	
	if (u1+u2 <= MAXLOOP) {
	    I = exp_E_IntLoop(u1,u2, type, rtype[type2],
			      MCm->S1[(int)(i+1)],
			      MCm->S1[(int)(j-1)],
			      MCm->S1[ip-1],
			      MCm->S1[jp+1],
			      MCm->pf_params)
		* scale_[u1+u2+2];
	}
	
	//calculating Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	
	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase_[frag_len(i+1,ip-1)] * MCm->get_qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += MCm->get_qm(i+1,ip-1) * expMLbase_[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += MCm->get_qm(i+1,ip-1) * MCm->get_qm(jp+1,j-1);
	}
	
	// multiply with factor for inner base pair
	M *= exp_E_MLstem(type2,
			  MCm->S1[ip-1],
			  MCm->S1[jp+1],
			  MCm->pf_params);
	
	// multiply with factors for closing base pair
	M *= MCm->pf_params->expMLclosing
	    * exp_E_MLstem(rtype[type],
			   MCm->S1[j-1],
			   MCm->S1[i+1],
			   MCm->pf_params)
	    * scale_[2];
	
	FLT_OR_DBL Qtotal = I+M;
	
	Qtotal *= MCm->get_qb(ip,jp);
	
	
	return Qtotal/MCm->get_qb(i,j)
	    * MCm->get_bppm(i,j);
    }

    double RnaData::prob_basepair_external(size_type i,size_type j) const {

	if (!in_loop_probs_available_) return 1.0;
	
	size_t n=sequence_.length();

	assert(1<=i);
	assert(i<j);
	assert(j<=n);
	assert(frag_len_geq(i,j,TURN+2));
	
	// immediately return 0.0 when i and j cannot pair
	if (get_arc_prob(i,j)==0.0 || McCmat_->get_qb(i,j)==0.0) {
	    return 0.0;
	}
	
	FLT_OR_DBL extloop;
	
	if (!used_alifold_) {
	    McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);
	    extloop = exp_E_ExtLoop(MCm->get_ptype(i,j),
				    i>1 ? MCm->S1[i-1] : -1, 
				    j<n ? MCm->S1[j+1] : -1, 
				    MCm->pf_params);
	} else {
	    McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t *>(this->McCmat_);
	    
	    size_t n_seq=sequence_.row_number();
	    
	    extloop=1.0;
	    
	    for (size_t s=0; s<n_seq; s++) {
		int type = pair[MCm->S[s][i]][MCm->S[s][j]];
		if (type==0) type=7;
		
		extloop *= exp_E_ExtLoop(type, i>1 ? MCm->S5[s][i] : -1,
					 j<n ? MCm->S3[s][j] : -1,
					 MCm->pf_params);
	    }
	}

	return
	    (McCmat_->q1k[i-1]
	     * McCmat_->get_qb(i,j)
	     * extloop
	     * McCmat_->qln[j+1]
	     )
	    / McCmat_->qln[1];

    }
    
    std::ostream &
    RnaData::write_unpaired_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
	
	// write lines for loops closed by base pairs
	for(arc_prob_matrix_t::size_type i=1; i<=sequence_.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence_.length(); ++j) {
		if (arc_probs_(i,j)>threshold1) {
		    bool had_entries=false;
		    for(arc_prob_matrix_t::size_type k=i+1; k<=j-1; ++k) {
			double p=prob_unpaired_in_loop(k,i,j);
			if (p>threshold2) {
			    if (!had_entries) {out << i << " " << j; had_entries=true;}
			    out << " " << k << " " << p;
			}
		    }
		    if (had_entries) out << std::endl;
		}
	    }
	}
	
	// write lines for external loop
	
	bool had_entries=false;
	for(arc_prob_matrix_t::size_type k=1; k<=sequence_.length(); ++k) {
	    double p=prob_unpaired_external(k);
	    if (p>threshold2) {
		if (!had_entries) {out << 0 << " " << (sequence_.length()+1); had_entries=true;}
		out << " " << k << " " << p;
	    }
	}
	if (had_entries) out << std::endl;
    	return out;
    }
    
	
    std::ostream &
    RnaData::write_basepair_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
	// write lines for loops closed by base pairs
	for(arc_prob_matrix_t::size_type i=1; i<=sequence_.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence_.length(); ++j) {
		if (arc_probs_(i,j)>threshold1) {
		    bool had_entries=false;
		    for(arc_prob_matrix_t::size_type ip=i+1; ip<=j-1; ++ip) {
			for(arc_prob_matrix_t::size_type jp=ip+1; jp<=j-1; ++jp) {
			    double p=prob_basepair_in_loop(ip,jp,i,j);
			    if (p>threshold2) {
				if (!had_entries) {out << i << " " << j; had_entries=true;}
				out << " " << ip << " " << jp << " " << p;
			    }
			}
		    }
		    if (had_entries) out << std::endl;
		}
	    }
	}
	
	// write lines for external loop
	bool had_entries=false;
	for(arc_prob_matrix_t::size_type ip=1; ip<=sequence_.length(); ++ip) {
	    for(arc_prob_matrix_t::size_type jp=ip+1; jp<=sequence_.length(); ++jp) {
		double p=prob_basepair_external(ip,jp);
		if (p>threshold2) {
		    if (!had_entries) {out << 0 << " " << (sequence_.length()+1); had_entries=true;}
		    out << " " << ip << " " << jp << " " << p;
		}
	    }
	}
	if (had_entries) out << std::endl;
    	return out;
    }
    
    
    std::ostream &
    RnaData::write_basepair_and_in_loop_probs(std::ostream &out,double threshold1,double threshold2,double threshold3, bool write_probs, bool diff_encoding) const {
	
	size_t i=0;
	size_t j=sequence_.length()+1;

	arc_prob_matrix_t::size_type last_i=i;
	arc_prob_matrix_t::size_type last_j=j;
	

	// write line for external loop
	
	out << (diff_encoding?(int)i-(int)last_i:(int)i) << " " << (diff_encoding?(int)last_j-(int)j:(int)j)<< " 1 ;";
	
	arc_prob_matrix_t::size_type last_k=i;
	for(arc_prob_matrix_t::size_type k=1; k<=sequence_.length(); ++k) {
	    double p=prob_unpaired_external(k);
	    if (p>threshold2) {
		out << " " << (diff_encoding?(k-last_k):k);
		if (write_probs) out << " " << p;
		last_k=k;
	    }
	}
	out << ";";
	
	arc_prob_matrix_t::size_type last_ip=i;
	arc_prob_matrix_t::size_type last_jp=j;
	for(arc_prob_matrix_t::size_type ip=1; ip<=sequence_.length(); ++ip) {
	    for(arc_prob_matrix_t::size_type jp=sequence_.length(); jp>ip; --jp) {
		if (arc_probs_(ip,jp)>threshold1) {
		    double p=prob_basepair_external(ip,jp);
		    if (p>threshold3) {
			out << " " << (diff_encoding?(int)ip-(int)last_ip:(int)ip) << " " << (diff_encoding?(int)last_jp-(int)jp:(int)jp);
			if (write_probs) out << " " << p;
			last_ip=ip;
			last_jp=jp;
		    }
		}
	    }
	}
	out << std::endl;

	
	// write lines for internal loops
	for(arc_prob_matrix_t::size_type i=1; i<=sequence_.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence_.length(); ++j) {
		
		if (arc_probs_(i,j)>threshold1) {
		    
		    out << (diff_encoding?(int)i-(int)last_i:(int)i) << " " << (diff_encoding?(int)last_j-(int)j:(int)j);
		    last_i=i; 
		    last_j=j;
		    
		    // write base pair and stacking probability
		    out << " " << arc_probs_(i,j);
		    if (arc_2_probs_(i,j)>threshold1) {
			out << " " << arc_2_probs_(i,j);
		    }
		    out << " ;";
		    
		    // write unpaired in loop
		    arc_prob_matrix_t::size_type last_k=i;
		    for(arc_prob_matrix_t::size_type k=i+1; k<=j-1; ++k) {
			double p=prob_unpaired_in_loop(k,i,j);
			if (p>threshold2) {
			    out << " " << (diff_encoding?(k-last_k):k);
			    if (write_probs) out << " " << p;
			    last_k=k;
			}
		    }
		    
		    out << " ;";
		    
		    arc_prob_matrix_t::size_type last_ip=i;
		    arc_prob_matrix_t::size_type last_jp=j;
		    
		    for(arc_prob_matrix_t::size_type ip=i+1; ip<=j-1; ++ip) {
			for(arc_prob_matrix_t::size_type jp=j-1; jp>ip ; --jp) {
			    double p=prob_basepair_in_loop(ip,jp,i,j);
			    if (p>threshold3) {
				
				out << " " << (diff_encoding?(int)ip-(int)last_ip:(int)ip) << " " << (diff_encoding?(int)last_jp-(int)jp:(int)jp);
				if (write_probs) out << " " << p;
				last_ip=ip;
				last_jp=jp;
			    }
			}
		    }
		    out << std::endl;
		}
	    }
	}

	return out;	
    }


#endif // HAVE_LIBRNA

    std::ostream &
    RnaData::write_basepair_probs(std::ostream &out,double threshold) const {
	for(arc_prob_matrix_t::size_type i=1; i<=sequence_.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence_.length(); ++j) {
		
		if (arc_probs_(i,j)>threshold) {
		    out << i << " " << j;
		    
		    // write base pair and stacking probability
		    out << " " << arc_probs_(i,j);
		    if (arc_2_probs_(i,j)>threshold) {
			out << " " << arc_2_probs_(i,j);
		    }
		}
	    }
	}
	return out;
    }
    
    std::string RnaData::seqname_from_filename(const std::string &s) const {
	size_type i;
	size_type j;
    
	assert(s.length()>0);
    
	for (i=s.length(); i>0 && s[i-1]!='/'; i--)
	    ;

	for (j=i; j<s.length() && s[j]!='.'; j++)
	    ;

	std::string name=s.substr(i,j-i);

	if (name.length()>3 && name.substr(name.length()-3,3) == "_dp") {
	    return name.substr(0,name.length()-3);
	}
	
	return name;
    }

    std::ostream &
    RnaData::write_pp(std::ostream &out,
		      int width,
		      double thresh1,
		      double thresh2,
		      double thresh3) const
    {
    
	size_type length=sequence_.length();
    
	// write sequence
	for(size_type k=1; k<=length; k+=width) {
	    sequence_.write( out, k, std::min(length,k+width-1) );
	    out <<std::endl;
	}

	// write constraints
	if (seq_constraints_ != "") {
	    out << "#C "<<seq_constraints_<<std::endl;
	}
    
	// write separator
	out << std::endl;
	out << "#" << std::endl;
    
	// write probabilities
    
	for (size_type i=1; i<=length; i++) {
	    for (size_type j=i+1; j<=length; j++) {
		double p=get_arc_prob(i,j);
		if (p > thresh1) {
		    out << i << " " << j << " " << p;
		
		    // write joint probability if above threshold 1
		    double p2=get_arc_2_prob(i,j);
		    if ( p2 > thresh1 ) {
			out << " " << p2;
		    }
		    
		    // write positions of bases with unpaired in loop probabilities above threshold 2
		
		    for (size_type k=i+1; k<j; ++k) {
		    }
		    
		    // write positions of base pairs with in loop probabilities above threshold 3
		    
		    out << std::endl;
		}
	    } // end for j
	} // end for i
	
	std::cerr << "Warning: rna_data::write_pp not fully implemented!"<<std::endl;

	return out;
    }
    
}
    
