#include <fstream>
#include <sstream>
#include <map>

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

#endif // HAVE_LIBRNA



// // for getrusage()
// #include <sys/resource.h>
// #include <sys/types.h>
// // for gettimeofday()
// #include <sys/time.h>

namespace LocARNA {

    // ------------------------------------------------------------
    // implementation of classes McC_matrices_base, McC_matrices_t, McC_ali_matrices_t
#ifdef HAVE_LIBRNA

    void *
    space_memcpy(void *from,size_t size) {
	if (from==NULL) return from;
	void *p = space(size);
	memcpy(p,from,size);
	return p;
    }

    McC_matrices_base::McC_matrices_base()
	: length(0),local_copy(false),qb(0),qm(0),bppm(0),iindx(0),q1k(0),qln(0)	
    {}

    void
    McC_matrices_base::init(size_t length_) {
	length=length_;
	
	qb=0;
	qm=0;
	bppm=0;
	q1k=0;
	qln=0;
	
	iindx = get_iindx(length);
    }
    
    void
    McC_matrices_base::deep_copy(const McC_matrices_base &McCmat) {
	local_copy=true;

	length=McCmat.length;

	size_t size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
	
	qb= (FLT_OR_DBL *) space_memcpy(McCmat.qb,size);
	qm= (FLT_OR_DBL *) space_memcpy(McCmat.qm,size);
	bppm= (FLT_OR_DBL *) space_memcpy(McCmat.bppm,size);
	q1k= (FLT_OR_DBL *) space_memcpy(McCmat.q1k,sizeof(FLT_OR_DBL)*(length+1));
	qln= (FLT_OR_DBL *) space_memcpy(McCmat.qln,sizeof(FLT_OR_DBL)*(length+2));
	pf_params= (pf_paramT *) space_memcpy(McCmat.pf_params,sizeof(pf_paramT));

	iindx= get_iindx(length);
    }

    McC_matrices_base::~McC_matrices_base() {
	if (local_copy) {
	    free_all();
	} else {
	    free(iindx);
	}
    }

    void McC_matrices_base::free_all() {
	if (qb) free(qb);
	if (qm) free(qm);
	if (q1k) free(q1k);
	if (qln) free(qln);
	if (iindx) free(iindx);
	if (pf_params) free(pf_params);
    }
    
    // ----------------------------------------

    McC_matrices_t::McC_matrices_t(char *sequence, bool local_copy) {
	
	if (local_copy) {
	    McC_matrices_t McCmat_tmp(sequence,false);
	    deep_copy(McCmat_tmp);
	} else {
	    McC_matrices_base::init(strlen(sequence));
	    
	    this->sequence=sequence;

	    // get pointers to McCaskill matrices
	    get_pf_arrays(&S,
			  &S1,
			  &ptype,
			  &qb,
			  &qm,
			  &q1k,
			  &qln);
	    
	    // get pointer to McCaskill base pair probabilities
	    bppm = export_bppm();	    
	    
	    pf_params = get_scaled_pf_parameters();
	}
    }

    void
    McC_matrices_t::deep_copy(const McC_matrices_t &McCmat) {
	McC_matrices_base::deep_copy(McCmat);
	
	sequence = (char *) space_memcpy(McCmat.sequence,sizeof(char)*(length+1));
	S = (short *) space_memcpy(McCmat.S,sizeof(short)*(length+2));
	S1 = (short *) space_memcpy(McCmat.S1,sizeof(short)*(length+2));
	ptype= (char *) space_memcpy(McCmat.ptype,sizeof(char)*((length+1)*(length+2)/2));
    }

    McC_matrices_t::~McC_matrices_t() {
	if (local_copy) {
	    free_all();
	}
    }


    void McC_matrices_t::free_all() {
	free(sequence);
	free(S);
	free(S1);
	free(ptype);
    }

    // ----------------------------------------
    McC_ali_matrices_t::McC_ali_matrices_t(size_t n_seq_, size_t length_, bool local_copy_)
	: n_seq(n_seq_)
    {	
	if (local_copy_) {
	    McC_ali_matrices_t McCmat_tmp(n_seq,length_,false);
	    deep_copy(McCmat_tmp);
	} else {
	    McC_matrices_base::init(length_);
	    
	    // get pointers to McCaskill matrices
	    get_alipf_arrays(&S,
			     &S5,
			     &S3,
			     &a2s,
			     &Ss,
			     &qb,
			     &qm,
			     &q1k,
			     &qln,
			     &pscore);

	    // get pointer to McCaskill base pair probabilities
	    bppm = alipf_export_bppm();	    
	    
	    pf_params = get_scaled_alipf_parameters(n_seq);
	}
    }
    

    void
    McC_ali_matrices_t::deep_copy(const McC_ali_matrices_t &McCmat) {
	McC_matrices_base::deep_copy(McCmat);
		
	n_seq = McCmat.n_seq;

	S    = (short **)          space(n_seq * sizeof(short *));
	S5   = (short **)          space(n_seq * sizeof(short *));
	S3   = (short **)          space(n_seq * sizeof(short *));
	a2s  = (unsigned short **) space(n_seq * sizeof(unsigned short *));
	Ss   = (char **)           space(n_seq * sizeof(char *));

	for (size_t i=0; i<n_seq; i++) {
	    S[i]   = (short *)          space_memcpy(McCmat.S[i],  (length+2) * sizeof(short));
	    S5[i]  = (short *)          space_memcpy(McCmat.S5[i], (length+2) * sizeof(short));
	    S3[i]  = (short *)          space_memcpy(McCmat.S3[i], (length+2) * sizeof(short));
	    a2s[i] = (unsigned short *) space_memcpy(McCmat.a2s[i],(length+2) * sizeof(unsigned short));
	    Ss[i]  = (char *)           space_memcpy(McCmat.Ss[i], (length+2) * sizeof(char));
	}

	pscore = (short*) space_memcpy(McCmat.pscore,
				       ((length+1)*(length+2))/2 * sizeof(short));
    }

    McC_ali_matrices_t::~McC_ali_matrices_t() {
	if (local_copy) {
	    free_all();
	}
    }


    void McC_ali_matrices_t::free_all() {
	free_sequence_arrays(n_seq,&S,&S5,&S3,&a2s,&Ss);
	if (pscore) free(pscore);
    }

#   endif
    

    // ------------------------------------------------------------
    // implementation of class RnaData
    //
    RnaData::RnaData(const std::string &file,
		     bool readPairProbs,
		     bool readStackingProbs,
		     bool readInLoopProbs)
	:
	sequence(),
	arc_probs_(0),
	arc_2_probs_(0),
	seq_constraints_(""),
	used_alifold(false)
    {
	initFromFile(file, readPairProbs, readStackingProbs, readInLoopProbs);
    }
        
    RnaData::RnaData(const Sequence &sequence_)
	: sequence(sequence_),
	  pair_probs_available(false),	  
	  arc_probs_(0),
	  arc_2_probs_(0),
	  seq_constraints_(""),
	  used_alifold(false)
    {
    }
    
    RnaData::~RnaData() {
#ifdef HAVE_LIBRNA
	if (McCmat) {delete McCmat;} 
#endif //HAVE_LIBRNA
    }

#ifdef HAVE_LIBRNA

    void
    RnaData::computeEnsembleProbs(const PFoldParams &params,bool inLoopProbs, bool use_alifold) {
	
	stopwatch.start("bpp");
	
	assert(use_alifold || sequence.row_number()==1);

	used_alifold=use_alifold;

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
	if (!used_alifold) {
	    free_pf_arrays();
	} else {
	    free_alipf_arrays();
	}

	pair_probs_available=true;
	in_loop_probs_available=inLoopProbs;

	stopwatch.stop("bpp");
    }
    
    void
    RnaData::compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs) {
	assert(sequence.row_number()==1);

	// use MultipleAlignment to get pointer to c-string of the
	// first (and only) sequence in object sequence.
	//
	size_t length = sequence.length();
	
	char c_sequence[length+1];
	std::string seqstring = MultipleAlignment(sequence).seqentry(0).seq().to_string();
	strcpy(c_sequence,seqstring.c_str());
	
	char c_structure[length+1];
	
	// ----------------------------------------
	// set folding parameters
	if (params.noLP) {noLonelyPairs=1;}
	
	
	// ----------------------------------------
	// call fold for setting the pf_scale
	double en = fold(c_sequence,c_structure);
	// std::cout << c_structure << std::endl;
	free_arrays();
	
	// set pf_scale
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	pf_scale = exp(-en/kT/length);

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
	McCmat=new McC_matrices_t(c_sequence,inLoopProbs); // makes local copy (if needed)
	
	// precompute further tables (expMLbase, scale, qm2) for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	if (inLoopProbs) {
	    expMLbase.resize(length+1);
	    scale.resize(length+1);
	    
	    // ----------------------------------------
	    // from scale_pf_params
	    //
	    kT = McCmat->pf_params->kT;   /* kT in cal/mol  */
	    
	    /* scaling factors (to avoid overflows) */
	    if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
		pf_scale = exp(-(-185+(McCmat->pf_params->temperature-37.)*7.27)/kT);
		if (pf_scale<1) pf_scale=1;
	    }
	    
	    scale[0] = 1.;
	    scale[1] = 1./pf_scale;
	    expMLbase[0] = 1;
	    expMLbase[1] = McCmat->pf_params->expMLbase * scale[1];
	    for (size_t i=2; i<=sequence.length(); i++) {
		scale[i] = scale[i/2]*scale[i-(i/2)]; // scale[i] = pow(scale[1],(double)i)
		expMLbase[i] = pow(McCmat->pf_params->expMLbase, (double)i) * scale[i];
	    }
	
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2();
	}
    }


    void
    RnaData::compute_McCaskill_alifold_matrices(const PFoldParams &params, bool inLoopProbs) {
	
	size_t length = sequence.length();
	size_t n_seq = sequence.row_number();
	
	// ----------------------------------------
	// write sequences to array of C-strings
	MultipleAlignment ma(sequence);
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
	double en = alifold(c_sequences,c_structure);
	// std::cout << c_structure << std::endl;
	free_alifold_arrays();
	
	// set pf_scale
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	pf_scale = exp(-en/kT/length);
	
	
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
	McCmat=new McC_ali_matrices_t(n_seq,length,inLoopProbs); // makes local copy (if needed)
	
	// precompute further tables (expMLbase, scale, qm2) for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	if (inLoopProbs) {
	    expMLbase.resize(length+1);
	    scale.resize(length+1);
	    
	    // ----------------------------------------
	    // from scale_pf_params
	    //
	    double scaling_factor=McCmat->pf_params->pf_scale;
	    kT = McCmat->pf_params->kT / n_seq;   /* kT in cal/mol  */
	    
	    
	    /* scaling factors (to avoid overflows) */
	    if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
		scaling_factor = exp(-(-185+(McCmat->pf_params->temperature-37.)*7.27)/kT);
		if (scaling_factor<1) scaling_factor=1;
		McCmat->pf_params->pf_scale=scaling_factor;
	    }
	    scale[0] = 1.;
	    scale[1] = 1./scaling_factor;

	    expMLbase[0] = 1;
	    expMLbase[1] = McCmat->pf_params->expMLbase/scaling_factor;
	    for (size_t i=2; i<=sequence.length(); i++) {
		scale[i] = scale[i/2]*scale[i-(i/2)];
		expMLbase[i] = pow(McCmat->pf_params->expMLbase, (double)i) * scale[i];
	    }
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2_ali();
	}

	//free c_sequences and c_structure
	for (size_t i=0; i<n_seq; i++) {
	    delete c_sequences[i];
	}
	delete c_sequences;
	delete c_structure;
    }
    
#endif // HAVE_LIBRNA

    void 
    RnaData::initFromFile(const std::string &filename,
			  bool readPairProbs,
			  bool readStackingProbs,
			  bool readInLoopProbs			       
			  ) {
	assert(!readStackingProbs || readPairProbs);
	assert(!readInLoopProbs || readPairProbs);

	std::ifstream in(filename.c_str());
	if (! in.good()) {
	    std::cerr << "Cannot read "<<filename<<std::endl;
	    exit(-1);
	}

	// set to true if probs become available
	in_loop_probs_available=false;
	pair_probs_available=false;
	stacking_probs_available=false;

	std::string s;
	// read first line and decide about file-format
	in >> s;
	in.close();
	if (s == "%!PS-Adobe-3.0") {
	    // try reading as dot.ps file (as generated by RNAfold)
	    readPS(filename,readPairProbs,readStackingProbs);
	} else if (s.substr(0,7) == "CLUSTAL" || s[0]=='>') {
	    //read to multiple alignment object
	    MultipleAlignment ma(filename,
				 (s[0]=='>')
				 ?MultipleAlignment::FASTA
				 :MultipleAlignment::CLUSTAL);
	    // convert to sequence
	    sequence = Sequence(ma);
	} else {
	    // try reading as PP-file (proprietary format, which is easy to read and contains pair probs)
	    readPP(filename,readPairProbs,readStackingProbs,readInLoopProbs);
	}
	
	// DUMP for debugging
	//std::cout << arc_probs_ << std::endl;
	//std::cout << arc_2_probs_ << std::endl;
    }

    void RnaData::readPS(const std::string &filename,
			 bool readPairProbs,
			 bool readStackingProbs
			 ) {
	assert(!readStackingProbs || readPairProbs);
		
	in_loop_probs_available=false;
	pair_probs_available=readPairProbs;
	
	bool contains_stacking_probs=false;
	
	std::ifstream in(filename.c_str());
    	
	std::string s;
	while (in >> s && s!="/sequence") {
	    if (s=="stacked") contains_stacking_probs=true;
	}

	stacking_probs_available = readStackingProbs && contains_stacking_probs;
	
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
        
	sequence.append_row(seqname,seqstr);
            
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
		
		    if (! (1<=i && i<j && j<=sequence.length())) {
			std::cerr << "WARNING: Input dotplot "<<filename<<" contains invalid line " << line << " (indices out of range)" << std::endl;
			//exit(-1);
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


    void RnaData::readPP(const std::string &filename,
			 bool readPairProbs,
			 bool readStackingProbs,
			 bool readInLoopProbs
			 ) {
	
	assert(!readStackingProbs || readPairProbs);
	assert(!readInLoopProbs || readPairProbs);
	
	in_loop_probs_available=false; // this is not (yet) supported!
	pair_probs_available=readPairProbs;
	
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
	    sequence.append_row(it->first,it->second);
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
		std::cerr << "Error in PP input line \""<<line<<"\" (i>=j).\n"<<std::endl;
		exit(-1);
	    }
      
	    set_arc_prob(i,j,p);
      
	    if (readStackingProbs) {
		double p2;
	  
		if (in >> p2) {
		    set_arc_2_prob(i,j,p2); // p2 is joint prob of (i,j) and (i+1,j-1)
		    stacking_probs_available = true;
		}
	    }
	}
    }
    
    void RnaData::clear_arc_probs() {
	arc_probs_.clear();
	arc_2_probs_.clear();
	pair_probs_available=false;
	stacking_probs_available=false;
    }


#ifdef HAVE_LIBRNA
    void
    RnaData::set_arc_probs_from_McCaskill_bppm(double threshold, bool stacking) {
	clear_arc_probs();
	
	for( size_t i=1; i <= sequence.length(); i++ ) {
	    for( size_t j=i+1; j <= sequence.length(); j++ ) {
		
		double p= McCmat->get_bppm(i,j);
		
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
	
	pair_probs_available=true;
	stacking_probs_available=stacking;
    }

    void
    RnaData::compute_Qm2(){
	assert(!used_alifold);
	McC_matrices_t *McCmat = static_cast<McC_matrices_t *>(this->McCmat);

	size_type len = sequence.length();
	

	std::vector<FLT_OR_DBL> qqm(len+2,0);
	std::vector<FLT_OR_DBL> qqm1(len+2,0);
	
	//qm1.resize((len+1)*(len+2)/2);
	qm2.resize((len+1)*(len+2)/2);
	
	// initialize qqm1
	for (size_type i=1; i<=len; i++) {
	    qqm[i]=0;
	    qqm1[i]=0;
	}
	
	for(size_type j=TURN+2; j<=len; j++) {
	    // --------------------
	    // one column of Qm1, which will be needed in the calculation of Qm2  
	    for(size_type i=j-TURN-1; i>=1; i--) {
		char type=McCmat->get_ptype(i,j);
		qqm[i]= qqm1[i]*expMLbase[1];
		if(type) {
		    qqm[i] +=
			McCmat->get_qb(i,j)
			* exp_E_MLstem(type,
				       (i>1)   ? McCmat->S1[i-1] : -1, 
				       (j<len) ? McCmat->S1[j+1] : -1,  
				       McCmat->pf_params);
		}
		
		//qm1[McCmat->iidx(i,j)]=qqm[i];

		assert(qqm[i] <= McCmat->get_qm(i,j));
		assert((!frag_len_geq(i,j-1,TURN+2)) || qqm1[i] <= McCmat->get_qm(i,j-1));
	    }
	    	    
	    // --------------------
	    // calculates column "j" of the Qm2 matrix
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*(TURN+2)+1; i>=1; i--) {
		    qm2[McCmat->iidx(i,j)] = 0;
		    for(size_type k = i + TURN+1; (k+1)+TURN+1 <= j; k++) {
			qm2[McCmat->iidx(i,j)] +=
			    McCmat->get_qm(i,k)*qqm[k+1];
		    }
		    assert(qm2[McCmat->iidx(i,j)] <= McCmat->get_qm(i,j));
		}
	    }
	    	    
	    // --------------------
	    // swap qqm and qqm1 (in constant time)
	    qqm1.swap(qqm);
	}
    }


    void
    RnaData::compute_Qm2_ali(){
	assert(used_alifold);
	assert(McCmat);
	McC_ali_matrices_t *McCmat = static_cast<McC_ali_matrices_t *>(this->McCmat);

	size_type len   = sequence.length();
	size_type n_seq = sequence.row_number();
	
	std::vector<FLT_OR_DBL> qqm(len+2,0);
	std::vector<FLT_OR_DBL> qqm1(len+2,0);
	std::vector<int> type(n_seq);
	
	qm2.resize((len+1)*(len+2)/2);
	
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
		    type[s] = pair[McCmat->S[s][i]][McCmat->S[s][j]];
		    if (type[s]==0) type[s]=7;
		}
		
		qqm[i]= qqm1[i]*expMLbase[1];
		
		FLT_OR_DBL qbt1=1.0; // collects contribution "inner basepair of multiloop"
		for (size_t s=0; s<n_seq; s++) {
		    qbt1 *= exp_E_MLstem(type[s], 
					 i>1 ? McCmat->S5[s][i] : -1,
					 j<len ? McCmat->S3[s][j] : -1,
					 McCmat->pf_params);
		}
		qqm[i] += McCmat->get_qb(i,j) * qbt1;

	    }
	    
	    // --------------------
	    // calculate a row of the matrix Qm2
	    //
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*TURN-3; i>=1; i--) {
		    qm2[McCmat->iidx(i+1,j-1)] = 0;
		    for(size_type k = i+TURN+2; k< j-TURN-2; k++) {
			qm2[McCmat->iidx(i+1,j-1)] +=
			    McCmat->get_qm(i+1,k)*qqm1[k+1];
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
	assert(!used_alifold);
	McC_matrices_t *McCmat = static_cast<McC_matrices_t *>(this->McCmat);

    	int type = McCmat->get_ptype(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if ((type==0)
	    || (((type==3)||(type==4))&&no_closingGU)
	    || (McCmat->get_qb(i,j)==0.0)
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
	
	if (!in_loop_probs_available) return 1.0;
	
	McC_ali_matrices_t *McCmat = static_cast<McC_ali_matrices_t*>( this->McCmat );
	
        size_t n_seq = sequence.row_number();
	
	// immediately return 0.0 if i and j do not pair
	if (get_arc_prob(i,j)==0.0 || McCmat->get_qb(i,j)==0.0) {return 0.0;}
	
	// get base pair types for i,j of all sequences
	std::vector<int> type(n_seq);
	
	for (size_t s=0; s<n_seq; ++s) {
	    type[s] = pair[McCmat->S[s][i]][McCmat->S[s][j]];
	    if (type[s]==0) type[s]=7;
	}

	// ------------------------------------------------------------
	// hairpin contribution
        //
	
	FLT_OR_DBL H=1.0;
	
	for (size_t s=0; s<n_seq; s++) {
	    size_t u = McCmat->a2s[s][j-1]-McCmat->a2s[s][i];
	    if (McCmat->a2s[s][i]<1) continue;
	    char loopseq[10];
	    if (u<7){
		strncpy(loopseq, McCmat->Ss[s]+McCmat->a2s[s][i]-1, 10);
	    }
	    H *= exp_E_Hairpin(u, type[s],
			       McCmat->S3[s][i], McCmat->S5[s][j], 
			       loopseq, 
			       McCmat->pf_params);
        }
        H *= scale[j-i+1];
    
	// ------------------------------------------------------------
	// interior loop contributions
	//
	
	FLT_OR_DBL I = 0.0;
	
	// case 1: i<k<i´<j´<j
	for (size_t ip=k+1; ip <= std::min(i+MAXLOOP+1,j-TURN-2); ip++) {
	    for (size_t jp = std::max(ip+TURN+1 + MAXLOOP ,j-1 + ip-i-1 ) - MAXLOOP; jp<j; jp++) {
		
		FLT_OR_DBL qloop=1.0;
		
		if (McCmat->get_qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = McCmat->a2s[s][ip-1] - McCmat->a2s[s][i];
		    size_t u2 = McCmat->a2s[s][j-1] - McCmat->a2s[s][jp];
		    
		    int type_2 = pair[McCmat->S[s][jp]][McCmat->S[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    McCmat->S3[s][i],
					    McCmat->S5[s][j],
					    McCmat->S5[s][ip],
					    McCmat->S3[s][jp],
					    McCmat->pf_params
					    );
		}
		
		I += McCmat->get_qb(ip,jp) * scale[ip-i+j-jp] * qloop;
	    }
	}
 
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip <= std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    for (size_t jp = std::max(ip+TURN+1 + MAXLOOP,j-1+ ip-i-1) - MAXLOOP; jp<k; jp++) {
		
		FLT_OR_DBL qloop=1.0;
		
		if (McCmat->get_qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = McCmat->a2s[s][ip-1] - McCmat->a2s[s][i];
		    size_t u2 = McCmat->a2s[s][j-1] - McCmat->a2s[s][jp];
		    
		    int type_2 = pair[McCmat->S[s][jp]][McCmat->S[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    McCmat->S3[s][i],
					    McCmat->S5[s][j],
					    McCmat->S5[s][ip],
					    McCmat->S3[s][jp],
					    McCmat->pf_params
					    );
		}
		
		I += McCmat->get_qb(ip,jp) * scale[ip-i+j-jp] * qloop;
	    }
	}
    
    
	// ------------------------------------------------------------
	// multiloop contributions
        //

	FLT_OR_DBL M = 0.0;
	
	// no base pair <= k:   i....k-----qm2-------j
	// valid entries of qm2 have space for 2 inner base pairs,
	// i.e. at least length of "(...)(...)" (for TURN=3)
	if ( frag_len_geq(k+1, j-1, 2*(TURN+2)) ) {
	    M += qm2[McCmat->iidx(k+1,j-1)] * expMLbase[k-i];
	}
	
	// no base pair >= k
	if ( frag_len_geq(i+1,k-1,2*(TURN+2)) ) {
	    M += qm2[McCmat->iidx(i+1,k-1)] * expMLbase[j-k];
	}
	
	// base pairs <k and >k
	if ( frag_len_geq(i+1,k-1,TURN+2) && frag_len_geq(k+1,j-1,TURN+2) ) {
	    M += McCmat->get_qm(i+1,k-1) * expMLbase[1] *  McCmat->get_qm(k+1,j-1);
	}
	
	// multiply with contribution for closing of multiloop
	
	for (size_t s=0; s<n_seq; s++) {
	    int tt = rtype[type[s]];
	    
	    M *= McCmat->pf_params->expMLclosing 
		* exp_E_MLstem(tt,McCmat->S5[s][j],McCmat->S3[s][i], McCmat->pf_params);
	}
	M *= scale[2];
	
	FLT_OR_DBL Qtotal=H+I+M;
	
	double kTn   = McCmat->pf_params->kT/10.;   /* kT in cal/mol  */
	
	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	Qtotal *= exp(McCmat->get_pscore(i,j)/kTn);
	
	FLT_OR_DBL p_k_cond_ij = Qtotal/McCmat->get_qb(i,j); 
	
	FLT_OR_DBL res = p_k_cond_ij * McCmat->get_bppm(i,j);
	
	return res;
    }
    
    double RnaData::prob_unpaired_in_loop(size_type k,size_type i,size_type j) const {
	assert(i+TURN+1 <= j);
	assert(i<k);
	assert(k<j);
	
	if (!in_loop_probs_available) return 1.0;
	
	
	if (used_alifold) {
	    return prob_unpaired_in_loop_ali(k, i, j);
	}
	
	assert(!used_alifold);
	McC_matrices_t *McCmat = static_cast<McC_matrices_t *>(this->McCmat);

	const char *c_sequence = McCmat->sequence;
	
	
	int type = ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}

	// ------------------------------------------------------------
	// Hairpin loop energy contribution
	
	size_t u=j-i-1;
	FLT_OR_DBL H = exp_E_Hairpin(u, type, McCmat->S1[i+1], McCmat->S1[j-1],
				     c_sequence+i-1, McCmat->pf_params) * scale[u+2];
	
	// ------------------------------------------------------------
	// Interior loop energy contribution
	FLT_OR_DBL I = 0.0;

	// case 1: i<k<i´<j´<j
	for (size_t ip=k+1; ip <= std::min(i+MAXLOOP+1,j-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp = std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<j; jp++) {
		int type2 = McCmat->get_ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += McCmat->get_qb(ip,jp) 
			* (scale[u1+j-jp+1] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 McCmat->S1[i+1],McCmat->S1[j-1],
					 McCmat->S1[ip-1],McCmat->S1[jp+1], McCmat->pf_params));
		}
	    }
	}
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip<=std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp=std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<k; jp++) {
		int type2 = McCmat->get_ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += McCmat->get_qb(ip,jp)
			* (scale[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 McCmat->S1[i+1],McCmat->S1[j-1],
					 McCmat->S1[ip-1],McCmat->S1[jp+1], McCmat->pf_params));
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
	    M1 = expMLbase[frag_len(i+1,k)] * qm2[McCmat->iidx(k+1,j-1)];
	}
	
	// bases >=k unpaired
	if ( frag_len_geq(i+1, k-1, 2*(TURN+2)) ) {
	    M2 = qm2[McCmat->iidx(i+1,k-1)] * expMLbase[frag_len(k,j-1)];
	}
	
	// innner base pairs left and right of k
	if ( frag_len_geq(i+1,k-1,TURN+2) && frag_len_geq(k+1,j-1,TURN+2) ) {
	    M3 = McCmat->get_qm(i+1,k-1) * expMLbase[1] *  McCmat->get_qm(k+1,j-1);
	}
	
	M=M1+M2+M3;
	
	// multiply with contribution for closing of multiloop
	M *= McCmat->pf_params->expMLclosing 
	    * exp_E_MLstem(rtype[type],McCmat->S1[j-1],McCmat->S1[i+1], McCmat->pf_params)
	    * scale[2];
	
	FLT_OR_DBL Qtotal = H+I+M;

	FLT_OR_DBL p_k_cond_ij = Qtotal/McCmat->get_qb(i,j); 

	FLT_OR_DBL res = p_k_cond_ij * McCmat->get_bppm(i,j);
	
	return res;
    }

    double RnaData::prob_unpaired_external(size_type k) const {
	assert(1<=k);
	assert(k<=sequence.length());
	
	if (!in_loop_probs_available) return 1.0;
	
	return (McCmat->q1k[k-1] * scale[1] * McCmat->qln[k+1]) / McCmat->qln[1];
    }


    double
    RnaData::prob_basepair_in_loop_ali(size_type ip,
				       size_type jp,
				       size_type i,
				       size_type j) const {

	if (!in_loop_probs_available) return 1.0;
	
	McC_ali_matrices_t *McCmat = static_cast<McC_ali_matrices_t *>(this->McCmat);
	
	size_t n_seq = sequence.row_number();
	
	// note: the following tests cover the case that the distances of i,j and ip,jp are too small

	// immediately return 0.0 if i and j do not pair
	if (get_arc_prob(i,j)==0.0 || McCmat->get_qb(i,j)==0.0) {return 0.0;}
	
	// immediately return 0.0 when ip and jp cannot pair
	if (get_arc_prob(ip,jp)==0.0 || McCmat->get_qb(ip,jp)==0.0) {return 0.0;}
	
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
	    type[s] = pair[McCmat->S[s][i]][McCmat->S[s][j]];
	    if (type[s]==0) type[s]=7;

	    type2[s] = pair[McCmat->S[s][ip]][McCmat->S[s][jp]];
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
		
		size_t u1 = McCmat->a2s[s][ip-1] - McCmat->a2s[s][i];
		size_t u2 = McCmat->a2s[s][j-1] - McCmat->a2s[s][jp];
		
		I *= exp_E_IntLoop(u1,u2,
				   type[s], rtype[type2[s]],
				   McCmat->S3[s][i],
				   McCmat->S5[s][j],
				   McCmat->S5[s][ip],
				   McCmat->S3[s][jp],
				   McCmat->pf_params);
	    }
	    I *= scale[ip-i+j-jp];
	}
	
	// ------------------------------------------------------------
	// Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase[frag_len(i+1,ip-1)] * McCmat->get_qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += McCmat->get_qm(i+1,ip-1) * expMLbase[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += McCmat->get_qm(i+1,ip-1) * McCmat->get_qm(jp+1,j-1);
	}
	
	for (size_t s=0; s<n_seq; s++) {
	    // multiply with factor for inner base pair
	    M *= exp_E_MLstem(type2[s],
			      McCmat->S5[s][ip],
			      McCmat->S3[s][jp],
			      McCmat->pf_params);
	    
	    // multiply with factors for closing base pair
	    M *= McCmat->pf_params->expMLclosing
		* exp_E_MLstem(rtype[type[s]],
			       McCmat->S5[s][j],
			       McCmat->S3[s][i],
			       McCmat->pf_params);
	}
	
	M *= scale[2]; // scale for closing base pair
	
	// ------------------------------------------------------------
	FLT_OR_DBL Qtotal= (I+M);
	
	Qtotal *= McCmat->get_qb(ip,jp);

	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	double kTn   = McCmat->pf_params->kT/10.;   /* kT in cal/mol  */
	Qtotal *= exp(McCmat->get_pscore(i,j)/kTn);

	
	return
	    (Qtotal/McCmat->get_qb(i,j))
	    *McCmat->get_bppm(i,j);
    }

    double
    RnaData::prob_basepair_in_loop(size_type ip,
				   size_type jp,
				   size_type i,
				   size_type j) const {
	
	if (!in_loop_probs_available) return 1.0;

	if (used_alifold) {
	    return prob_basepair_in_loop_ali(ip, jp, i, j);
	}
	
	assert(!used_alifold);
	McC_matrices_t *McCmat = static_cast<McC_matrices_t *>(this->McCmat);
	

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
			      McCmat->S1[(int)(i+1)],
			      McCmat->S1[(int)(j-1)],
			      McCmat->S1[ip-1],
			      McCmat->S1[jp+1],
			      McCmat->pf_params)
		* scale[u1+u2+2];
	}
	
	//calculating Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	
	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase[frag_len(i+1,ip-1)] * McCmat->get_qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += McCmat->get_qm(i+1,ip-1) * expMLbase[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += McCmat->get_qm(i+1,ip-1) * McCmat->get_qm(jp+1,j-1);
	}
	
	// multiply with factor for inner base pair
	M *= exp_E_MLstem(type2,
			  McCmat->S1[ip-1],
			  McCmat->S1[jp+1],
			  McCmat->pf_params);
	
	// multiply with factors for closing base pair
	M *= McCmat->pf_params->expMLclosing
	    * exp_E_MLstem(rtype[type],
			   McCmat->S1[j-1],
			   McCmat->S1[i+1],
			   McCmat->pf_params)
	    * scale[2];
	
	FLT_OR_DBL Qtotal = I+M;
	
	Qtotal *= McCmat->get_qb(ip,jp);
	
	
	return Qtotal/McCmat->get_qb(i,j)
	    * McCmat->get_bppm(i,j);
    }

    double RnaData::prob_basepair_external(size_type i,size_type j) const {

	if (!in_loop_probs_available) return 1.0;
	
	size_t n=sequence.length();

	assert(1<=i);
	assert(i<j);
	assert(j<=n);
	assert(frag_len_geq(i,j,TURN+2));
	
	// immediately return 0.0 when i and j cannot pair
	if (get_arc_prob(i,j)==0.0 || McCmat->get_qb(i,j)==0.0) {
	    return 0.0;
	}
	
	FLT_OR_DBL extloop;
	
	if (!used_alifold) {
	    McC_matrices_t *McCmat = static_cast<McC_matrices_t *>(this->McCmat);
	    extloop = exp_E_ExtLoop(McCmat->get_ptype(i,j),
				    i>1 ? McCmat->S1[i-1] : -1, 
				    j<n ? McCmat->S1[j+1] : -1, 
				    McCmat->pf_params);
	} else {
	    McC_ali_matrices_t *McCmat = static_cast<McC_ali_matrices_t *>(this->McCmat);
	    
	    size_t n_seq=sequence.row_number();
	    
	    extloop=1.0;
	    
	    for (size_t s=0; s<n_seq; s++) {
		int type = pair[McCmat->S[s][i]][McCmat->S[s][j]];
		if (type==0) type=7;
		
		extloop *= exp_E_ExtLoop(type, i>1 ? McCmat->S5[s][i] : -1,
					 j<n ? McCmat->S3[s][j] : -1,
					 McCmat->pf_params);
	    }
	}

	return
	    (McCmat->q1k[i-1]
	     * McCmat->get_qb(i,j)
	     * extloop
	     * McCmat->qln[j+1]
	     )
	    / McCmat->qln[1];

    }
    
    std::ostream &
    RnaData::write_unpaired_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
	
	// write lines for loops closed by base pairs
	for(arc_prob_matrix_t::size_type i=1; i<=sequence.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence.length(); ++j) {
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
	for(arc_prob_matrix_t::size_type k=1; k<=sequence.length(); ++k) {
	    double p=prob_unpaired_external(k);
	    if (p>threshold2) {
		if (!had_entries) {out << 0 << " " << (sequence.length()+1); had_entries=true;}
		out << " " << k << " " << p;
	    }
	}
	if (had_entries) out << std::endl;
    }
    
	
    std::ostream &
    RnaData::write_basepair_in_loop_probs(std::ostream &out,double threshold1,double threshold2) const {
	// write lines for loops closed by base pairs
	for(arc_prob_matrix_t::size_type i=1; i<=sequence.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence.length(); ++j) {
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
	for(arc_prob_matrix_t::size_type ip=1; ip<=sequence.length(); ++ip) {
	    for(arc_prob_matrix_t::size_type jp=ip+1; jp<=sequence.length(); ++jp) {
		double p=prob_basepair_external(ip,jp);
		if (p>threshold2) {
		    if (!had_entries) {out << 0 << " " << (sequence.length()+1); had_entries=true;}
		    out << " " << ip << " " << jp << " " << p;
		}
	    }
	}
	if (had_entries) out << std::endl;
    }


    std::ostream &
    RnaData::write_basepair_and_in_loop_probs(std::ostream &out,double threshold1,double threshold2,double threshold3) const {
	for(arc_prob_matrix_t::size_type i=1; i<=sequence.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence.length(); ++j) {
		
		if (arc_probs_(i,j)>threshold1) {
		    out << i << " " << j;
		    
		    // write base pair and stacking probability
		    out << " " << arc_probs_(i,j);
		    if (arc_2_probs_(i,j)>threshold1) {
			out << " " << arc_2_probs_(i,j);
		    }
		    out << " ;";
		    
		    // write unpaired in loop
		    for(arc_prob_matrix_t::size_type k=i+1; k<=j-1; ++k) {
			double p=prob_unpaired_in_loop(k,i,j);
			if (p>threshold2) {
			    out << " " << k << " " << p;
			}
		    }
		    
		    out << " ;";
		    
		    for(arc_prob_matrix_t::size_type ip=i+1; ip<=j-1; ++ip) {
			for(arc_prob_matrix_t::size_type jp=ip+1; jp<=j-1; ++jp) {
			    double p=prob_basepair_in_loop(ip,jp,i,j);
			    if (p>threshold3) {
				out << " " << ip << " " << jp << " " << p;
			    }
			}
		    }
		    out << std::endl;
		}
	    }
	}

	// write lines for external loop

	out << 0 << " " << (sequence.length()+1) << "1 1 ;";
	for(arc_prob_matrix_t::size_type k=1; k<=sequence.length(); ++k) {
	    double p=prob_unpaired_external(k);
	    if (p>threshold2) {
		out << " " << k << " " << p;
	    }
	}
	out << ";";
	
	for(arc_prob_matrix_t::size_type ip=1; ip<=sequence.length(); ++ip) {
	    for(arc_prob_matrix_t::size_type jp=ip+1; jp<=sequence.length(); ++jp) {
		double p=prob_basepair_external(ip,jp);
		if (p>threshold3) {
		    out << " " << ip << " " << jp << " " << p;
		}
	    }
	}
	out << std::endl;
	
    }


#endif // HAVE_LIBRNA

    std::ostream &
    RnaData::write_basepair_probs(std::ostream &out,double threshold) const {
	for(arc_prob_matrix_t::size_type i=1; i<=sequence.length(); ++i) {
	    for(arc_prob_matrix_t::size_type j=i+1; j<=sequence.length(); ++j) {
		
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

	if (name.length()>3 && name.substr(name.length()-3,3) == "_dp")
	    name=name.substr(0,name.length()-3);
    
	return name;
    }

}

