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
}
#endif // HAVE_LIBRNA


// // for getrusage()
// #include <sys/resource.h>
// #include <sys/types.h>
// // for gettimeofday()
// #include <sys/time.h>

namespace LocARNA {

    // ------------------------------------------------------------
    // implementation of class McC_matrices_t
#ifdef HAVE_LIBRNA	
    McC_matrices_t::McC_matrices_t(size_t length_, bool local_copy_): length(length_),local_copy(local_copy_) {
	    
	if (local_copy_) {
	    McC_matrices_t McCmat_tmp(length,false);
	    deep_copy(McCmat_tmp);
	} else {
	    
	    // get pointers to McCaskill matrices
	    get_pf_arrays(&S_p,
			  &S1_p,
			  &ptype_p,
			  &qb_p,
			  &qm_p,
			  &q1k_p,
			  &qln_p);
	    // get pointer to McCaskill base pair probabilities
	    bppm = export_bppm();
	    
	    iindx = get_iindx(length);
	    
	    pf_params_p = get_scaled_pf_parameters();
	}
    }
    
    void *
    space_memcpy(void *from,size_t size) {
	void *p = space(size);
	memcpy(p,from,size);
	return p;
    }

    void
    McC_matrices_t::deep_copy(const McC_matrices_t &McCmat) {
	length=McCmat.length;

	size_t size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
	
	S_p = (short *) space_memcpy(McCmat.S_p,sizeof(short)*(length+2));
	S1_p = (short *) space_memcpy(McCmat.S1_p,sizeof(short)*(length+2));
	ptype_p= (char *) space_memcpy(McCmat.ptype_p,sizeof(char)*((length+1)*(length+2)/2));
	qb_p= (FLT_OR_DBL *) space_memcpy(McCmat.qb_p,size);
	qm_p= (FLT_OR_DBL *) space_memcpy(McCmat.qm_p,size);
	bppm= (FLT_OR_DBL *) space_memcpy(McCmat.bppm,size);
	q1k_p= (FLT_OR_DBL *) space_memcpy(McCmat.q1k_p,sizeof(FLT_OR_DBL)*(length+1));
	qln_p= (FLT_OR_DBL *) space_memcpy(McCmat.qln_p,sizeof(FLT_OR_DBL)*(length+2));
	pf_params_p= (pf_paramT *) space_memcpy(McCmat.pf_params_p,sizeof(pf_paramT));
		
	iindx= get_iindx(length);
    }

    McC_matrices_t::~McC_matrices_t() {
	    if (local_copy) {
		free_all();
	    } else {
		free(iindx);
	    }
	}

    void McC_matrices_t::free_all() {
	free(S_p);
	free(S1_p);
	free(ptype_p);
	free(qb_p);
	free(qm_p);
	free(q1k_p);
	free(qln_p);
	free(bppm);
	free(iindx);
	free(pf_params_p);
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
	seq_constraints_("")
    {
	initFromFile(file, readPairProbs, readStackingProbs, readInLoopProbs);
	
	consensus_sequence = MultipleAlignment(sequence).consensus_sequence();	
    }
        
    RnaData::RnaData(const Sequence &sequence_)
	: sequence(sequence_),
	  pair_probs_available(false),	  
	  consensus_sequence(MultipleAlignment(sequence).consensus_sequence()),
	  arc_probs_(0),
	  arc_2_probs_(0),
	  seq_constraints_("")
    {
    }
    
    RnaData::~RnaData() {
#ifdef HAVE_LIBRNA
	if (McCmat) {delete McCmat;} 
#endif //HAVE_LIBRNA
    }

#ifdef HAVE_LIBRNA

    void
    RnaData::computeEnsembleProbs(const PFoldParams &params,bool inLoopProbs) {
	// run McCaskill and get access to results
	// in McCaskill_matrices
	compute_McCaskill_matrices(params,inLoopProbs);
	
	// initialize the object from base pair probabilities
	// Use the same proability threshold as in RNAfold -p !
	set_arc_probs_from_McCaskill_bppm(10e-6,params.stacking);	
	
	// since we either have local copies of all McCaskill pf arrays
	// or don't need them anymore,
	// we can free the ones of the Vienna lib
	free_pf_arrays();
    }
    
    void
    RnaData::compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs) {	
	if (sequence.row_number()!=1) {
	    std::cerr << "McCaskill computation with multi-row Sequence object is not implemented." << std::endl;
	    exit(-1);
	}
	
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
	McCmat=new McC_matrices_t(length,inLoopProbs); // here is the local copy (if needed)
	
	// precompute further tables (expMLbase, scale, qm2) for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	if (inLoopProbs) {
	    expMLbase.resize(length+1);
	    scale.resize(length+1);
	    
	    kT = McCmat->pf_params_p->kT;   /* kT in cal/mol  */
	    
	    /* scaling factors (to avoid overflows) */
	    if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
		pf_scale = exp(-(-185+(McCmat->pf_params_p->temperature-37.)*7.27)/kT);
		if (pf_scale<1) pf_scale=1;
	    }
	    scale[0] = 1.;
	    scale[1] = 1./pf_scale;
	    expMLbase[0] = 1;
	    expMLbase[1] = McCmat->pf_params_p->expMLbase/pf_scale;
	    for (size_t i=2; i<=sequence.length(); i++) {
		scale[i] = scale[i/2]*scale[i-(i/2)];
		expMLbase[i] = pow(McCmat->pf_params_p->expMLbase, (double)i) * scale[i];
	    }
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2();
	}
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
	size_type len = sequence.length();

	std::vector<FLT_OR_DBL> qm1(len+2,0);
	std::vector<FLT_OR_DBL> qqm1(len+2,0);
	
	qm2.resize((len+1)*(len+2)/2);
	
      
	for(size_type index_j= TURN+2; index_j<=len; index_j++) {
	    // --------------------
	    // the first inner loop calculates one row of Qm1 that will be needed in the calculation of Qm2  
	    for(size_type index_i=index_j-TURN-1; index_i>=1; index_i--) {
		char type=McCmat->get_ptype(index_i,index_j);
		qm1[index_i]= qqm1[index_i]*expMLbase[1];
		if(type){
		    qm1[index_i] +=
			(McCmat->get_qb(index_i,index_j))
			* exp_E_MLstem(type,
				       (index_i>1) ? McCmat->S1_p[index_i-1] : -1, 
				       (index_j<len) ? McCmat->S1_p[index_j+1] : -1,  
				       McCmat->pf_params_p);
		}
	    }
	    // --------------------
	    // this part calculates the Qm2 matrix
	    if(index_j >= (2*(TURN+2))) {
		for(size_type index_i = index_j-2*TURN-3; index_i>=1; index_i--) {
		    qm2[McCmat->idx(index_i,index_j)] = 0;
		    for(size_type index_k = index_i+2; index_k< index_j-2; index_k++) {
			qm2[McCmat->idx(index_i+1,index_j-1)] +=
			    McCmat->get_qm(index_i+1,index_k)*qqm1[index_k+1];
		    }
		}
	    }
	    // --------------------
	    // swap row qm1 and qqm1 (in constant time)
	    qqm1.swap(qm1);
	}
    }

    
    char
    RnaData::ptype_of_admissible_basepair(size_type i,size_type j) const {
    	char type = McCmat->get_ptype(i,j);
	
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

    double RnaData::prob_unpaired_in_loop(size_type k,size_type i,size_type j) const {
	
	const char *c_sequence=consensus_sequence.c_str(); //!< @todo: check! This looks wrong!
	
	FLT_OR_DBL H,I,M;
	//calculating the Hairpin loop energy contribution
	
	char type = ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}

	H = exp_E_Hairpin((int)(j-i-1), type, McCmat->S1_p[i+1], McCmat->S1_p[j-1],
			  c_sequence+i-1, McCmat->pf_params_p) * scale[j-i+1];
	
	I = 0.0;
	//calculating the Interior loop energy contribution

	int u1;
	// case 1: i<k<i´<j´<j
	for (int ip=k+1; 
	     ip<=MIN2((int)(i+MAXLOOP+1),(int)(j-TURN-2));
	     ip++) {
	    u1 = (int)(ip-i-1);
	    for (int jp = MAX2((int)(ip+TURN+1),(int)(j-1-MAXLOOP+u1)); 
		 jp<(int)j;
		 jp++) {
		char type2 = McCmat->get_ptype(ip,jp); //!< @todo check! admissible type?
		if (type2) {
		    type2 = rtype[type2];
		    I += McCmat->get_qb(ip,jp) 
			* (scale[u1+j-jp+1] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 McCmat->S1_p[i+1],McCmat->S1_p[j-1],
					 McCmat->S1_p[ip-1],McCmat->S1_p[jp+1], McCmat->pf_params_p));
		}
	    }
	}
	//case 2: i<i´<j´<k<j
	for (int ip=i+1;
	     ip<=MIN2((int)(i+MAXLOOP+1),(int)(k-TURN-2)); 
	     ip++) {
	    u1 = ip-i-1;
	    for (int jp=MAX2((int)(ip+TURN+1),(int)(j-1-MAXLOOP+u1));
		 jp<(int)k;
		 jp++) {
		char type2 = McCmat->get_ptype(ip,jp); //!< @todo check! admissible type?
		if (type2) {
		    type2 = rtype[type2];
		    I += McCmat->get_qb(ip,jp)
			* (scale[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 McCmat->S1_p[i+1],McCmat->S1_p[j-1],
					 McCmat->S1_p[ip-1],McCmat->S1_p[jp+1], McCmat->pf_params_p));
		}
	    }
	}
	
	//calculating Multiple loop energy contribution
	M = 0.0;
	
	M += qm2[McCmat->idx(k+1,j-1)] * expMLbase[k-i];

	M += qm2[McCmat->idx(i+1,k-1)] * expMLbase[j-k];
	    
	M += McCmat->get_qm(i+1,k-1) * expMLbase[1] *  McCmat->get_qm(k+1,j-1);
	    
	// multiply with contribution for closing of multiloop
	M *= McCmat->pf_params_p->expMLclosing 
	    * exp_E_MLstem(rtype[type],McCmat->S1_p[j-1],McCmat->S1_p[i+1], McCmat->pf_params_p)
	    * scale[2];
	
	return ((H+I+M)/McCmat->get_qb(i,j))*get_arc_prob(i,j);
    }

    double RnaData::prob_unpaired_external(size_type k) const {
	return (McCmat->q1k_p[k-1] * scale[1] * McCmat->qln_p[k+1]) / McCmat->qln_p[1];
    }

    double
    RnaData::prob_basepair_in_loop(size_type ip,
				   size_type jp,
				   size_type i,
				   size_type j) const {
	FLT_OR_DBL Ipp;
	FLT_OR_DBL Mpp;
	
	// note: Ipp and Mpp are computed without factor get_qb(ip,jp),
	// which is multiplied only in the end.
	
	int type=ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}
	
	int type2=ptype_of_admissible_basepair(ip,jp);
	
	// immediately return 0.0 when ip and jp cannot pair
	if (type2==0) {return 0.0;}
	
	//calculating the Interior loop energy contribution
	//
	Ipp=0.0;
	
	int u1 =(int)(ip-i-1);
	int u2 =(int)(j-jp-1);
	
	Ipp = exp_E_IntLoop(u1,u2, type, rtype[type2],
			    McCmat->S1_p[(int)(i+1)],
			    McCmat->S1_p[(int)(j-1)],
			    McCmat->S1_p[ip-1],
			    McCmat->S1_p[jp+1],
			    McCmat->pf_params_p)
	    * scale[u1+u2+2];
	
	//calculating Multiple loop energy contribution
	//
	Mpp = 0.0;

	// inner base pairs only right of (ip,jp)
	Mpp += expMLbase[ip-i-1] * McCmat->get_qm(jp+1,j-1);
	
	// inner base pairs only left of (ip,jp)
	Mpp += McCmat->get_qm(i+1,ip-1) * expMLbase[j-jp-1];
	
	// inner base pairs left and right of (ip,jp)
	Mpp += McCmat->get_qm(i+1,ip-1) * McCmat->get_qm(jp+1,j-1);
	
	// multiply with factor for inner base pair
	Mpp *= exp_E_MLstem(type2,
			    McCmat->S1_p[ip-1],
			    McCmat->S1_p[jp+1],
			    McCmat->pf_params_p);
	
	// multiply with factors for closing base pair
	Mpp *= McCmat->pf_params_p->expMLclosing
	    * exp_E_MLstem(rtype[type],
			   McCmat->S1_p[j-1],
			   McCmat->S1_p[i+1],
			   McCmat->pf_params_p)
	    * scale[2];
	
	return 
	    (McCmat->get_qb(ip,jp)
	     *(Ipp+Mpp)/McCmat->get_qb(i,j))
	    *get_arc_prob(i,j);
    }

    double RnaData::prob_basepair_external(size_type i,size_type j) const {
	// immediately return 0.0 when i and j cannot pair
	if (ptype_of_admissible_basepair(i,j)==0) {return 0.0;}
	
	return 
	    (McCmat->q1k_p[i-1]
	     * McCmat->get_qb(i,j)
	     * McCmat->qln_p[j+1]
	     )
	    / McCmat->qln_p[1];
    }

#endif // HAVE_LIBRNA
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
