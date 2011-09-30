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
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/loop_energies.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/pair_mat.h>
}
#endif // HAVE_LIBRNA

#include <math.h>
#include <string.h>

namespace LocARNA {

    RnaData::RnaData(const std::string &file, bool stacking_, bool keepMcM):
	sequence(),
	stacking(stacking_),
	arc_probs_(0),
	arc_2_probs_(0),
	seq_constraints_("")
    {
	read(file, keepMcM);

	consensus_sequence = MultipleAlignment(sequence).consensus_sequence();	
    }
    
#ifdef HAVE_LIBRNA

    
    RnaData::RnaData(const Sequence &sequence_, bool keepMcM, bool stacking_)
	: sequence(sequence_),
	  stacking(stacking_),
	  consensus_sequence(MultipleAlignment(sequence).consensus_sequence()),
	  arc_probs_(0),
	  arc_2_probs_(0),
	  seq_constraints_("")
    {
	if (sequence.row_number()!=1) {
	    std::cerr << "Construction with multi-row Sequence object is not implemented." << std::endl;
	    exit(-1);
	}
	
	// run McCaskill and get access to results
	// in McCaskill_matrices
	compute_McCaskill_matrices();
	
	// initialize the object from base pair probabilities
	// Use the same proability threshold as in RNAfold -p !
	init_from_McCaskill_bppm();	
	
	// since we have local copies of all McCaskill pf arrays,
	// we can free the ones of the Vienna lib
	free_pf_arrays();
	
	// optionally deallocate McCaskill matrices
	if (!keepMcM) {
	    free_McCaskill_matrices();
	}
    }
    
    void
    RnaData::compute_McCaskill_matrices() {	
	if (sequence.row_number()!=1) {
	    std::cerr << "McCaskill computation with multi-row Sequence object is not implemented." << std::endl;
	    exit(-1);
	}
	
	assert(sequence.row_number()==1);
	
	// use MultipleAlignment to get pointer to c-string of the
	// first (and only) sequence in object sequence.
	//
	size_t length = sequence.length();
	
	char c_sequence[length+1];
	std::string seqstring = MultipleAlignment(sequence).seqentry(0).seq().to_string();
	strcpy(c_sequence,seqstring.c_str());
	
	char c_structure[length+1];
	//p_sequence= c_sequence;
	
	// std::cout <<"Call fold(" << c_sequence << "," << "c_structure" << ")"<< std::endl;
	
	// call fold for setting the pf_scale
	
	double en = fold(c_sequence,c_structure);
	// std::cout << c_structure << std::endl;
	free_arrays();
	
	// set pf_scale
	double kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
	pf_scale = exp(-en/kT/length);
	
	// std::cout <<"Call pf_fold(" << c_sequence << "," << "NULL" << ")"<< std::endl;
	
	// call pf_fold
	pf_fold(c_sequence,c_structure);
	
	
	McC_matrices_t McCmat;

	// get pointers to McCaskill matrices
	get_pf_arrays(&McCmat.S_p,
		      &McCmat.S1_p,
		      &McCmat.ptype_p,
		      &McCmat.qb_p,
		      &McCmat.qm_p,
		      &McCmat.q1k_p,
		      &McCmat.qln_p);
	
	// get pointer to McCaskill base pair probabilities
	McCmat.bppm = export_bppm();
	
	// since the space referenced by pointers in McCmat will be
	// overwritten by the next call to pf_fold, we have to copy
	// the data structures if we want to keep them.
	//
	
	iindx= get_iindx(sequence.length());
	
	unsigned int size;
	size  = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
	S_p   = encode_sequence(c_sequence, 0);
	S1_p  = encode_sequence(c_sequence, 1);
	ptype_p= (char *) space(sizeof(char)*((length+1)*(length+2)/2));
	qb_p= (FLT_OR_DBL *) space(size);
	qm_p= (FLT_OR_DBL *) space(size);
	q1k_p= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
	qln_p= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
	bppm= (FLT_OR_DBL *) space(size);
	
	int i,j;
	for (j=TURN+2;j<=(int)sequence.length(); j++) {
	  for (i=j-TURN-1; i>=1; i--) {
	    ptype_p[iindx[i]-j]= McCmat.ptype_p[iindx[i]-j];
	    qb_p[iindx[i]-j]= McCmat.qb_p[iindx[i]-j];
	    qm_p[iindx[i]-j]= McCmat.qm_p[iindx[i]-j];   
	    bppm[iindx[i]-j]= McCmat.bppm[iindx[i]-j];
	  }
	}
	
	
	for (size_type k=1; k<=sequence.length(); k++) {
	  q1k_p[k]= McCmat.q1k_p[k];
	  qln_p[k]= McCmat.qln_p[k];
	}
	q1k_p[0] = 1.0;
	qln_p[sequence.length()+1] = 1.0;
	
	// copying of McCaskill pf matrices done
	
	
	// precompute tables for computations
	// of probabilities unpaired / basepair in loop or external
	// as they are required for Exparna P functionality
	//
	
	pf_params_p= get_scaled_pf_parameters();

	expMLbase_p= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
	scale_p= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
	
	kT = pf_params_p->kT;   /* kT in cal/mol  */

	/* scaling factors (to avoid overflows) */
	if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
	  pf_scale = exp(-(-185+(pf_params_p->temperature-37.)*7.27)/kT);
	  if (pf_scale<1) pf_scale=1;
	}
	scale_p[0] = 1.;
	scale_p[1] = 1./pf_scale;
	expMLbase_p[0] = 1;
	expMLbase_p[1] = pf_params_p->expMLbase/pf_scale;
	for (size_t i=2; i<=sequence.length(); i++) {
	  scale_p[i] = scale_p[i/2]*scale_p[i-(i/2)];
	  expMLbase_p[i] = pow(pf_params_p->expMLbase, (double)i) * scale_p[i];
	}
	
	compute_Qm2();
	
    }

    
#endif // HAVE_LIBRNA

    RnaData::~RnaData() {
	free_McCaskill_matrices();
    }

    void
    RnaData::free_McCaskill_matrices() {
#     ifdef HAVE_LIBRNA	
	
	if(S_p) free(S_p);          
	if(S1_p) free(S1_p);	 	
	if(ptype_p) free(ptype_p);	 				
	if(qb_p) free(qb_p);					
	if(qm_p) free(qm_p);	 				
	if(q1k_p) free(q1k_p);	 
	if(qln_p) free(qln_p);	  
 	if(bppm) free(bppm);   
	if(iindx) free(iindx);
  	if(qm2) free(qm2);
	if(expMLbase_p) free(expMLbase_p); 
        if(scale_p) free(scale_p); 
	if(pf_params_p) free(pf_params_p);
        
	S_p= NULL;
	S1_p= NULL;
	ptype_p= NULL;
	qb_p= NULL;
	qm_p= NULL;
	q1k_p= NULL;
	qln_p= NULL;
 	bppm= NULL;
	iindx= NULL;
	qm2= NULL;
	pf_params_p= NULL;
	scale_p= NULL;
	expMLbase_p= NULL;
	
#     endif
	
    }

    
    // decide on file format and call either readPS or readPP
    void RnaData::read(const std::string &filename, bool keepMcM) {
  
	std::ifstream in(filename.c_str());
	if (! in.good()) {
	    std::cerr << "Cannot read "<<filename<<std::endl;
	    exit(-1);
	}
    
	std::string s;
	// read first line and decide about file-format
	in >> s;
	in.close();
	if (s == "%!PS-Adobe-3.0") {
	    // try reading as dot.ps file (as generated by RNAfold)
	    readPS(filename);
	    
	    //} else if (s == "<ppml>") {
	    // proprietary PPML file format
	    //readPPML(filename);

#ifdef HAVE_LIBRNA
	} else if (s.substr(0,7) == "CLUSTAL" || s[0]=='>') {
	    // assume multiple alignment format: read and compute base pair probabilities
	    readMultipleAlignment(filename, keepMcM, stacking);
#endif
	} else {
	    // try reading as PP-file (proprietary format, which is easy to read and contains pair probs)
	    readPP(filename);
	}

    
	// DUMP for debugging
	//std::cout << arc_probs_ << std::endl;
	//std::cout << arc_2_probs_ << std::endl;

    }

    void RnaData::readPS(const std::string &filename) {
	std::ifstream in(filename.c_str());
    
	bool contains_stacking_probs=false; // does the file contain probabilities for stacking
    
	std::string s;
	while (in >> s && s!="/sequence") {
	    if (s=="stacked") contains_stacking_probs=true;
	}
    
	if (stacking && ! contains_stacking_probs) {
	    std::cerr << "WARNING: Stacking requested, but no stacking probabilities in dot plot!" << std::endl;
	}

    
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
			else if (stacking && contains_stacking_probs && type=="lbox") { // read a stacking probability
			    set_arc_2_prob(i,j,p); // we store the joint probability of (i,j) and (i+1,j-1)
			}
		    }
		}
	    }
	}
    }

#ifdef HAVE_LIBRNA
    //bool flag;
    void RnaData::readMultipleAlignment(const std::string &filename, bool keepMcM, bool stacking) {
	//read to multiple alignment object
	MultipleAlignment ma(filename,MultipleAlignment::CLUSTAL); // accept clustal input
	// MultipleAlignment ma(filename,MultipleAlignment::FASTA); // accept fasta input
	
	// convert to sequence
	sequence = Sequence(ma);
	
	if (sequence.row_number()!=1) {
	    std::cerr << "ERROR: Cannot handle input from "<<filename<<"."<<std::endl
		      <<"        Base pair computation from multiple sequence alignment is not implemented." << std::endl;
	    exit(-1);
	}
	
	
	// run McCaskill and get access to results
	// in McCaskill_matrices
	compute_McCaskill_matrices();
	
	// initialize the object from base pair probabilities
	// Use the same proability threshold as in RNAfold -p !
	init_from_McCaskill_bppm();
	
	// since we have local copies of all McCaskill pf arrays,
	// we can free the ones of the Vienna lib
	free_pf_arrays();

	if (!keepMcM) {
	    free_McCaskill_matrices();
	}
	
    }

    void
    RnaData::init_from_McCaskill_bppm(double threshold) {
	for( size_t i=1; i <= sequence.length(); i++ ) {
	    for( size_t j=i+1; j <= sequence.length(); j++ ) {
		
		double p= get_bppm(i,j);
		
		if (p >= threshold) { // apply very relaxed filter 
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
    }

    void
    RnaData::compute_Qm2(){

      int len= sequence.length();
      qm1= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(len+2));
      qm2= (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * ((len+1)*(len+2)/2));
      FLT_OR_DBL *qqm1= (FLT_OR_DBL*) space(sizeof(FLT_OR_DBL)*(len+2));
      
      int index_i,index_j,index_k,type;
      for (index_i=1; index_i<=len; index_i++)
	qm1[index_i]=qqm1[index_i]=0;
     
      
      
      for(index_j= TURN+2; index_j<=len; index_j++){
      /*
       *the first inner loop calculates one row of Qm1 that will be needed in the calculation of Qm2  
       */
	for(index_i= index_j-TURN-1; index_i>=1; index_i--){
	 type=get_ptype(index_i,index_j);
	 qm1[index_i]= qqm1[index_i]*expMLbase_p[1];
	 if(type){
 	  qm1[index_i]+= (get_qb(index_i,index_j))* exp_E_MLstem(type, (index_i>1) ? 
									S1_p[index_i-1] : -1, (index_j<len) ? S1_p[index_j+1] : -1,  pf_params_p);

	  }
	  //qqm1=qm1;
	}
	//this part calculates the Qm2 matrix
	if(index_j >= (2*(TURN+2))){
	  for(index_i= index_j-2*TURN-3; index_i>=1; index_i--){
	    qm2[iindx[index_i]-index_j]= 0;
	      for(index_k= index_i+2; index_k< index_j-2; index_k++){
		qm2[iindx[index_i+1]-(index_j-1)]+= get_qm(index_i+1,index_k)*qqm1[index_k+1];
	      
	      }
	  }
	}
	  for(index_i= index_j-TURN-1; index_i>=1; index_i--){
	    qqm1[index_i]=qm1[index_i];
	  }
	}
	free(qm1);
	free(qqm1);
	qm1= NULL;
	qqm1= NULL;
    }

    
    int 
    RnaData::ptype_of_admissible_basepair(size_type i,size_type j) const {
    	int type = (int)(get_ptype(i,j));
	
	// immediately return 0.0 when i and j cannot pair
	if ((type==0)
	    || (((type==3)||(type==4))&&no_closingGU)
	    || (get_qb(i,j)==0.0)
	    || (get_arc_prob(i,j)==0.0))
	    {
		return 0;
	    }
	
	return type;
    }

    double RnaData::prob_unpaired_in_loop(size_type k,size_type i,size_type j) const {
	
	const char *c_sequence=consensus_sequence.c_str();
	
	FLT_OR_DBL H,I,M;
	//calculating the Hairpin loop energy contribution
	
	int type = ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}

	H = exp_E_Hairpin((int)(j-i-1), type, S1_p[(int)(i+1)], S1_p[(int)(j-1)],
			  c_sequence+i-1, pf_params_p) * scale_p[(int)(j-i+1)];
	
	I = 0.0;
	//calculating the Interior loop energy contribution

	int u1;
	// case 1: i<k<i´<j´<j
	for (int ip=(int)(k+1); 
	     ip<=(int)MIN2((int)(i+MAXLOOP+1),(int)(j-TURN-2));
	     ip++) {
	    u1 = (int)(ip-i-1);
	    for (int jp=(int)MAX2((int)(ip+TURN+1),(int)(j-1-MAXLOOP+u1)); 
		 jp<(int)j;
		 jp++) {
		int type2 =(int)(get_ptype(ip,jp));
		if (type2) {
		    type2 = rtype[type2];
		    I += get_qb(ip,jp) 
			* (scale_p[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 S1_p[(int)(i+1)],S1_p[(int)(j-1)],
					 S1_p[ip-1],S1_p[jp+1], pf_params_p));
		}
	    }
	}
	//case 2: i<i´<j´<k<j
	for (int ip=(int)(i+1);
	     ip<=(int)MIN2((int)(i+MAXLOOP+1),(int)(k-TURN-2)); 
	     ip++) {
	    u1 = (int)(ip-i-1);
	    for (int jp=(int)MAX2((int)(ip+TURN+1),(int)(j-1-MAXLOOP+u1));
		 jp<(int)k;
		 jp++) {
		int type2 =(int)(get_ptype(ip,jp)) ;
		if (type2) {
		    type2 = rtype[type2];
		    I += get_qb(ip,jp)
			* (scale_p[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 S1_p[(int)(i+1)],S1_p[(int)(j-1)],
					 S1_p[ip-1],S1_p[jp+1], pf_params_p));
		}
	    }
	}
	
	//calculating Multiple loop energy contribution
	M = 0.0;
	
	M += qm2[iindx[(int)(k+1)]-((int)(j-1))] * expMLbase_p[(int)(k-i)];

	M += qm2[iindx[(int)(i+1)]-((int)(k-1))] * expMLbase_p[(int)(j-k)];
	    
	M += get_qm(i+1,k-1) * expMLbase_p[1] *  get_qm(k+1,j-1);
	    
	// multiply with contribution for closing of multiloop
	M *= pf_params_p->expMLclosing 
	    * exp_E_MLstem(rtype[type],S1_p[(int)(j-1)],S1_p[(int)(i+1)], pf_params_p)
	    * scale_p[2]; 
	
	return ((H+I+M)/get_qb(i,j))*get_arc_prob(i,j);
    }

    double RnaData::prob_unpaired_external(size_type k) const {
	return (q1k_p[k-1] * scale_p[1] * qln_p[k+1]) / qln_p[1];
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
			    S1_p[(int)(i+1)],S1_p[(int)(j-1)],
			    S1_p[ip-1],S1_p[jp+1], pf_params_p)
	    * scale_p[u1+u2+2];
	
	//calculating Multiple loop energy contribution
	//
	Mpp = 0.0;

	// inner base pairs only right of (ip,jp)
	Mpp += expMLbase_p[ip-i-1] * get_qm(jp+1,j-1);
	
	// inner base pairs only left of (ip,jp)
	Mpp += get_qm(i+1,ip-1) * expMLbase_p[j-jp-1];
	
	// inner base pairs left and right of (ip,jp)
	Mpp += get_qm(i+1,ip-1) * get_qm(jp+1,j-1);
	
	// multiply with factor for inner base pair
	Mpp *= exp_E_MLstem(type2, S1_p[ip-1], S1_p[jp+1], pf_params_p);
	
	// multiply with factors for closing base pair
	Mpp *= pf_params_p->expMLclosing
	    * exp_E_MLstem(rtype[type],S1_p[(int)(j-1)],S1_p[(int)(i+1)], pf_params_p)
	    * scale_p[2];
	
	return (get_qb(ip,jp)*(Ipp+Mpp)/get_qb(i,j))*get_arc_prob(i,j);
    }

    double RnaData::prob_basepair_external(size_type i,size_type j) const {
	// immediately return 0.0 when i and j cannot pair
	if (ptype_of_admissible_basepair(i,j)==0) {return 0.0;}
	
	return (q1k_p[i-1] * get_qb(i,j) *  qln_p[j+1]) / qln_p[1];
    }

#endif // HAVE_LIBRNA

    void RnaData::readPP(const std::string &filename) {
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
      
	    if (stacking) {
		double p2;
	  
		if (in >> p2) {
		    set_arc_2_prob(i,j,p2); // p2 is joint prob of (i,j) and (i+1,j-1)
		}
	    }      
	}
    }

	

    /*
      void RnaData::readPPML(const std::string &filename) {

      std::ifstream in(filename.c_str());
    
      std::string tag;
      while (in>>tag) {
      if (tag == "<score>") {
      readScore(in);
      } else if  (tag == "<alignment>") {
      readAlignment(in);
      } else if(tag == "<bpp>") {
      readBPP(in);
      } else if(tag == "<constraints>") {
      readConstraints(in);
      }
      }

      std::string name;
      std::string seqstr;
    
      }
    */

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
