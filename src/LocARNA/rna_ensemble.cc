#include <fstream>
#include <sstream>
#include <map>
#include <limits>

#include "aux.hh"

#include "rna_ensemble_impl.hh"

#include "alphabet.hh"

#include "multiple_alignment.hh"

#include "global_stopwatch.hh"


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
    // implementation of class RnaEnsemble
    //

    RnaEnsemble::RnaEnsemble(const Sequence &sequence,
			     const PFoldParams &params,
			     bool inLoopProbs, 
			     bool use_alifold)
	: pimpl_(new RnaEnsembleImpl(this,
				     sequence,
				     params,
				     inLoopProbs,
				     use_alifold)) {
    }
    
    RnaEnsembleImpl::RnaEnsembleImpl(RnaEnsemble *self,
				     const Sequence &sequence,
				     const PFoldParams &params,
				     bool inLoopProbs, 
				     bool use_alifold=true)
	: self_(self),
	  sequence_(sequence),
	  pair_probs_available_(false),
	  stacking_probs_available_(false),
	  in_loop_probs_available_(false),
	  McCmat_(0L), // 0 pointer
	  used_alifold_(false),
	  min_free_energy_(std::numeric_limits<double>::infinity()),
	  min_free_energy_structure_("") 
    {
	compute_ensemble_probs(params,inLoopProbs,use_alifold);
    }
    
    RnaEnsemble::~RnaEnsemble() {
	delete pimpl_;
    }

    RnaEnsembleImpl::~RnaEnsembleImpl() {
#     ifdef HAVE_LIBRNA
	if (McCmat_) {
	    delete McCmat_;
	}
#     endif
    }

    bool
    RnaEnsemble::has_base_pair_probs() const {
	return pimpl_->pair_probs_available_;
    }

    bool
    RnaEnsemble::has_stacking_probs() const {
	return pimpl_->stacking_probs_available_;
    }

    bool
    RnaEnsemble::has_in_loop_probs() const {
	return pimpl_->in_loop_probs_available_;
    }


    const Sequence &
    RnaEnsemble::sequence() const {
	return pimpl_->sequence_;
    }

    double
    RnaEnsemble::min_free_energy() const {
	return pimpl_->min_free_energy_;
    }

    std::string
    RnaEnsemble::min_free_energy_structure() const {
	return pimpl_->min_free_energy_structure_;
    }
    
    //! \brief get length of sequence
    //! \return sequence length
    size_type
    RnaEnsemble::length() const {
	return pimpl_->sequence_.length();
    }
    
    

#ifdef HAVE_LIBRNA

    void
    RnaEnsembleImpl::compute_ensemble_probs(const PFoldParams &params,bool inLoopProbs, bool use_alifold) {
	
	stopwatch.start("bpp");
	
	assert(use_alifold || pimpl_->sequence_.row_number()==1);

	used_alifold_=use_alifold;

	// run McCaskill and get access to results
	// in McCaskill_matrices
	if (!use_alifold) {
	    compute_McCaskill_matrices(params,inLoopProbs);
	} else {
	    make_pair_matrix();
	    compute_McCaskill_alifold_matrices(params,inLoopProbs);
	}
	
	throw failure("ACHTUNG DAS KLINGT GEFAEHRLICH");
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
    RnaEnsembleImpl::compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs) {
	assert(sequence_.row_number()==1);

	// use MultipleAlignment to get pointer to c-string of the
	// first (and only) sequence in object sequence.
	//
	size_t length = sequence_.length();
	
	char c_sequence[length+1];
	std::string seqstring = 
	    MultipleAlignment(sequence_).seqentry(0).seq().to_string();
	
	strcpy(c_sequence,seqstring.c_str());
	
	char c_structure[length+1];
	
	// ----------------------------------------
	// set folding parameters
	if (params.noLP()) {noLonelyPairs=1;}
	
	
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
	McCmat_ = 
	    new McC_matrices_t(c_sequence,inLoopProbs); // makes local copy (if needed)
	
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
	    kT = McCmat_->pf_params_->kT;   /* kT in cal/mol  */
	    
	    /* scaling factors (to avoid overflows) */
	    if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
		pf_scale = 
		    exp(-(-185+(McCmat_->pf_params_->temperature-37.)*7.27)/kT);
		if (pf_scale<1) pf_scale=1;
	    }
	    
	    scale_[0] = 1.;
	    scale_[1] = 1./pf_scale;
	    expMLbase_[0] = 1;
	    expMLbase_[1] = 
		McCmat_->pf_params_->expMLbase * scale_[1];
	    
	    for (size_t i=2; i<=sequence_.length(); i++) {
		
		// due to the folowing: scale_[i] = pow(scale_[1],(double)i)
		scale_[i] = 
		    scale_[i/2]*scale_[i-(i/2)];
 
		expMLbase_[i] = 
		    pow(McCmat_->pf_params_->expMLbase, (double)i) * scale_[i];
	    }
	
	    
	    // ----------------------------------------
	    // compute the Qm2 matrix
	    compute_Qm2();
	}
    }


    void
    RnaEnsembleImpl::compute_McCaskill_alifold_matrices(const PFoldParams &params, bool inLoopProbs) {
	
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
	if (params.noLP()) {noLonelyPairs=1;}
	
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
	McCmat_ =
	    new McC_ali_matrices_t(n_seq,length,inLoopProbs); // makes local copy (if needed)
	
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
	    double scaling_factor=McCmat_->pf_params_->pf_scale;
	    kT = McCmat_->pf_params_->kT / n_seq;   /* kT in cal/mol  */
	    
	    
	    /* scaling factors (to avoid overflows) */
	    if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
		scaling_factor = 
		    exp(-(-185+(McCmat_->pf_params_->temperature-37.)*7.27)/kT);
		if (scaling_factor<1) scaling_factor=1;
		McCmat_->pf_params_->pf_scale=scaling_factor;
	    }
	    scale_[0] = 1.;
	    scale_[1] = 1./scaling_factor;

	    expMLbase_[0] = 1;
	    expMLbase_[1] = 
		McCmat_->pf_params_->expMLbase/scaling_factor;
	    for (size_t i=2; i<=sequence_.length(); i++) {
		scale_[i] = 
		    scale_[i/2]*scale_[i-(i/2)];
		expMLbase_[i] = 
		    pow(McCmat_->pf_params_->expMLbase, (double)i) * scale_[i];
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
    
    void
    RnaEnsembleImpl::compute_Qm2(){
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
		char type=MCm->ptype(i,j);
		qqm[i]= qqm1[i]*expMLbase_[1];
		if(type) {
		    qqm[i] +=
			MCm->qb(i,j)
			* exp_E_MLstem(type,
				       (i>1)   ? MCm->S1_[i-1] : -1, 
				       (j<len) ? MCm->S1_[j+1] : -1,  
				       MCm->pf_params_);
		}
		
		//qm1[McCmat->iidx(i,j)]=qqm[i];

		assert(qqm[i] <= MCm->qm(i,j));
		assert((!frag_len_geq(i,j-1,TURN+2)) || qqm1[i] <= MCm->qm(i,j-1));
	    }
	    	    
	    // --------------------
	    // calculates column "j" of the Qm2 matrix
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*(TURN+2)+1; i>=1; i--) {
		    qm2_[MCm->iidx(i,j)] = 0;
		    for(size_type k = i + TURN+1; (k+1)+TURN+1 <= j; k++) {
			qm2_[MCm->iidx(i,j)] +=
			    MCm->qm(i,k)*qqm[k+1];
		    }
		    assert(qm2_[MCm->iidx(i,j)] <= MCm->qm(i,j));
		}
	    }
	    	    
	    // --------------------
	    // swap qqm and qqm1 (in constant time)
	    qqm1.swap(qqm);
	}
    }


    void
    RnaEnsembleImpl::compute_Qm2_ali(){
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
		    type[s] = pair[MCm->S_[s][i]][MCm->S_[s][j]];
		    if (type[s]==0) type[s]=7;
		}
		
		qqm[i]= qqm1[i]*expMLbase_[1];
		
		FLT_OR_DBL qbt1=1.0; // collects contribution "inner basepair of multiloop"
		for (size_t s=0; s<n_seq; s++) {
		    qbt1 *= exp_E_MLstem(type[s], 
					 i>1 ? MCm->S5_[s][i] : -1,
					 j<len ? MCm->S3_[s][j] : -1,
					 MCm->pf_params_);
		}
		qqm[i] += MCm->qb(i,j) * qbt1;

	    }
	    
	    // --------------------
	    // calculate a row of the matrix Qm2
	    //
	    if(j >= (2*(TURN+2))) {
		for(size_type i = j-2*TURN-3; i>=1; i--) {
		    qm2_[MCm->iidx(i+1,j-1)] = 0;
		    for(size_type k = i+TURN+2; k< j-TURN-2; k++) {
			qm2_[MCm->iidx(i+1,j-1)] +=
			    MCm->qm(i+1,k)*qqm1[k+1];
		    }
		}
	    }
	    
	    // --------------------
	    // swap row qqm and qqm1 (in constant time)
	    qqm1.swap(qqm);
	}
    }

    int
    RnaEnsembleImpl::ptype_of_admissible_basepair(size_type i,size_type j) const {
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);

    	int type = MCm->ptype(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if ((type==0)
	    || (((type==3)||(type==4))&&no_closingGU)
	    || (MCm->qb(i,j)==0.0)
	    || (self_->arc_prob(i,j)==0.0))
	    {
		return 0;
	    }
	
	return type;
    }

    double 
    RnaEnsemble::arc_prob(size_type i, size_type j) const {
	return pimpl_->McCmat_->bppm(i,j);
    }
    
    double
    RnaEnsemble::arc_2_prob(size_type i, size_type j) const {
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(pimpl_->McCmat_);

	if ( MCm->qb(i+1,j-1) == 0.0 ) {
	    return 0.0;
	}
	
	FLT_OR_DBL p = arc_prob(i,j);
	p *= MCm->qb(i+1,j-1)/MCm->qb(i,j);
	p *= exp_E_IntLoop(0,0,MCm->ptype(i,j),
			   MCm->rev_ptype(i+1,j-1),
			   0,0,0,0, MCm->pf_params_)*pimpl_->scale_[2];
	
	return p;
    }

    double
    RnaEnsembleImpl::unpaired_in_loop_prob_ali(size_type k,size_type i,size_type j) const {
    	assert(frag_len_geq(i,j,TURN+2));
	assert(i<k);
	assert(k<j);
	assert(in_loop_probs_available_);
	
	McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t*>(this->McCmat_);
	
        size_t n_seq = sequence_.row_number();
	
	// immediately return 0.0 if i and j do not pair
	if (self_->arc_prob(i,j)==0.0 || MCm->qb(i,j)==0.0) {return 0.0;}
	
	// get base pair types for i,j of all sequences
	std::vector<int> type(n_seq);
	
	for (size_t s=0; s<n_seq; ++s) {
	    type[s] = pair[MCm->S_[s][i]][MCm->S_[s][j]];
	    if (type[s]==0) type[s]=7;
	}

	// ------------------------------------------------------------
	// hairpin contribution
        //
	
	FLT_OR_DBL H=1.0;
	
	for (size_t s=0; s<n_seq; s++) {
	    size_t u = MCm->a2s_[s][j-1]-MCm->a2s_[s][i];
	    if (MCm->a2s_[s][i]<1) continue;
	    char loopseq[10];
	    if (u<7){
		strncpy(loopseq, MCm->Ss_[s]+MCm->a2s_[s][i]-1, 10);
	    }
	    H *= exp_E_Hairpin(u, type[s],
			       MCm->S3_[s][i], MCm->S5_[s][j], 
			       loopseq, 
			       MCm->pf_params_);
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
		
		if (MCm->qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = MCm->a2s_[s][ip-1] - MCm->a2s_[s][i];
		    size_t u2 = MCm->a2s_[s][j-1] - MCm->a2s_[s][jp];
		    
		    int type_2 = pair[MCm->S_[s][jp]][MCm->S_[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    MCm->S3_[s][i],
					    MCm->S5_[s][j],
					    MCm->S5_[s][ip],
					    MCm->S3_[s][jp],
					    MCm->pf_params_
					    );
		}
		
		I += MCm->qb(ip,jp) * scale_[ip-i+j-jp] * qloop;
	    }
	}
 
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip <= std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    for (size_t jp = std::max(ip+TURN+1 + MAXLOOP,j-1+ ip-i-1) - MAXLOOP; jp<k; jp++) {
		
		FLT_OR_DBL qloop=1.0;
		
		if (MCm->qb(ip,jp)==0) {
		    continue;
		}
		
		for (size_t s=0; s<n_seq; s++) {
		    size_t u1 = MCm->a2s_[s][ip-1] - MCm->a2s_[s][i];
		    size_t u2 = MCm->a2s_[s][j-1] - MCm->a2s_[s][jp];
		    
		    int type_2 = pair[MCm->S_[s][jp]][MCm->S_[s][ip]]; 
		    if (type_2 == 0) type_2 = 7;
		    
		    qloop *= exp_E_IntLoop( u1, u2,
					    type[s], type_2,
					    MCm->S3_[s][i],
					    MCm->S5_[s][j],
					    MCm->S5_[s][ip],
					    MCm->S3_[s][jp],
					    MCm->pf_params_
					    );
		}
		
		I += MCm->qb(ip,jp) * scale_[ip-i+j-jp] * qloop;
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
	    M += MCm->qm(i+1,k-1) * expMLbase_[1] *  MCm->qm(k+1,j-1);
	}
	
	// multiply with contribution for closing of multiloop
	
	for (size_t s=0; s<n_seq; s++) {
	    int tt = rtype[type[s]];
	    
	    M *= MCm->pf_params_->expMLclosing 
		* exp_E_MLstem(tt,MCm->S5_[s][j],MCm->S3_[s][i], MCm->pf_params_);
	}
	M *= scale_[2];
	
	FLT_OR_DBL Qtotal=H+I+M;
	
	double kTn   = MCm->pf_params_->kT/10.;   /* kT in cal/mol  */
	
	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	Qtotal *= exp(MCm->pscore(i,j)/kTn);
	
	FLT_OR_DBL p_k_cond_ij = Qtotal/MCm->qb(i,j); 
	
	FLT_OR_DBL res = p_k_cond_ij * MCm->bppm(i,j);
	
	return res;
    }
    
    double 
    RnaEnsemble::unpaired_in_loop_prob(size_type k,size_type i,size_type j) const {
	assert(i+TURN+1 <= j);
	assert(i<k);
	assert(k<j);
	
	if (!pimpl_->in_loop_probs_available_) return 1.0;
	
	
	if (pimpl_->used_alifold_) {
	    return pimpl_->unpaired_in_loop_prob_ali(k, i, j);
	} else {
	    return pimpl_->unpaired_in_loop_prob_noali(k, i, j);
	}
    }
    
    double 
    RnaEnsembleImpl::unpaired_in_loop_prob_noali(size_type k,size_type i,size_type j) const {
	assert(!used_alifold_);
	assert(in_loop_probs_available_);

	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(McCmat_);

	const char *c_sequence = MCm->sequence_;
	
	
	int type = ptype_of_admissible_basepair(i,j);
	
	// immediately return 0.0 when i and j cannot pair
	if (type==0) {return 0.0;}

	// ------------------------------------------------------------
	// Hairpin loop energy contribution
	
	size_t u=j-i-1;
	FLT_OR_DBL H = exp_E_Hairpin(u, type, MCm->S1_[i+1], MCm->S1_[j-1],
				     c_sequence+i-1, MCm->pf_params_) * scale_[u+2];
	
	// ------------------------------------------------------------
	// Interior loop energy contribution
	FLT_OR_DBL I = 0.0;

	// case 1: i<k<i´<j´<j
	for (size_t ip=k+1; ip <= std::min(i+MAXLOOP+1,j-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp = std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<j; jp++) {
		int type2 = MCm->ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += MCm->qb(ip,jp) 
			* (scale_[u1+j-jp+1] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 MCm->S1_[i+1],MCm->S1_[j-1],
					 MCm->S1_[ip-1],MCm->S1_[jp+1], MCm->pf_params_));
		}
	    }
	}
	//case 2: i<i´<j´<k<j
	for (size_t ip=i+1; ip<=std::min(i+MAXLOOP+1,k-TURN-2); ip++) {
	    size_t u1 = ip-i-1;
	    for (size_t jp=std::max(ip+TURN+1+MAXLOOP,j-1+u1)-MAXLOOP; jp<k; jp++) {
		int type2 = MCm->ptype(ip,jp);
		if (type2) {
		    type2 = rtype[type2];
		    I += MCm->qb(ip,jp)
			* (scale_[(int)(u1+j-jp+1)] *
			   exp_E_IntLoop(u1,(int)(j-jp-1), type, type2,
					 MCm->S1_[i+1],MCm->S1_[j-1],
					 MCm->S1_[ip-1],MCm->S1_[jp+1], MCm->pf_params_));
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
	    M3 = MCm->qm(i+1,k-1) * expMLbase_[1] *  MCm->qm(k+1,j-1);
	}
	
	M=M1+M2+M3;
	
	// multiply with contribution for closing of multiloop
	M *= MCm->pf_params_->expMLclosing 
	    * exp_E_MLstem(rtype[type],MCm->S1_[j-1],MCm->S1_[i+1], MCm->pf_params_)
	    * scale_[2];
	
	FLT_OR_DBL Qtotal = H+I+M;

	FLT_OR_DBL p_k_cond_ij = Qtotal/MCm->qb(i,j); 

	FLT_OR_DBL res = p_k_cond_ij * MCm->bppm(i,j);
	
	return res;
    }

    double 
    RnaEnsemble::unpaired_external_prob(size_type k) const {
	assert(1<=k);
	assert(k<=pimpl_->sequence_.length());
	
	if (!pimpl_->in_loop_probs_available_) return 1.0;
	
	return 
	    (pimpl_->McCmat_->q1k_[k-1] * pimpl_->scale_[1] * pimpl_->McCmat_->qln_[k+1])
	    / pimpl_->McCmat_->qln_[1];
    }


    double
    RnaEnsembleImpl::arc_in_loop_prob_ali(size_type ip,
						size_type jp,
						size_type i,
						size_type j) const {
	assert(in_loop_probs_available_);
	
	McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t *>(this->McCmat_);
	
	size_t n_seq = sequence_.row_number();
	
	// note: the following tests cover the case that the distances of i,j and ip,jp are too small

	// immediately return 0.0 if i and j do not pair
	if (self_->arc_prob(i,j)==0.0 || MCm->qb(i,j)==0.0) {return 0.0;}
	
	// immediately return 0.0 when ip and jp cannot pair
	if (self_->arc_prob(ip,jp)==0.0 || MCm->qb(ip,jp)==0.0) {return 0.0;}
	
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
	    type[s] = pair[MCm->S_[s][i]][MCm->S_[s][j]];
	    if (type[s]==0) type[s]=7;

	    type2[s] = pair[MCm->S_[s][ip]][MCm->S_[s][jp]];
	    if (type2[s]==0) type2[s]=7;
	}


	// note: I and M are computed without factor qb(ip,jp),
	// which is multiplied only in the end.
	
	// ------------------------------------------------------------
	// Interior loop energy contribution
	//
	FLT_OR_DBL I=0.0;
	
	if ((frag_len(i,ip)+frag_len(jp,j))<=MAXLOOP) {
	    I = 1.0;
	    for (size_t s=0; s<n_seq; s++) {
		
		size_t u1 = MCm->a2s_[s][ip-1] - MCm->a2s_[s][i];
		size_t u2 = MCm->a2s_[s][j-1] - MCm->a2s_[s][jp];
		
		I *= exp_E_IntLoop(u1,u2,
				   type[s], rtype[type2[s]],
				   MCm->S3_[s][i],
				   MCm->S5_[s][j],
				   MCm->S5_[s][ip],
				   MCm->S3_[s][jp],
				   MCm->pf_params_);
	    }
	    I *= scale_[ip-i+j-jp];
	}
	
	// ------------------------------------------------------------
	// Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase_[frag_len(i+1,ip-1)] * MCm->qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += MCm->qm(i+1,ip-1) * expMLbase_[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += MCm->qm(i+1,ip-1) * MCm->qm(jp+1,j-1);
	}
	
	for (size_t s=0; s<n_seq; s++) {
	    // multiply with factor for inner base pair
	    M *= exp_E_MLstem(type2[s],
			      MCm->S5_[s][ip],
			      MCm->S3_[s][jp],
			      MCm->pf_params_);
	    
	    // multiply with factors for closing base pair
	    M *= MCm->pf_params_->expMLclosing
		* exp_E_MLstem(rtype[type[s]],
			       MCm->S5_[s][j],
			       MCm->S3_[s][i],
			       MCm->pf_params_);
	}
	
	M *= scale_[2]; // scale for closing base pair
	
	// ------------------------------------------------------------
	FLT_OR_DBL Qtotal= (I+M);
	
	Qtotal *= MCm->qb(ip,jp);

	// multiply with pscore contribution for closing base pair (i,j),
	// like in the calculation of Qb(i,j)
	double kTn   = MCm->pf_params_->kT/10.;   /* kT in cal/mol  */
	Qtotal *= exp(MCm->pscore(i,j)/kTn);

	
	return
	    (Qtotal/MCm->qb(i,j))
	    *MCm->bppm(i,j);
    }

    double
    RnaEnsemble::arc_in_loop_prob(size_type ip,
				  size_type jp,
				  size_type i,
				  size_type j) const {
	
	if (!pimpl_->in_loop_probs_available_) return 1.0;

	if (pimpl_->used_alifold_) {
	    return pimpl_->arc_in_loop_prob_ali(ip, jp, i, j);
	} else {
	    return pimpl_->arc_in_loop_prob_noali(ip, jp, i, j);
	}
    }
    

    double
    RnaEnsembleImpl::arc_in_loop_prob_noali(size_type ip,
					    size_type jp,
					    size_type i,
					    size_type j) const {
	
	assert(in_loop_probs_available_);
	assert(!used_alifold_);
	McC_matrices_t *MCm = static_cast<McC_matrices_t *>(this->McCmat_);
	

	// note: I and M are computed without factor qb(ip,jp),
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
			      MCm->S1_[(int)(i+1)],
			      MCm->S1_[(int)(j-1)],
			      MCm->S1_[ip-1],
			      MCm->S1_[jp+1],
			      MCm->pf_params_)
		* scale_[u1+u2+2];
	}
	
	//calculating Multiple loop energy contribution
	//
	FLT_OR_DBL M = 0.0;

	
	// inner base pairs only right of (ip,jp)
	if ( frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += expMLbase_[frag_len(i+1,ip-1)] * MCm->qm(jp+1,j-1);
	}
	
	// inner base pairs only left of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) ) {
	    M += MCm->qm(i+1,ip-1) * expMLbase_[frag_len(jp+1,j-1)];
	}
	
	// inner base pairs left and right of (ip,jp)
	if ( frag_len_geq(i+1, ip-1, TURN+2) && frag_len_geq(jp+1, j-1, TURN+2) ) {
	    M += MCm->qm(i+1,ip-1) * MCm->qm(jp+1,j-1);
	}
	
	// multiply with factor for inner base pair
	M *= exp_E_MLstem(type2,
			  MCm->S1_[ip-1],
			  MCm->S1_[jp+1],
			  MCm->pf_params_);
	
	// multiply with factors for closing base pair
	M *= MCm->pf_params_->expMLclosing
	    * exp_E_MLstem(rtype[type],
			   MCm->S1_[j-1],
			   MCm->S1_[i+1],
			   MCm->pf_params_)
	    * scale_[2];
	
	FLT_OR_DBL Qtotal = I+M;
	
	Qtotal *= MCm->qb(ip,jp);
	
	
	return Qtotal/MCm->qb(i,j)
	    * MCm->bppm(i,j);
    }

    double RnaEnsemble::arc_external_prob(size_type i,size_type j) const {

	if (!pimpl_->in_loop_probs_available_) return 1.0;
	
	size_t n=pimpl_->sequence_.length();

	assert(1<=i);
	assert(i<j);
	assert(j<=n);
	assert(pimpl_->frag_len_geq(i,j,TURN+2));
	
	// immediately return 0.0 when i and j cannot pair
	if (arc_prob(i,j)==0.0 || pimpl_->McCmat_->qb(i,j)==0.0) {
	    return 0.0;
	}
	
	FLT_OR_DBL extloop;
	
	if (!pimpl_->used_alifold_) {
	    McC_matrices_t *MCm = static_cast<McC_matrices_t *>(pimpl_->McCmat_);
	    extloop = exp_E_ExtLoop(MCm->ptype(i,j),
				    i>1 ? MCm->S1_[i-1] : -1, 
				    j<n ? MCm->S1_[j+1] : -1, 
				    MCm->pf_params_);
	} else {
	    McC_ali_matrices_t *MCm = static_cast<McC_ali_matrices_t *>(pimpl_->McCmat_);
	    
	    size_t n_seq=pimpl_->sequence_.row_number();
	    
	    extloop=1.0;
	    
	    for (size_t s=0; s<n_seq; s++) {
		int type = pair[MCm->S_[s][i]][MCm->S_[s][j]];
		if (type==0) type=7;
		
		extloop *= exp_E_ExtLoop(type, i>1 ? MCm->S5_[s][i] : -1,
					 j<n ? MCm->S3_[s][j] : -1,
					 MCm->pf_params_);
	    }
	}

	return
	    (pimpl_->McCmat_->q1k_[i-1]
	     * pimpl_->McCmat_->qb(i,j)
	     * extloop
	     * pimpl_->McCmat_->qln_[j+1]
	     )
	    / pimpl_->McCmat_->qln_[1];

    }
    

#endif // HAVE_LIBRNA
    
} // end namespace LocARNA 
