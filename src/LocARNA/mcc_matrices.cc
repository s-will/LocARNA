#include <iostream>
#include <cstdlib> // import free()

#include "aux.hh"


#include "mcc_matrices.hh"

extern "C" {
#include <ViennaRNA/fold_vars.h>
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




namespace LocARNA {
    // ------------------------------------------------------------
    // implementation of classes McC_matrices_base, McC_matrices_t, McC_ali_matrices_t

    void *
    space_memcpy(void *from,size_t size) {
	if (from==NULL) return from;
	void *p = space(size);
	memcpy(p,from,size);
	return p;
    }

    McC_matrices_base::McC_matrices_base()
	: length_(0),
	  local_copy_(false),
	  qb_(0),
	  qm_(0),
	  bppm_(0),
	  iindx_(0),
	  q1k_(0),
	  qln_(0),
	  pf_params_(0)
    {}

    void
    McC_matrices_base::init(size_t length) {
	length_=length;
	
	qb_=0;
	qm_=0;
	bppm_=0;
	q1k_=0;
	qln_=0;
	pf_params_=0;
	
	iindx_ = get_iindx(length);
    }
    
    void
    McC_matrices_base::deep_copy(const McC_matrices_base &McCmat) {
	local_copy_=true;

	length_=McCmat.length_;

	size_t size = sizeof(FLT_OR_DBL) * ((length_+1)*(length_+2)/2);
	
	qb_= (FLT_OR_DBL *) space_memcpy(McCmat.qb_,size);
	qm_= (FLT_OR_DBL *) space_memcpy(McCmat.qm_,size);
	bppm_= (FLT_OR_DBL *) space_memcpy(McCmat.bppm_,size);
	q1k_= (FLT_OR_DBL *) space_memcpy(McCmat.q1k_,sizeof(FLT_OR_DBL)*(length_+1));
	qln_= (FLT_OR_DBL *) space_memcpy(McCmat.qln_,sizeof(FLT_OR_DBL)*(length_+2));
	pf_params_= (pf_paramT *) space_memcpy(McCmat.pf_params_,sizeof(pf_paramT));

	iindx_= get_iindx(length_);
    }

    McC_matrices_base::~McC_matrices_base() {
	if (local_copy_) {
	    free_all_local();
	}
	if (iindx_) free(iindx_);
	if (pf_params_) free(pf_params_);
    }

    void McC_matrices_base::free_all_local() {
	if (qb_) free(qb_);
	if (qm_) free(qm_);
	if (q1k_) free(q1k_);
	if (qln_) free(qln_);
	if (bppm_) free(bppm_);
    }
    
    // ----------------------------------------


    McC_matrices_t::McC_matrices_t(char *sequence, bool local_copy)
	:McC_matrices_base()
    {	
	if (local_copy) {
	    McC_matrices_t McCmat_tmp(sequence,false);
	    deep_copy(McCmat_tmp);
	} else {
	    McC_matrices_base::init(strlen(sequence));
	    
	    this->sequence_=sequence;

	    // get pointers to McCaskill matrices
	    get_pf_arrays(&S_,
			  &S1_,
			  &ptype_,
			  &qb_,
			  &qm_,
			  &q1k_,
			  &qln_);
	    
	    // get pointer to McCaskill base pair probabilities
	    bppm_ = export_bppm();	    
	    
	    pf_params_ = get_scaled_pf_parameters(); //allocates space for pf_params!
	}
    }

    char 
    McC_matrices_t::rev_ptype(size_t i, size_t j) const {
	return rtype[(size_t)ptype(i,j)];
    }

    void
    McC_matrices_t::deep_copy(const McC_matrices_t &McCmat) {
	McC_matrices_base::deep_copy(McCmat);
	
	sequence_ = (char *) space_memcpy(McCmat.sequence_,sizeof(char)*(length_+1));
	S_ = (short *) space_memcpy(McCmat.S_,sizeof(short)*(length_+2));
	S1_ = (short *) space_memcpy(McCmat.S1_,sizeof(short)*(length_+2));
	ptype_= (char *) space_memcpy(McCmat.ptype_,sizeof(char)*((length_+1)*(length_+2)/2));
    }

    McC_matrices_t::~McC_matrices_t() {
	if (local_copy_) {
	    free_all();
	}
    }


    void McC_matrices_t::free_all() {
	free(sequence_);
	free(S_);
	free(S1_);
	free(ptype_);
    }

    // ----------------------------------------
    McC_ali_matrices_t::McC_ali_matrices_t(size_t n_seq, size_t length, bool local_copy)
	: n_seq_(n_seq)
    {	
	if (local_copy) {
	    McC_ali_matrices_t McCmat_tmp(n_seq,length,false);
	    deep_copy(McCmat_tmp);
	} else {
	    McC_matrices_base::init(length);
	    
	    // get pointers to McCaskill matrices
	    get_alipf_arrays(&S_,
			     &S5_,
			     &S3_,
			     &a2s_,
			     &Ss_,
			     &qb_,
			     &qm_,
			     &q1k_,
			     &qln_,
			     &pscore_);

	    // get pointer to McCaskill base pair probabilities
	    bppm_ = alipf_export_bppm();	    
	    
	    pf_params_ = get_scaled_alipf_parameters(n_seq); // allocates space
	}
    }
    

    void
    McC_ali_matrices_t::deep_copy(const McC_ali_matrices_t &McCmat) {
	McC_matrices_base::deep_copy(McCmat);
		
	n_seq_ = McCmat.n_seq_;

	S_    = (short **)          space(n_seq_ * sizeof(short *));
	S5_   = (short **)          space(n_seq_ * sizeof(short *));
	S3_   = (short **)          space(n_seq_ * sizeof(short *));
	a2s_  = (unsigned short **) space(n_seq_ * sizeof(unsigned short *));
	Ss_   = (char **)           space(n_seq_ * sizeof(char *));

	for (size_t i=0; i<n_seq_; i++) {
	    S_[i]   = (short *)          space_memcpy(McCmat.S_[i],  (length_+2) * sizeof(short));
	    S5_[i]  = (short *)          space_memcpy(McCmat.S5_[i], (length_+2) * sizeof(short));
	    S3_[i]  = (short *)          space_memcpy(McCmat.S3_[i], (length_+2) * sizeof(short));
	    a2s_[i] = (unsigned short *) space_memcpy(McCmat.a2s_[i],(length_+2) * sizeof(unsigned short));
	    Ss_[i]  = (char *)           space_memcpy(McCmat.Ss_[i], (length_+2) * sizeof(char));
	}

	pscore_ = (short*) space_memcpy(McCmat.pscore_,
				       ((length_+1)*(length_+2))/2 * sizeof(short));
    }

    McC_ali_matrices_t::~McC_ali_matrices_t() {
	if (local_copy_) {
	    free_all();
	}
    }


    void McC_ali_matrices_t::free_all() {
	free_sequence_arrays(n_seq_,&S_,&S5_,&S3_,&a2s_,&Ss_);
	if (pscore_) free(pscore_);
    }

} // end namespace LocARNA
