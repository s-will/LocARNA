#ifdef HAVE_LIBRNA

#include <iostream>

#include "aux.hh"


#include "mcc_matrices.hh"

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

} // end namespace LocARNA


#   endif // HAVE_LIBRNA
