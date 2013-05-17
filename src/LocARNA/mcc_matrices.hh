#ifndef LOCARNA_MCC_MATRICES_HH
#define LOCARNA_MCC_MATRICES_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/params.h>
}

namespace LocARNA {

    class McC_matrices_base {
    protected:
    	size_t length;     //!< sequence length
	bool local_copy; //!< whether pointers point to local copies of data structures
	
	FLT_OR_DBL *qb;  //!< Q<sup>B</sup> matrix					
	FLT_OR_DBL *qm;  //!< Q<sup>M</sup> matrix					
	
	FLT_OR_DBL *bppm;  //!< base pair probability matrix
	
    	int* iindx;        //!< iindx from librna's get_iindx()
    	/** 
	 * @brief construct empty
	 */
	McC_matrices_base();

	/** 
	 * @brief initialize
	 * 
	 * @param length length of sequence 
	 */
	void
	init(size_t length_);

    public:
	FLT_OR_DBL *q1k; //!< 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)	
	FLT_OR_DBL *qln; //!< 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)      
	
	pf_paramT *pf_params; //!< parameters for pf folding
	
		
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual
	~McC_matrices_base();

	
	//! \brief index in triagonal matrix
	size_t iidx(size_t i,size_t j) const {
	    assert(1<=i);
	    assert(i<=j);
	    assert(j<=length);

	    return iindx[i]-j;
	}

	/** 
	 * @brief Access matrix bppm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_bppm(size_t i, size_t j) const { return bppm[iidx(i,j)]; }

	/** 
	 * @brief Access matrix qb
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_qb(size_t i, size_t j) const { return qb[iidx(i,j)]; }

	/** 
	 * @brief Access matrix qm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL get_qm(size_t i, size_t j) const { return qm[iidx(i,j)]; }
	    
    protected:

	void free_all_local();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_base &McCmat);
    };
    
    //! @brief  structure for McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_pf_arrays() and get_bppm() of Vienna librna
    class McC_matrices_t : public McC_matrices_base {
	char *ptype;	   //!< pair type matrix					

    public:

	char *sequence;  //!< 0-terminated sequence string
	short *S;        //!< 'S' array (integer representation of nucleotides)	
	short *S1;	   //!< 'S1' array (2nd integer representation of nucleotides)	
	
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param sequence the sequence as 0-terminated C-string 
	 * @param local_copy  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_matrices_t(char *sequence, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual 
	~McC_matrices_t();


	/** 
	 * @brief Access matrix ptype
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	char get_ptype(size_t i, size_t j) const { return ptype[iidx(i,j)]; }

    protected:

	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_t &McCmat);
    };

     //! @brief  structure for Alifold-McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_alipf_arrays() and alipf_export_bppm() of Vienna librna
    class McC_ali_matrices_t : public McC_matrices_base {
    protected:
	size_t n_seq;     //!< sequence length
	
    public:

	short **S;       //!< 'S' array (integer representation of nucleotides)	
	short **S5;	   //!< 'S5' array
	short **S3;	   //!< 'S3' array
	unsigned short  **a2s;  //!< 'a2s' array
	char **Ss;	   //!< 'Ss' array
    protected:
	short *pscore; //!< alifold covariance/conservation scores
    public:
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param length_ length of sequence 
	 * @param local_copy_  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_ali_matrices_t(size_t n_seq, size_t length, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual
	~McC_ali_matrices_t();


	/** 
	 * @brief Access matrix pscore
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	short get_pscore(size_t i, size_t j) const { return pscore[iidx(i,j)]; }


    protected:

	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_ali_matrices_t &McCmat);
    };
    
} // end namespace LocARNA


#endif // LOCARNA_MCC_MATRICES_HH
