#include "aligner_p.hh"

#include "arc_matches.hh"
#include "trace_controller.hh"

#include <math.h>
#include <assert.h>
#include <iomanip>
// #include <queue>

using namespace std;

// ------------------------------------------------------------
// AlignerP: compute partition function and probabilities of arc matchs and base matchs
//


// allocate space for the inside matrices 
void
AlignerP::alloc_inside_matrices() {
    Dmat.resize(bpsA.num_bps(), bpsB.num_bps());
    Dmat.fill((pf_score_t )0); // this is essential, such that we can avoid to test validity of arc matches 
        
    //std::cout << "Size of Dmat:" << sizeof(Dmat)+bpsA.num_bps()*bpsB.num_bps()*sizeof(pf_score_t) << std::endl;
  
    M.resize(seqA.length()+1, seqB.length()+1);
    M.fill((pf_score_t )0);
    
    //std::cout << "Size of M:" << sizeof(M)+(seqA.length()+1)*(seqB.length()+1)*sizeof(pf_score_t) << std::endl;
    
    
    E.resize(seqB.length()+1); // size: one row of M/Mprime matrix
}


// allocate space for the inside matrices 
void
AlignerP::alloc_outside_matrices() {

    Dmatprime.resize(bpsA.num_bps(), bpsB.num_bps());
    Dmatprime.fill((pf_score_t )0);
  

    Mprime.resize(seqA.length()+1, seqB.length()+1);
    Mprime.fill((pf_score_t )0);
  
    Eprime.resize(seqB.length()+1); // size: one row of M/Mprime matrix
    


    Mrev.resize(seqA.length()+1, seqB.length()+1);
    Erev.resize(seqB.length()+1); // size: one row of M/Mprime matrix
    
    Erev_mat.resize(seqA.length()+1, seqB.length()+1); // size as Mrev
    Frev_mat.resize(seqA.length()+1, seqB.length()+1); // size as Mrev

}


//===========================================================================
AlignerP::AlignerP(const Sequence &seqA_,
		   const Sequence &seqB_,
		   const ArcMatches &arc_matches_,
		   const AlignerParams *ap_, 
		   const Scoring *s_,
		   const pf_score_t pf_scale_
		   ) :
    scoring(s_),
    params(ap_),
    seqA(seqA_),
    bpsA(arc_matches_.get_base_pairsA()),
    seqB(seqB_),
    bpsB(arc_matches_.get_base_pairsB()),
    arc_matches(arc_matches_), 
    r(1, 1, seqA_.length(), seqB_.length()),
    pf_scale(pf_scale_),
    //Dmat((pf_score_t )0),
    //Dmatprime((pf_score_t )0),
    am_prob(0.0),
    bm_prob(0.0),
    D_created(false),
    Dprime_created(false)
{
    
}


// ====================================================================================================
// ====================================================================================================
// INSIDE ALGORITHM
// ====================================================================================================
// ====================================================================================================


//===========================================================================
// Initialization of inside matrices  M,E,F
//

// Initialize the "0th" row and column of the M matrix, i.e. actually
// this initializes row al and column bl of the M matrix for the
// alignment below of the arcs starting at al and bl. Additionally,
// initialize all (according to trace controller) invalid matrix
// entries that can be accessed from valid ones with 0.
//
void AlignerP::init_M(size_type al, size_type ar, size_type bl, size_type br) {
    
    M(al, bl)=((pf_score_t)1)/pf_scale; // empty alignment

    pf_score_t indel_score;

    // initialize column bl of M
    indel_score = scoring->exp_indel_opening()/pf_scale;
    size_type i;
    for (i=al+1; i<ar; i++) {
	if (params->trace_controller.min_col(i)>bl) break; // fill only as long as column bl is accessible

	indel_score *= scoring->exp_gapA(i, bl);
	M(i, bl) = indel_score;
    }

    // fill entries left of valid entries 
    for ( ; i<ar; i++) {
	assert(params->trace_controller.min_col(i)>bl);
	M(i,params->trace_controller.min_col(i)-1) = 0;
    }
    
    // initialize row al of M
    size_type max_col = std::min(br-1,params->trace_controller.max_col(al));
    indel_score = scoring->exp_indel_opening()/pf_scale;
    size_type j;
    for (j=bl+1; j<=max_col; j++) {
	indel_score *= scoring->exp_gapB(al, j);
	M(al, j) = indel_score;
    }
    // fill entries above valid entries 
    // here j points to one position right of the last initialized entry in row al
    for (size_type i=al+1; i<ar; i++) {
	for (; j<min(br,params->trace_controller.max_col(i)+1); ++j) {
	    M(i-1,j)=0;
	}
    }
}

void AlignerP::init_E(size_type al, size_type ar, size_type bl, size_type br) {
    //
    // all entries are 0 initially, since there is no alignment
    // with empty subsequence of seqA that ends with a gapped base of seqA 
    //
    for (size_type j=bl; j<br; j++) {
	E[j] = (pf_score_t)0;
    }
}

//===========================================================================	
// recursive computation of single entries in inside matrices M,E,F
//

// compute entry in E-vector for (i,j)
inline
pf_score_t
AlignerP::comp_E_entry(size_type al, size_type bl, size_type i, size_type j) {
    return 
	E[j] * scoring->exp_gapA(i, j)
	+ 
	(M(i-1,j)-E[j]) * scoring->exp_gapA(i, j) * scoring->exp_indel_opening();
}

// compute entry in F for (i,j)
inline
pf_score_t
AlignerP::comp_F_entry(size_type al, size_type bl, size_type i, size_type j) {
    return 
	F * scoring->exp_gapB(i, j)
	+
	(M(i,j-1)-F) * scoring->exp_gapB(i, j) * scoring->exp_indel_opening();
}


// recursively calculate the entry M(i,j) handling the
// standard cases in alignment: base match, base in/del, arc match
// for the alignment below of arcs starting at (al,bl)
// pre: E and F entry is already computed
inline
pf_score_t
AlignerP::comp_M_entry(size_type al, size_type bl, size_type i, size_type j) {
    
    pf_score_t pf;
    pf = 
	// base match
	M(i-1, j-1) * scoring->exp_basematch(i, j)
    
	// base del
	+ E[j]
	
	// base ins
	+ F;
    
    // --------------------
    // arc match
    
    // standard case for arc match (without restriction to lonely pairs)
    
    const BasePairs::RightAdjList &adjlA = bpsA.right_adjlist(i);
    const BasePairs::RightAdjList &adjlB = bpsB.right_adjlist(j);
    
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    for (BasePairs::RightAdjList::const_iterator arcA=adjlA.begin(); 
	 arcA !=adjlA.end() && arcA->left() > al; ++arcA) {
	for (BasePairs::RightAdjList::const_iterator arcB=adjlB.begin(); 
	     arcB !=adjlB.end() && arcB->left() > bl; ++arcB) {
	    
	    // consider score for match of basepairs
	    //assert(M(arcA->left()-1, arcB->left()-1) > 0);
      	    
	    pf += M(arcA->left()-1, arcB->left()-1) * D(*arcA, *arcB) * pf_scale;
	    // note: disallowed arc matchs (due to heuristic) are
	    // handled correctly, since there D(arcA->idx(), arcB->idx()) was set to 0
	}
    }
    
    return pf;
}


//===========================================================================
// computation of all entries in inside matrices M,E,F inside of one arc pair
//

void AlignerP::align_inside_arcmatch(size_type al, size_type ar, size_type bl, size_type br) {
    
    //initialize M matrix
    init_M(al, ar, bl, br);
    
    //initialize E vector
    init_E(al, ar, bl, br);

    for (size_type i=al+1; i<ar; i++) {
	F = (pf_score_t)0; // init F
    
	// limit entries due to trace controller
	size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
    
	for (size_type j=min_col; j<=max_col; j++) {
	    E[j]   = comp_E_entry(al,bl,i,j);
	    F      = comp_F_entry(al,bl,i,j);
	    M(i,j) = comp_M_entry(al,bl,i,j);
	}
    }
}

//===========================================================================
// Get the values in D from M for some al,bl
//

// compute the entries in the D matrix that
// can be computed from the matrix/matrices M
// for the subproblem al,bl,max_ar,max_br
//
// pre: M matrix is computed by a call to 
//      align_inside_arcmatch(al,max_ar,bl,max_br)
//
void AlignerP::fill_D(size_type al, size_type bl,
		      size_type max_ar, size_type max_br) {
    
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();
	arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();
	
	size_type ar = arcA.right();
	size_type br = arcB.right();
	
	//
	// if right ends ar,br exceed the limits max_ar,max_br resp.
	// the corresponding D entry is set to 0
	// in order to dissalow the arc match completely.
	// This occurs only due to an am heuristic.
	if (ar>max_ar || br>max_br) {	    
	    D(am) = (pf_score_t)0;
	} else {
	    D(am) = M(ar-1, br-1) * scoring->exp_arcmatch(am);
	}
    }
}

//===========================================================================
// Compute all entries of D

void AlignerP::align_D() {
    // ------------------------------------------------------------
    // General workflow:
    //
    // for all combinations of left arc ends al and ar
    // 1.) determie for which arc-pairs the D entries can be computed   
    // in one run, 2.) call align_inside_arcmatch 3.) call fill_D
    // ------------------------------------------------------------
    
    // ------------------------------------------------------------
    // traverse the left ends al,bl of arcs in descending order
    //
    for (size_type al=r.get_endA(); al>=r.get_startA(); al--) {
	
	size_type min_bl=r.get_startB();
	size_type max_bl=r.get_endB();
	
	// restrict range for left ends of bl due to trace controller
	min_bl = std::max(min_bl,params->trace_controller.min_col(al));
	max_bl = std::min(max_bl,params->trace_controller.max_col(al));
    
	
	for (size_type bl=max_bl; bl>=min_bl; bl--) {
	    
	    
	    // ------------------------------------------------------------
	    // get maximal right ends of possible arc matchs with left ends al,bl 
	    //
	    
	    size_type max_ar=al;
	    size_type max_br=bl;
	    
	    // ------------------------------------------------------------
	    // get the maximal right ends of any arc match with left ends (al,bl)
	    arc_matches.get_max_right_ends(al,bl,&max_ar,&max_br,false); 
	    
	    // ------------------------------------------------------------
	    // align under the maximal pair of arcs
	    //
	    align_inside_arcmatch(al, max_ar, bl, max_br);
	    
	    // ------------------------------------------------------------
	    // fill D matrix entries
	    //
	    fill_D(al, bl, max_ar, max_br);
	}
    }
    
    D_created=true; // now the matrix D is built up
}

//===========================================================================
// Do the complete inside phase of the partition function computation
//

pf_score_t AlignerP::align_inside() {
    // ------------------------------------------------------------
    // computes D matrix and then
    // does the alignment on the top level
    // ------------------------------------------------------------

    if (!D_created) {// std::cout<< "D is going to be created "<<endl; 
	alloc_inside_matrices();
      
	align_D();
    }

    assert(r.get_startA()>0);
    assert(r.get_startB()>0);

    align_inside_arcmatch(r.get_startA()-1, r.get_endA()+1, r.get_startB()-1,r.get_endB()+1);

    partFunc = M(r.get_endA(), r.get_endB());

    //assert(partFunc>0);

    return partFunc;
}


// ================================================================================
// REVERSE INSIDE
// ================================================================================

//#define printMrev(x,i,j) std::cout << x << " " << "Mrev(" << i << "," << j << ") := " << Mrev(i,j) << std::endl


// =================================================================
// init alignment matrix for aligning the fragments seqA(al..ar) and seqB(bl..br) 
void
AlignerP::init_Mrev(size_type al, size_type ar, size_type bl, size_type br) {
    // std::cout << "init_Mrev " << al << " " << ar << " " << bl << " " << br << std::endl;
    
    assert(al>=1);
    assert(bl>=1);
    
    Mrev(ar,br)=((pf_score_t)1)/pf_scale; // empty sequences
    // printMrev(1,ar,br);

    // initialize column br, subsequence B empty
    pf_score_t indel_score;
    indel_score = scoring->exp_indel_opening()/pf_scale;
    size_type i;
    for (i=ar; i>=al; ) { i--;
	if (params->trace_controller.max_col(i)<br) {
	    ++i;
	    break; // fill only as long as column bl is accessible
	}
	indel_score *= scoring->exp_gapA(i+1,br);
	Mrev(i,br) = indel_score;
	// printMrev(2,i,br);
    }
    // fill entries right of valid entries
    for ( ; i>=al; ) { i--;
	Mrev(i,params->trace_controller.max_col(i)+1) = 0;
	// printMrev(3,i,params->trace_controller.max_col(i)+1);
    }
    
    // initialize row ar, subsequence A empty
    indel_score = scoring->exp_indel_opening()/pf_scale;
    
    size_type min_col = std::max( bl-1, params->trace_controller.min_col(ar) );
    size_type j;
    for (j=br; j>min_col; ) { j--;
	indel_score *= scoring->exp_gapB(ar,j+1);
	Mrev(ar,j)= indel_score;
	// printMrev(4,ar,j);
    }
    // fill entries below valid entries
    for (size_type i=ar; i>=al; ) { i--;
	for (; j>std::max(bl-1,params->trace_controller.min_col(i)); ) { --j;
	    Mrev(i+1,j)=0;
	    // printMrev(5,i+1,j);
	}
    }
}

void
AlignerP::init_Erev(size_type al, size_type ar, size_type bl, size_type br) {
    for (size_type j=br; j>=bl; ) { --j;
	Erev[j]= (pf_score_t)0;
    }
}

inline
pf_score_t
AlignerP::comp_Erev_entry(size_type i, size_type j) {
    return 
	Erev[j] * scoring->exp_gapA(i+1, j)
	+
	(Mrev(i+1,j)-Erev[j]) * scoring->exp_gapA(i+1, j) * scoring->exp_indel_opening();
}

inline
pf_score_t
AlignerP::comp_Frev_entry(size_type i, size_type j) {
    return 
	Frev * scoring->exp_gapB(i, j+1)
	+
	(Mrev(i,j+1)-Frev) * scoring->exp_gapB(i, j+1) * scoring->exp_indel_opening();
}

// compute reversed M matrix entry; cases: base match, base in/del, arc match
// compute pf of alignments i+1..ar and j+1..br 
inline
pf_score_t
AlignerP::comp_Mrev_entry(size_type i, size_type j,size_type ar, size_type br) {
    
    pf_score_t pf;
    
    pf =
	// base match
	Mrev(i+1,j+1) * scoring->exp_basematch(i+1,j+1)
	
	// base del
	+ Erev[j]
	
	// base ins
	+ Frev;
    
    // arc match
    // standard case for arc match (without restriction to lonely pairs)
    const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(i+1);
    const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(j+1);

    // for all pairs of arcs in A and B that have right ends i+1 and j+1, respectively
    //
    for (BasePairs::LeftAdjList::const_iterator arcA=adjlA.begin();
	 arcA!=adjlA.end() && arcA->right() <= ar; ++arcA) {
	for (BasePairs::LeftAdjList::const_iterator arcB=adjlB.begin();
	     arcB!=adjlB.end() && arcB->right() <= br; ++arcB) {
	    
	    pf +=
		D(*arcA,*arcB) * Mrev(arcA->right(),arcB->right()) * pf_scale;
	}
    }
    return pf;
}


void
AlignerP::align_reverse(size_type al, size_type ar, size_type bl, size_type br, bool copy){
    assert(al>0);
    assert(bl>0);

    // Mrev.fill with -1 for debugging!
    for (size_t i=al; i<=ar; ++i) {
	for (size_t j=bl; j<=br; ++j) {
	    Mrev(i,j) = -1;
	}
    }



    init_Mrev(al,ar,bl,br);
    init_Erev(al,ar,bl,br);
    
    for(size_type i=ar; i>=al; ) { --i;//i from ar-1 downto al-1!
	Frev = (pf_score_t)0;
	
	// limit entries due to trace controller
	size_type min_col = std::max(bl,params->trace_controller.min_col(i)+1)-1;
	size_type max_col = std::min(br,params->trace_controller.max_col(i)+1)-1;
	
	for(size_type j=max_col+1; j>min_col; ) { --j;  //j from max_col downto min_col
	    Erev[j]   = comp_Erev_entry(i,j);
	    Frev      = comp_Frev_entry(i,j);
	    if (copy) { 
		Erev_mat(i,j) = Erev[j];
		Frev_mat(i,j) = Frev;
	    }

	    Mrev(i,j) = comp_Mrev_entry(i,j,ar,br);
	}
    }

    
    // cout << "Mrev " << al << " " << ar << " " << bl << " " << br <<std::endl;
    // for (size_t i=al; i<=ar; ++i) {
    // 	for (size_t j=bl; j<=br; ++j) {
    // 	    std::cout << Mrev(i,j) << " ";
    // 	}
    // 	std::cout << std::endl;
    // }
}


// ====================================================================================================
// ====================================================================================================
// OUTSIDE ALGORITHM
// ====================================================================================================
// ====================================================================================================


// helper functions for optimization

AlignerP::size_type
AlignerP::leftmost_covering_arc(size_type start,const BasePairs &bps,size_type l,size_type r) const {
    for (size_type i=start; i<l; ++i) {
	const BasePairs::LeftAdjList &adjl = bps.left_adjlist(i);
	for (BasePairs::LeftAdjList::const_iterator arc=adjl.begin(); arc!=adjl.end(); ++arc) {
	    if ( arc->right() > r ) {
		return i;
	    }
	}
    }
    return l;
}

AlignerP::size_pair
AlignerP::leftmost_covering_arcmatch(size_type al,size_type bl,size_type ar,size_type br) const {
    // implement in a fast but possibly under-estimating way
    
    size_pair sp(leftmost_covering_arc(r.get_startA(),bpsA,al,ar),
		 leftmost_covering_arc(r.get_startB(),bpsB,bl,br)
		 );
    return sp;
}


AlignerP::size_type
AlignerP::rightmost_covering_arc(const BasePairs &bps,size_type l,size_type r,size_type stop) const {
    assert(r>=0);
  
    for (size_type i=stop+1; i>r+1; ) { --i;
	const BasePairs::RightAdjList &adjl = bps.right_adjlist(i);
	for (BasePairs::RightAdjList::const_iterator arc=adjl.begin(); arc!=adjl.end(); ++arc) {
	    if ( arc->left() < l ) {
		return i;
	    }
	}
    }
    return r;
}

AlignerP::size_pair
AlignerP::rightmost_covering_arcmatch(size_type al,size_type bl,size_type ar,size_type br) const {
    // implement in a fast but possibly under-estimating way
    size_pair sp(rightmost_covering_arc(bpsA,al,ar,r.get_endA()),
		 rightmost_covering_arc(bpsB,bl,br,r.get_endB())
		 );    
    return sp;
}


//===========================================================================
// Initialization of outside matrices Mprime, Eprime


// // Initialize the last row and column of Mprime for aligning outside of (a,b).
// // Initialize for the computation of outside alignments for all holes that start with
// // left ends al,bl and end with right ends ar<=i<=lenA and br<=j<=lenB
// // The initialization of Mprime relies on the M matrix, this matrix has to be filled
// // with prefix alignment scores of the complete sequences.
// void
// AlignerP::init_Mprime(size_type al, size_type ar, size_type bl, size_type br) {   
//     assert(al>0 && bl>0);
    
//     Mprime(r.get_endA(),r.get_endB()) = M(al-1,bl-1);
//     //std::cout<<"in init_Mprime :  with al,ar,bl,br "<<al<<" "<<ar<< " : "<<bl<<" "<<br<<endl;
    
//     pf_score_t pf;
    
//     pf = M(al-1,bl-1) * scoring->exp_indel_opening();
//     for (int i=(int)r.get_endA()-1;i>=(int)ar; i--) {
// 	pf *= scoring->exp_gapA(i+1,r.get_endB());
// 	Mprime(i,r.get_endB()) = pf;
//     }
    
//     pf = M(al-1,bl-1) * scoring->exp_indel_opening();
//     for (int j=(int)r.get_endB()-1; j>=(int)br; j--) {
// 	pf *= scoring->exp_gapB(r.get_endA(),j+1);
// 	Mprime(r.get_endA(),j) = pf; 
//     }
// }

// inline
// void
// AlignerP::init_Eprime(size_type al, size_type ar, size_type bl, size_type br) {
//     for (int j=(int)r.get_endB()-1; j>=(int)br; j--) {
// 	Eprime[j] = 0;
//     }
// }


// compute a single entry of Eprime
inline 
pf_score_t 
AlignerP::comp_Eprime_entry(size_type al, size_type bl, size_type i, size_type j) {
    
    return 
	Eprime[j] * scoring->exp_gapA(i+1, j)
	+
	(Mprime(i+1,j)-Eprime[j]) * scoring->exp_gapA(i+1, j) * scoring->exp_indel_opening();
}


// compute a single entry of Fprime
inline
pf_score_t 
AlignerP::comp_Fprime_entry(size_type al, size_type bl, size_type i, size_type j) {
    return 
	Fprime * scoring->exp_gapB(i, j+1)
	+
	(Mprime(i,j+1)-Fprime) * scoring->exp_gapB(i, j+1) * scoring->exp_indel_opening();
}


pf_score_t 
AlignerP::virtual_Mprime(size_type al, size_type bl, size_type i, size_type j, size_type max_ar, size_type max_br) const {
    if (i>=max_ar || j>=max_br) {
	return M(al-1,bl-1)*Mrev(i,j)*pf_scale;
    }
    return Mprime(i,j);
}


// compute a single entry in Mprime
// standard cases in alignment: base match, base in/del, arc match
// computes partition function for alignments with hole al..i and bl..j (recursing to larger holes)
//
// pre: preceeding values in Mprime, Eprime, Fprime are computed
inline
pf_score_t 
AlignerP::comp_Mprime_entry(size_type al, size_type bl, size_type i, size_type j, size_type max_ar, size_type max_br) {
    
    //assert(params->trace_controller.is_valid(i,j));
    
    pf_score_t pf;
    
    // check proper intialization
    //assert(params->trace_controller.is_valid(i+1,j+1) || Mprime(i+1,j+1)==0 );
    //assert(params->trace_controller.is_valid(i+1,j  ) || Eprime[j]==0 );
    //assert(params->trace_controller.is_valid(i  ,j+1) || Fprime==0 );
    
    pf =
	// base match
	Mprime(i+1, j+1) * scoring->exp_basematch(i+1, j+1)
	
	// base del
	+ Eprime[j]
	
	// base ins
	+ Fprime;
    
    //std::cout<<"Max score of outside up to case 3: " << pf <<"  "<<al<<"  "<<bl<<"  "<<i<<"  "<<j<<endl;

    // arc match, case 4
    {
	const BasePairs::RightAdjList &adjlA = bpsA.right_adjlist(i+1);
	const BasePairs::RightAdjList &adjlB = bpsB.right_adjlist(j+1);

	// for all pairs of arcs in A and B that have right ends i+1 and j+1, respectively
	//

	for (BasePairs::RightAdjList::const_reverse_iterator arcA=adjlA.rbegin();
	     arcA!=adjlA.rend() && arcA->left() < al; ++arcA) {
	    for (BasePairs::RightAdjList::const_reverse_iterator arcB=adjlB.rbegin();
		 arcB!=adjlB.rend() && arcB->left() < bl; ++arcB) {
		// consider score for match of basepair
		
		// assert(Mrev(arcA->left(),arcB->left()) > 0);
		
		pf += Dprime(*arcA,*arcB) * Mrev(arcA->left(),arcB->left()) * pf_scale;
	    }
	}
	//std::cout<<"Max score of outside up to case 4: " << pf <<"  "<<al<<"  "<<bl<<"  "<<i<<"  "<<j<<endl;
    }
    
    // arc match, case 5
    {
		
	const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(i+1);
	const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(j+1);
	
	// for all pairs of arcs in A and B that have left ends i+1 and j+1, respectively
	for (BasePairs::LeftAdjList::const_iterator arcA=adjlA.begin(); arcA!=adjlA.end(); ++arcA) {
	    for (BasePairs::LeftAdjList::const_iterator arcB=adjlB.begin(); arcB!=adjlB.end(); ++arcB) {
		// consider score for match of basepairs
		//std::cout << *arcA << "." << *arcB << std::endl;
		
		//NOTE: if arcA, arcB cannot be matched due to heuristics, then D(*arcA,*arcB) is 0.

		pf += virtual_Mprime(al, bl, arcA->right(),arcB->right(),max_ar,max_br) * D(*arcA,*arcB) * pf_scale;
	    }
	}
    }

    //std::cout<<"Max score of outside up to case 5: " <<"  "<<al<<"  "<<bl<<"  "<<i<<"  "<<j<<endl;
    return pf;
}


// ===========================================================================
//
// NOTE: There is further potential for optimization:
// The latter optimization will have best effect, if we compute only
// the necessary entries in Mprime, Eprime, Fprime.
//
void
AlignerP::align_outside_arcmatch(size_type al,size_type ar,size_type max_ar,size_type bl,size_type br,size_type max_br) {
    assert(al>0);
    assert(bl>0);
    
    // fill the Mrev matrix 
    size_pair start = leftmost_covering_arcmatch(al,bl,ar,br);
    
    // cout << "Outside "
    //  	 <<al<<".."<<ar<<" "
    //  	 <<bl<<".."<<br <<" "
    //  	 <<"Start " << start.first<<" "<<start.second <<" "
    //  	 << max_ar << " " << max_br << std::endl;
    
    align_reverse(start.first+1,al-1,start.second+1,bl-1);
        
    // fill the outside matrices Mprime,Eprime,Fprime
    //
    // init Mprime (al,ar,bl,br);
    // init Eprime (al,ar,bl,br);
    
    Mprime.fill(-1); // only for debugging
    
    // initialize the valid entries in column max_br and row max_ar
    // note that max_ar,max_br is not necessarily valid!
    if (params->trace_controller.is_valid(max_ar,max_br)) {
	Mprime(max_ar,max_br) = M(al-1,bl-1)*Mrev(max_ar,max_br)*pf_scale;
    }
        
    // fill column max_br
    size_type i;
    for(i=max_ar; i>ar; ) { i--;
	if (params->trace_controller.max_col(i) < max_br) { i++; break; }
	if (params->trace_controller.is_valid(i,max_br)) {
	    Mprime(i,max_br) = M(al-1,bl-1)*Mrev(i,max_br)*pf_scale;
	}
    }

    // fill invalid entries right of valid entries 
    for ( ; i>ar; ) {
	i--;
	if (params->trace_controller.max_col(i)+1 <= max_br) {
	    Mprime(i,params->trace_controller.max_col(i)+1) = 0;
	}
    }
    
    // fill row max_ar ( Mprime and Eprime )
    size_type j;
    size_type min_col = std::max(br,params->trace_controller.min_col(max_ar));
    size_type max_col = std::min(max_br-1,params->trace_controller.max_col(max_ar));
    for(j=max_col+1; j>min_col;) { j--;
	Eprime[j]        = M(al-1,bl-1)*Erev_mat(max_ar,j)*pf_scale;
	Mprime(max_ar,j) = M(al-1,bl-1)*Mrev(max_ar,j)*pf_scale;
    }
    // fill invalid entries below valid entries 
    for (size_type i=max_ar; i>ar; ) { i--;
	for (; j>max(bl,params->trace_controller.min_col(i)); ) {
	    --j;
	    Mprime(i+1,j)=0;
	    Eprime[j]=0;
	}
    }
    
    for(size_type i=max_ar; i>ar; ) {
	i--;
	
	if (params->trace_controller.is_valid(i,max_br)) {
	    Fprime = M(al-1,bl-1)*Frev_mat(i,max_br)*pf_scale;
	} else {
	    Fprime=0;
	}

	size_type min_col = std::max(br,params->trace_controller.min_col(i));
	size_type max_col = std::min(max_br-1,params->trace_controller.max_col(i));
    	
	for(size_type j=max_col+1; j>min_col;) {
	    j--;
	    
	    Fprime      = comp_Fprime_entry(al,bl,i,j);
	    Eprime[j]   = comp_Eprime_entry(al,bl,i,j);
	    Mprime(i,j) = comp_Mprime_entry(al,bl,i,j,max_ar,max_br);
	}
    }
  
    // cout << "Mprime "<<al<<".."<<ar<<" "<<bl<<".."<<br <<std::endl;
    // cout << Mprime<<std::endl;
}

//===========================================================================
// fill all entries in Dprime for some (al,bl)
//

// compute the entries in the Dprime matrix that
// can be computed from the matrix/matrices Mprime
// for the subproblem al,bl,min_ar,min_br
//
// pre: Mprime matrices are computed by a call to 
//      align_out_arcmatch(al,min_ar,bl,min_br)
//
void
AlignerP::fill_Dprime(size_type al, size_type bl,
		      size_type min_ar, size_type min_br,
		      size_type max_ar, size_type max_br)
{
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();
	arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it); 
	
	size_type ar = am.arcA().right();
	size_type br = am.arcB().right();
	
	//
	// if right ends ar,br exceed the limits min_ar,min_br resp.
	// the corresponding Dprime entry is set to zero
	// in order to dissalow the arc match completely.
	// This occurs only due to the am heuristic.
	//
	if (ar<min_ar || br<min_br ) {
	    Dprime(am) = (pf_score_t)0;
	} else {
	    //std::cout << "Lookup Mprime("<<ar<<","<<br<<")="<<Mprime(ar,br)<<std::endl;
	    //assert( Mprime(ar,br)>0 );
	    Dprime(am) = virtual_Mprime(al, bl, ar, br, max_ar, max_br) * scoring->exp_arcmatch(am);
	}
    }
}


//===========================================================================

// compute all entries Dprime
void AlignerP::align_Dprime() {
    // ------------------------------------------------------------
    // General workflow:
    //
    // for all combinations of left arc ends al and bl
    // 1.) determine for which arc-pairs the Dprime entries can be computed   
    // in one run,
    // 2.) call align_outside_arcmatch
    // 3.) call fill_Dprime
    // ------------------------------------------------------------

    // ------------------------------------------------------------
    // traverse the left ends al,bl of arcs in ascending order
    //
    for (size_type al=r.get_startA(); al<=r.get_endA(); al++) {
	
	// restrict range for left ends of bl due to trace controller
	size_type min_bl = std::max(r.get_startB(), params->trace_controller.min_col(al));
	size_type max_bl = std::min(r.get_endB(),   params->trace_controller.max_col(al));
    
	for (size_type bl=min_bl; bl<=max_bl; bl++) {
	    
	    // ------------------------------------------------------------
	    // get minimal right ends of arc matchs with left ends al,bl 
	    //
	    size_type min_ar=r.get_endA()+1;
	    size_type min_br=r.get_endB()+1;
	    
	    arc_matches.get_min_right_ends(al,bl,&min_ar,&min_br); 
	    
	    // continue, when there is no arc match with left ends al,bl
	    // this is only a small optimization and not needed for correctness
	    if (min_ar > r.get_endA() || min_br > r.get_endB()) continue;
	    
	    // ------------------------------------------------------------
	    // get rightmost end of covering arc match.
	    // idea: for positions right of right_end, there is no dependency
	    // between the alignment of the fragments left and right
	    // of the hole.
	    //
	    size_pair max_r = rightmost_covering_arcmatch(al,bl,min_ar,min_br);
	    
	    // ------------------------------------------------------------
	    // align outside the arc
	    align_outside_arcmatch(al, min_ar, max_r.first, bl, min_br, max_r.second);
      
	    // ------------------------------------------------------------
	    // fill Dprime matrix entries
	    //
	    fill_Dprime(al, bl, min_ar, min_br, max_r.first, max_r.second);
	}
    }
    Dprime_created=true; // now the matrix Dprime is built up
}


//===========================================================================

// perform outside algo, compute Dprime
void AlignerP::align_outside() {
    if (!Dprime_created) { // computes Dprime matrix, only if not already created
    
	// pre-compute matrix Mrev00 that contains scores of all suffix alignments
	// this matrix is used in an optimization that
	// works in case that there is no arc-match that covers the whole 
	// and then just decomposes the fragment into prefix-alignment + suffix-alignment
	//
	// as a special trick, we can use the same 2D-matrix Mrev here,
	// as for the reverse matrix that we compute for each al,bl
	// argument: al,bl are visited in ascending order, for the
	// smaller Mrev one overwrites only values that are smaller or equal al,bl
	// for the optimization we need only entries that are larger al,bl

	alloc_outside_matrices();

	align_reverse(r.get_startA(),r.get_endA(),r.get_startB(),r.get_endB(),true);

	align_Dprime();
    }
}



// ====================================================================================================
// ====================================================================================================
// COMPUTING PROBABILITIES
// ====================================================================================================
// ====================================================================================================


//===========================================================================
//compute arc match probabilities

void 
AlignerP::compute_arcmatch_probabilities() {

    //cout.precision(2);
    
    /*
    cout << "D" << std::endl
     	 << Dmat << std::endl;
    
    cout << "Dprime" << std::endl
     	 << Dmatprime << std::endl;
    */
    
    // iterate over all arc matches
    for(ArcMatches::const_iterator it=arc_matches.begin(); arc_matches.end()!=it; ++it) {
	const Arc &arcA=it->arcA();
	const Arc &arcB=it->arcB();
	
	// trace_controller validity due to ArcMatches class
	assert(params->trace_controller.is_valid_match(arcA.left(),arcB.left()));
	assert(params->trace_controller.is_valid_match(arcA.right(),arcB.right()));
	
	am_prob(arcA.idx(),arcB.idx()) =
	    (D(arcA,arcB)/(long double)partFunc)
	    *  Dprime(arcA,arcB) * pf_scale / scoring->exp_arcmatch(*it);
	
	//std::cout << arcA << " " << arcB << ": " << D(arcA,arcB) << " " << Dprime(arcA,arcB) << " " <<  am_prob(arcA.idx(),arcB.idx()) <<  std::endl;  
	
	if (! (am_prob(arcA.idx(),arcB.idx())<=1) ) {
	    std::cerr << "ERROR: am prob " << arcA <<" " << arcB <<" " << am_prob(arcA.idx(),arcB.idx())<<std::endl;
	    exit(-1);
	}
    }
    //std::cout<<"Arc match probs calculated"<<endl;
    //std::cout<<am_prob<<std::endl;
}

//===========================================================================
// compute base match probabilities

// pre: arc match probabilites am_prob are already computed
void
AlignerP::compute_basematch_probabilities( bool basematch_probs_include_arcmatch )
{
    //std::cout<<"Skip compute_basematch_probabilities"<<std::endl;
    //return;

    // sum the conditional partition functions for all distinct cases
    
    // consider only arc matchs with a probability of at more than am_prob_threshold
    double am_prob_threshold=sqrt(params->min_am_prob); // use something quite conservative as threshold, such that user can still control this 
	
    // --------------------------------------------------
    // cases, where edge is enclosed by arc match
    //
    
    for(size_type al=r.get_startA();al<=r.get_endA();al++){

	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(al));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(al));
	for(size_type bl=min_col; bl<=max_col; bl++){
	    
	    // trace controller allows trace through (al,bl), but not
	    // necessarily match of al and bl
	    if (! params->trace_controller.is_valid_match(al,bl)) continue;
	    
	    const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(al);
	    const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(bl);

	    if(adjlA.size() >0 && adjlB.size()>0)
		{
		    
		    assert(D_created);assert(Dprime_created);
		    
		    // get max_ar and max_br, where am_prob larger than threshold
		    // (which implies that the arc match is valid!).
		    // This is used only for limiting the inside recomputation.
		    
		    size_type max_ar=al;
		    size_type max_br=bl;
		    
		    for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
			 arcA!=adjlA.end(); ++arcA) {
			for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
			     arcB!=adjlB.end(); ++arcB) {
			    size_type ar = arcA->right();
			    size_type br = arcB->right();
			    
			    // Note that the match of arcA and arcB
			    // may be illegal due to heuristics!
			    // However, in this case am_prob is 0.0,
			    // since am_prob is of type SparseMatrix
			    // with default 0.0 and we wrote values
			    // only for arc matches in the arc_matches
			    // object (see compute_arcmatch_probabilities).
			    
			    if ( am_prob(arcA->idx(),arcB->idx()) > am_prob_threshold ) {
				max_ar=max(max_ar, ar);
				max_br=max(max_br, br);
			    }
			}
		    }
		    
		    // Align inside limited by the determined maximal ar and br
		    align_inside_arcmatch(al,max_ar,bl,max_br);
		    
		    for (BasePairs::LeftAdjList::const_iterator arcA=adjlA.begin();
			 arcA!=adjlA.end(); ++arcA) {
			for (BasePairs::LeftAdjList::const_iterator arcB=adjlB.begin();
			     arcB!=adjlB.end(); ++arcB) {
			    
			    if (am_prob(arcA->idx(),arcB->idx()) > am_prob_threshold) {
				// again note that the above comparison is sufficient to guarantee the validity of
				// the arc match arcA~arcB
				
				size_type ar=arcA->right();
				size_type br=arcB->right();
				
				// compute the reverse matrix for all values below of the arc match (al,ar)~(bl,br)
				align_reverse(al+1,ar-1,bl+1,br-1);
				
				// a part of the pf-contrib can be computed outside of the loops
				pf_score_t arcmatch_outside_pf=
				    Dprime(*arcA,*arcB);
				
				// add contributions for all alignment edges enclosed by the arc match (arcA,arcB)
				for(size_type i=al+1;i<ar;i++){
				    
				    // limit entries due to trace controller
				    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
				    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
				    
				    for(size_type j=min_col;j<=max_col;j++){
					
					if ( ! params->trace_controller.is_valid_match(i,j) ) continue;
										
					bm_prob(i,j) += 
					    M(i-1,j-1)
					    * scoring->exp_basematch(i,j)
					    * Mrev(i,j)
					    * pf_scale
					    * arcmatch_outside_pf
					    * pf_scale;
				    }
				}
			    }
			}
		    }
		}
	    //	else{std::cout<<"left adj list empty " <<al<<"  "<<bl<<endl;}
	}
    }

    // --------------------------------------------------
    // extra case, where there is no enclosing arc match of alignment edge (i,j)
  
    align_inside_arcmatch(0,r.get_endA()+1,0,r.get_endB()+1);
    align_reverse(r.get_startA(),r.get_endA(),r.get_startB(),r.get_endB());
  
    for(size_type i=r.get_startA();i<=r.get_endA();i++){
	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(i));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(i));
	for(size_type j=min_col;j<=max_col;j++){
	    
	    if ( ! params->trace_controller.is_valid_match(i,j) ) continue;
	    
	    bm_prob(i,j) += M(i-1,j-1) * scoring->exp_basematch(i,j) * Mrev(i,j) * pf_scale;
	}
    }
  

    // --------------------------------------------------
    // divide the conditional partition functions by total partition function  
	
    for(size_type i=r.get_startA();i<=r.get_endA();i++){
	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(i));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(i));
	for(size_type j=min_col;j<=max_col;j++){
	    
	    if ( ! params->trace_controller.is_valid_match(i,j) ) continue;
	    
	    bm_prob(i,j)=bm_prob(i,j)/partFunc;
	    
	    //assert(bm_prob(i,j)<=1);
	    if (bm_prob(i,j)>1) {
		std::cerr << "ERROR: bm prob " << i <<" " << j <<" " << bm_prob(i,j)<<std::endl;
		exit(-1);
	    }
	}
    }

    // --------------------------------------------------
  
    if ( basematch_probs_include_arcmatch ) {
	// in this mode the base match probs cover the cases where
	// the bases are matched due to a structural match.
	// thus, we add these probabilities.
      
	// iterate over all arc matches
	for(ArcMatches::const_iterator it=arc_matches.begin(); arc_matches.end()!=it; ++it) {
	    const Arc &arcA=it->arcA();
	    const Arc &arcB=it->arcB();
	  
	    double amp = am_prob(arcA.idx(),arcB.idx());
	  
	    bm_prob(arcA.left(),arcB.left()) += amp;
	    bm_prob(arcA.right(),arcB.right()) += amp;
	}
    }

}



//===========================================================================
// write base match probabilities

void 
AlignerP::write_basematch_probabilities(std::ostream &out)
{
    for(size_type i=1;i<=r.get_endA();i++){
	for(size_type j=1;j<=r.get_endB();j++){
	    if (bm_prob(i,j)>=params->min_bm_prob) {
		out<<i<<" "<<j<<" "<<bm_prob(i,j);
		out<<std::endl;
	    }
	}
    }
}

//===========================================================================
// write arcmatch probabilities
//
void
AlignerP::write_arcmatch_probabilities(std::ostream &out)
{
    // iterate over all arc matches
    for(ArcMatches::const_iterator it=arc_matches.begin(); arc_matches.end()!=it; ++it) {
	const Arc &arcA=it->arcA();
	const Arc &arcB=it->arcB();
	
	if (am_prob(arcA.idx(),arcB.idx())>=params->min_am_prob) {
	    out <<arcA.left()<<" "<<arcA.right()<<" "
		<<arcB.left()<<" "<<arcB.right()<<" "
		<<am_prob(arcA.idx(),arcB.idx())
		<<std::endl;
	}
    }
}



//===========================================================================
// fragment match probabilities
//
double
AlignerP::compute_fragment_match_prob(size_type i,size_type j,size_type k,size_type l) {

    // computes the probability to match [i..j] and [k..l] 
    // as D_{(i-1,j+1);(k-1,l+1)}/arcmatch(i-1,j+1,k-1,l+1) * D'_{(i,j);(k,l)}/arcmatch(i,j,k,l) / Z
    
    // in general these values are not computed yet, so we will just recompute them
    // possible optimization: see whether values exist already in D or D'
    
    pf_score_t in;  // pf inside  of fragments [i..j] and [k..l]
    pf_score_t out; // pf outside of fragments [i..j] and [k..l]
    
    
    M.fill(0);
    align_inside_arcmatch(i-1,j+1,k-1,l+1); // arcs (i-1,j+1) and (k-1,l+1) enclose the fragments
    in = M(j,l);
    
    // ensure that pre-conditions are met for align_outside_arcmatch 
    align_inside_arcmatch(r.get_startA()-1, r.get_endA()+1, r.get_startB()-1,r.get_endB()+1);
    align_reverse(r.get_startA(),r.get_endA(),r.get_startB(),r.get_endB(),true);
    
    
    size_pair max_r = rightmost_covering_arcmatch(i,k,j,l);
    
    //Mprime.fill(0);
    
    align_outside_arcmatch(i, j, max_r.first, k, l, max_r.second);
    
    out = virtual_Mprime(i, k, j, l, max_r.first, max_r.second);
    
    // std::cout << "in: "<<in<<" out: "<<out<<std::endl;
    
    return in/((long double)partFunc)*out*pf_scale;
    
}
