#include "aligner_impl.hh"


#include "anchor_constraints.hh"
#include "trace_controller.hh"
#include "basepairs.hh"
#include "sequence.hh"

#include <math.h>
#include <assert.h>

#include <queue>

#include <iostream>


namespace LocARNA {

    /*
      SEMANTIC OF ANCHOR CONSTRAINTS
  
      Anchor constraints (i,j) enforce that positions i in A and j in B are matched;
      neither i nor j are deleted (for local alignment, this implies that
      both positions occur in the local alignment)
  
      Names that occur in only one sequence, do not impose constraints
    */



    /* SEMANTICS OF MAX DIFF HEURISTIC

       restrict matrix cells (i,j) to valid trace cells due to trace_controller
   
    */


    /*
      recursions for global alignment:

      the following recursions simplify the actual locarna algorithm
      most notably:
      * the computation of matrices M^ab and M^a'b' where al=a'l and bl=bl' is joined
      * no constraint handling is done
      * no heuristics (max-diff-am)

      standard recursion
      M^ab_ij = max (M^ab_i-1j-1 + s(i,j);
      M^ab_i-1j + g;
      M^ab_ij-1 + g;
      max M^ab_a'l-1b'l-1 + D_a'b': a'b' a'r=i b'r=j)
      D_ab    = M^ab_ar-1br-1 + s(a,b)


      no lonely pairs
      M^ab_ij = max (M^ab_i-1j-1 + s(i,j);
      M^ab_i-1j + g;
      M^ab_ij-1 + g;
      max M^ab_a'l-1b'l-1 + D_a'b': a'b' a'r=i b'r=j)
      D_ab    = max { M^a'b'_a'r-1b'r-1 + s(a,b) + s(a',b');
      D_a'b' + s(a,b) }
      where a'l-1=al,b'l-1=bl,a'r+1=ar,b'r+1=br


      stacking
      M^ab_ij = max ( M^ab_i-1j-1 + s(i,j);
      M^ab_i-1j + g;
      M^ab_ij-1 + g;
      max M^ab_a'l-1b'l-1 + D_a'b': a'b' a'r=i b'r=j )
      D_ab    = max ( M^ab_ar-1br-1 + s(a,b);
      D^a'b' + s'(a',b') )
      where a'l-1=al,b'l-1=bl,a'r+1=ar,b'r+1=br
    */





    // ------------------------------------------------------------
    // Aligner: align / compute similarity
    //

    Aligner::Aligner(const Aligner &aligner)
	:pimpl_(new AlignerImpl(*aligner.pimpl_))
    {}
    

    AlignerImpl::AlignerImpl(const AlignerImpl &a)
	: params_(a.params_),
	  scoring_(a.scoring_),
	  mod_scoring_(0),
	  seqA_(a.seqA_),
	  seqB_(a.seqB_),
	  arc_matches_(a.arc_matches_),
	  bpsA_(a.bpsA_),
	  bpsB_(a.bpsB_),
	  r_(a.r_),
	  Dmat_(a.Dmat_),
	  Ms_(a.Ms_),
	  Es_(a.Es_),
	  Fs_(a.Fs_),
	  min_i_(a.min_i_),
	  min_j_(a.min_j_),
	  max_i_(a.max_i_),
	  max_j_(a.max_j_),
	  D_created_(a.D_created_),
	  alignment_(a.alignment_),
	  def_scoring_view_(this),
	  mod_scoring_view_(this),
	  free_endgaps_(a.free_endgaps_)
    {}
    
    Aligner::Aligner(const AlignerParams &ap) 
	: pimpl_(new AlignerImpl(*(ap.seqA_),*(ap.seqB_),*(ap.arc_matches_),&ap,ap.scoring_))
    {} 
 

    AlignerImpl::AlignerImpl(const Sequence &seqA, 
			     const Sequence &seqB,
			     const ArcMatches &arc_matches,
			     const AlignerParams *ap,
			     const Scoring *s
			     )
	: params_(ap),
	  scoring_(s),
	  mod_scoring_(0),
	  seqA_(seqA), seqB_(seqB),
	  arc_matches_(arc_matches),
	  bpsA_(arc_matches_.get_base_pairsA()),
	  bpsB_(arc_matches_.get_base_pairsB()),
	  r_(1,1,seqA.length(),seqB.length()),
	  min_i_(1),
	  min_j_(1),
	  max_i_(seqA.length()),
	  max_j_(seqB.length()),
	  D_created_(false),
	  alignment_(seqA,seqB),
	def_scoring_view_(this),
	mod_scoring_view_(this),
	free_endgaps_(params_->free_endgaps_)
    {
	Ms_.resize(params_->struct_local_?8:1);
	Es_.resize(params_->struct_local_?4:1);
	Fs_.resize(params_->struct_local_?4:1);
    
	Dmat_.resize(bpsA_.num_bps(),bpsB_.num_bps());
	Dmat_.fill(infty_score_t::neg_infty);
    
	for (size_t k=0; k<(params_->struct_local_?8:1); k++) {
	    Ms_[k].resize(seqA_.length()+1,seqB_.length()+1);
	}
	for (size_t k=0; k<(params_->struct_local_?4:1); k++) {
	    Es_[k].resize(seqB_.length()+1);
	}
    }


    AlignerImpl::~AlignerImpl() {
    	if (mod_scoring_!=0) delete mod_scoring_;
    }

    Aligner::~Aligner() {
	delete(pimpl_);
    }

    
    Alignment const & 
    Aligner::get_alignment() const {return pimpl_->alignment_;} 

    void
    Aligner::set_alignment(const Alignment &alignment) {
	pimpl_->alignment_ = alignment;
    }
    

    // standard cases in alignment: base match, base in/del, arc match
    // (for structure local alignment this is extended by exclusion handling)
    // 
    // if lonely basepairs are disallowed, there is special treatment
    //
    // align_noex has a side effect: it computes entry in Es[state] and Fs[state]
    template<class ScoringView>
    infty_score_t
    AlignerImpl::align_noex(int state, pos_type al, pos_type bl, pos_type i, pos_type j,ScoringView sv) {
    
	assert(0<=state && state<4);
    
	assert(params_->trace_controller_->is_valid(i,j));
    
	M_matrix_t &M = Ms_[state];
	ScoreVector &E = Es_[state];
	infty_score_t &F = Fs_[state];
    
	// compute E entry
	if ( (! params_->constraints_->aligned_in_a(i)) ) {
	    // due to constraints, i can be deleted
	    E[j] = 
		std::max( E[j] + sv.scoring()->gapA(i),
			  M(i-1,j) + sv.scoring()->gapA(i) + sv.scoring()->indel_opening() );
	} else {
	    // due to constraints, i cannot be deleted
	    E[j] = infty_score_t::neg_infty;
	}
    
	// compute F entry
	if ( (! params_->constraints_->aligned_in_b(j)) ) {
	    // due to constraints, j can be inserted
	    F=std::max( F + sv.scoring()->gapB(j),
			M(i,j-1) + sv.scoring()->gapB(j) + sv.scoring()->indel_opening() );
	} else {
	    // due to constraints, j cannot be inserted
	    F = infty_score_t::neg_infty;
	}
    

	// use tainted type to save operations that normalize infinity
	tainted_infty_score_t max_score = infty_score_t::neg_infty;
    
	// base match
	if (params_->constraints_->allowed_edge(i,j)) {
	    max_score = M(i-1,j-1) + sv.scoring()->basematch(i,j);
	}
    
	// base del
	max_score=std::max(max_score, (tainted_infty_score_t)E[j]);
    
	// base ins
	max_score=std::max(max_score, (tainted_infty_score_t)F);

	// arc match
    
	// standard case for arc match (without restriction to lonely pairs)
	//    
    
	if ( params_->constraints_->allowed_edge(i,j) ) {
	    const BasePairs::RightAdjList &adjlA = bpsA_.right_adjlist(i);
	    const BasePairs::RightAdjList &adjlB = bpsB_.right_adjlist(j);
	
	    // for all pairs of arcs in A and B that have right ends i and j, respectively
	    //
	    for (BasePairs::RightAdjList::const_iterator arcA=adjlA.begin();
		 arcA!=adjlA.end() && arcA->left() > al  ; ++arcA) {
		for (BasePairs::RightAdjList::const_iterator arcB=adjlB.begin();
		     arcB!=adjlB.end() && arcB->left() > bl ; ++arcB) {
		
		    // no need to check (params_->constraints_->allowed_edge(arcA->left(),arcB->left()))
		    // or other "constraints"
		    // because for these arc matches holds that sv.D(*arcA,*arcB)==neg_infty
		
		
		    tainted_infty_score_t new_score =
			M(arcA->left()-1,arcB->left()-1)
			+ sv.D(*arcA,*arcB);
		
		    if (new_score > max_score) {
			//std::cout << *arcA << "-"<< *arcB << ": "<<M(arcA->left()-1,arcB->left()-1)<<"+"<<D(*arcA,*arcB)<<"="<<new_score<<std::endl;
			max_score=new_score;
		    }
		}
	    }
	}
    
	return max_score;

    
	// The following code turned out to be much slower than the above one

	//     const ArcMatchVec &right_adj_list = arc_matches.common_right_end_list(i,j);
    
	//     for(ArcMatchVec::const_iterator it=right_adj_list.begin(); right_adj_list.end() != it; ) {
	
	// 	// NOTES: *it is the arc match index
	// 	//        we iterate only over valid arc matches, i.e.
	// 	//        constraints (including anchor c. and heuristic ones) are satisified
	
	// 	const ArcMatch &am = *it;
	
	// 	const Arc &arcA=am.arcA();
	// 	const Arc &arcB=am.arcB();
	
	// 	//if ( arcA.left() <= al || arcB.left() <= bl ) {++it; continue;}
	
	
	// 	// These optimizations assume that the list is sorted
	// 	//  lexicographically descending by (arcA.left, arcB.left)
	// 	//
	// 	if ( arcA.left() <= al ) break;
	
	// 	if ( arcB.left() <= bl ) {
	    
	// 	    // iterate to the next different al.
	// 	    // this could be optimized further using a helper vector
	// 	    // that allows to jump directly to this entry
	// 	    do {
	// 		it++;
	// 	    } while (right_adj_list.end()!=it && it->arcA().left()==al);
	    
	// 	    continue;
	// 	}
	
	// 	//std::cerr << am.idx() << std::endl;
	
	// 	max_score = std::max( max_score, M(arcA.left()-1,arcB.left()-1) + D[am.idx()] );
	
	// 	++it;
	//     }
    
    }


    // generic initalization method.
    //
    // The method takes care of anchor constraints. Positions that are constraint cannot be deleted/inserted
    //
    // The method initializes Ms and Es. ATTENTION: Fs has to be re-initialized for each row i during recursion.  
    //
    // Init has to be aware of the trace controllers restriction that potentially limit the
    // computed entries in each matrix row
    // It is necessary to initialize all invalid entries that may be accessed from these valid entries with neg_infty
    //
    //
    // In case of global alignment without exclusions and without restriction due to trace controller initialization is easy: 
    //    M(al,bl)=0
    //    M(i,bl)=indel_opening+i*indel, for i>0
    //    M(al,j)=indel_opening+j*indel, for j>0
    //
    //    E[j]=-infty, (since there is no alignment, where subseq of a is empty that deletes the last position of the subseq of a)
    //    F is initialized with -infty for row i>0 
    //
    //

    // Initialization in case of local alignment 

    template <class ScoringView>
    void
    AlignerImpl::init_state(int state, pos_type al, pos_type ar, pos_type bl, pos_type br, 
			    bool globalA, bool exclA,
			    bool globalB, bool exclB, 
			    ScoringView sv) {
    
	// alignments that have empty subsequence in A (i=al) and
	// end with gap in alistr of B do not exist ==> -infty 
	if (state<4) {
	    ScoreVector   &E = Es_[state];
	    for (pos_type j=bl; j<br; j++) {
		E[j]=infty_score_t::neg_infty;
	    }
	}
    
	M_matrix_t &M = Ms_[state];
    
	// al,bl can only be reached in states, where this is legal with cost 0 for empty alignment
	M(al,bl) = (infty_score_t)0;
    
	// std::cout << "COL "<<bl<<" AL: "<<al<<" AR: "<<ar<<std::endl;
    
	// init first col bl
	//
	infty_score_t indel_score=(infty_score_t)sv.scoring()->indel_opening();
	if (exclA) {
	    indel_score = (infty_score_t)sv.scoring()->exclusion();
	} else if (!globalA) {
	    indel_score = (infty_score_t)0;
	}

	// handling of anchor constraints:
	// anchored positions must not be excluded, 
	// nor deleted
    
	pos_type i;
	for (i=al+1; i<ar; i++) {

	    if (params_->trace_controller_->min_col(i)>bl) break; // fill only as long as column bl is accessible

	    if (!indel_score.is_neg_infty()) {
		if (params_->constraints_->aligned_in_a(i)) {
		    indel_score=infty_score_t::neg_infty;
		}
		else if (!exclA && globalA) {
		    indel_score += sv.scoring()->gapA(i);
		}
	    }
	    M(i,bl) = (infty_score_t)indel_score;
	}

	// fill entries left of valid entries 
	for ( ; i<ar; i++) {
	    assert(params_->trace_controller_->min_col(i)>bl);
	    M(i,params_->trace_controller_->min_col(i)-1) = infty_score_t::neg_infty; 
	}
    
	// init first row al
	//
	indel_score=(infty_score_t)sv.scoring()->indel_opening();
	if (exclB) {
	    indel_score = (infty_score_t)sv.scoring()->exclusion();
	} else if (!globalB) {
	    indel_score = (infty_score_t)0;
	}
    
	pos_type j;
	for (j=bl+1 ; j < std::min(br, params_->trace_controller_->max_col(al)+1) ; j++) {
	    if (!indel_score.is_neg_infty()) {
		if (params_->constraints_->aligned_in_b(j)) {
		    indel_score=infty_score_t::neg_infty;
		}
		else if (!exclB && globalB && !indel_score.is_neg_infty()) {
		    indel_score += sv.scoring()->gapB(j);
		}
	    }
	    M(al,j) = (infty_score_t)indel_score;
	}

	// fill entries above valid entries 
	// here j points to one position right of the last initialized entry in row al
	for (i=al+1; i<ar; i++) {
	    for (;
		 j<std::min(br,params_->trace_controller_->max_col(i)+1); ++j) {
		M(i-1,j)=infty_score_t::neg_infty;
	    }
	}

    }

    // ----------------------------------------
    // recomputes M matrix/matrices
    // after the call the matrix is filled in the range [al..ar-1] x [bl..br-1]
    void
    AlignerImpl::align_in_arcmatch(pos_type al,pos_type ar,pos_type bl,pos_type br,
				   bool allow_exclusion) {
	
	assert(br>0); // if br<=0 we run into trouble below when computing br-1

	// cout << al << " " << ar <<" " << bl << " " << br <<endl;

	// When using RMatrix as type of the M matrices,
	// the M matrix/matrices can be restricted to the range of the arc match
	// for (pos_type state=0; state < ((allow_exclusion)?8:1); state++) {
	//std::cout <<state<<" "<<Ms_[state].sizes().first<<" "<<Ms_[state].sizes().second<<" "<<al<<" "<<ar<<" "<<bl<<" "<<br<<std::endl; 
	//Ms_[state].restrict(al,ar-1,bl,br-1);
	// }
    
    
	// if in a sequence the state is not open than gaps with cost scoring->gap() have to be introduced.
	// In open cases the maximal similarity is 0 in the first row/column,
	// since an exclusion can be introduced.
    
	init_state(E_NO_NO,al,ar,bl,br,true ,false,true ,false,def_scoring_view_);
    
	if (allow_exclusion) {
	    init_state(E_X_NO, al,ar,bl,br, true , true , true , false,def_scoring_view_);
	    init_state(E_NO_X, al,ar,bl,br, true , false, true , true ,def_scoring_view_);
	    init_state(E_X_X,  al,ar,bl,br, true , true , true , true ,def_scoring_view_);
	
	    // open states
	    init_state(E_OP_NO,al,ar,bl,br, false, false, true , false,def_scoring_view_);
	    init_state(E_NO_OP,al,ar,bl,br, true , false, false, false,def_scoring_view_);
	    init_state(E_X_OP, al,ar,bl,br, true , true , false, false,def_scoring_view_);
	    init_state(E_OP_X, al,ar,bl,br, false, false, true , true ,def_scoring_view_);
	}

	// ----------------------------------------
	// alignment for state E_NO_NO
	//
    
	for (pos_type i=al+1; i<ar; i++) {
	    Fs_[E_NO_NO]=infty_score_t::neg_infty;
	
	    // limit entries due to trace controller
	    pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
	    pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	
	    for (pos_type j=min_col; j<=max_col; j++) {
		Ms_[E_NO_NO](i,j)=align_noex(E_NO_NO,al,bl,i,j,def_scoring_view_);
	    }
	}
    
	//
	// end state E_NO_NO
	// ----------------------------------------
    
	if (allow_exclusion) {
	    int state;

	    state=E_OP_NO;
	    for (pos_type i=al+1; i<ar; i++) { 
	
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(
					      params_->constraints_->aligned_in_a(i)?infty_score_t::neg_infty:Ms_[state](i-1,j),
					      Ms_[E_NO_NO](i,j)
					      );
		}
	    }
	    
	    state=E_NO_OP;
	    for (pos_type i=al+1; i<ar; i++) {
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(
					      params_->constraints_->aligned_in_b(j)?infty_score_t::neg_infty:Ms_[state](i,j-1),
					      Ms_[E_NO_NO](i,j)
					      );
		}
	    }

	    state=E_NO_X;
	    for (pos_type i=al+1; i<ar; i++) {
		Fs_[state]=infty_score_t::neg_infty;
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(align_noex(state,al,bl,i,j,def_scoring_view_),
					      Ms_[E_NO_OP](i,j)+scoring_->exclusion());
		}
	    }
	
	    state=E_OP_X;
	    for (pos_type i=al+1; i<ar; i++) {
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(
					      params_->constraints_->aligned_in_a(i)?infty_score_t::neg_infty:Ms_[state](i-1,j),
					      Ms_[E_NO_X](i,j)
					      );
		}
	    }
	
	    state=E_X_NO;
	    for (pos_type i=al+1; i<ar; i++) {
		Fs_[state]=infty_score_t::neg_infty;
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(align_noex(state,al,bl,i,j,def_scoring_view_),
					      Ms_[E_OP_NO](i,j)+scoring_->exclusion());
		}
	    }

	    state=E_X_OP;
	    for (pos_type i=al+1; i<ar; i++) {
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = std::max(
					      params_->constraints_->aligned_in_b(j)?infty_score_t::neg_infty:Ms_[state](i,j-1),
					      Ms_[E_X_NO](i,j)
					      );
		}
	    }

	    state=E_X_X;
	    for (pos_type i=al+1; i<ar; i++) {
		Fs_[state]=infty_score_t::neg_infty;
		// limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params_->trace_controller_->min_col(i));
		pos_type max_col = std::min(br-1,params_->trace_controller_->max_col(i));
	    
		for (pos_type j=min_col; j<=max_col; j++) {
		    Ms_[state](i,j) = 
			std::max(align_noex(state,al,bl,i,j,def_scoring_view_),
				 std::max(Ms_[E_OP_X](i,j)+scoring_->exclusion(),
					  Ms_[E_X_OP](i,j)+scoring_->exclusion())); }
	    }
	}
    }

    // compute the entries in the D matrix that
    // can be computed from the matrix/matrices M
    // for the subproblem al,bl,max_ar,max_br
    //
    // pre: M matrices are computed by a call to 
    //      align_in_arcmatch(al,max_ar,bl,max_br,params_->struct_local_)
    //
    void 
    AlignerImpl::fill_D_entries(pos_type al, pos_type bl)
    {
	for(ArcMatchIdxVec::const_iterator it=arc_matches_.common_left_end_list(al,bl).begin();
	    arc_matches_.common_left_end_list(al,bl).end() != it; ++it ) {
	
	    const ArcMatch &am = arc_matches_.arcmatch(*it);
	
	    const Arc &arcA=am.arcA();
	    const Arc &arcB=am.arcB();

	
	    pos_type ar = arcA.right();
	    pos_type br = arcB.right();
	
	    infty_score_t m=Ms_[0](ar-1,br-1);
	    if (params_->struct_local_) {
		// if we align structure local
		// we need to determine the maximum of entries in Ms_[k],
		// where k is a non-open state
		for (size_t k=1; k<4; k++) {
		    m=std::max(m, Ms_[k](ar-1,br-1));
		}
	    }
	
	    D(am) = m + scoring_->arcmatch(am);
	
	    //std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << std::endl;

	    if (scoring_->stacking()) {
		if (arc_matches_.exists_inner_arc_match(am)
		    &&
		    scoring_->is_stackable_am(am)
		    ) {
		    const ArcMatch &inner_am = arc_matches_.inner_arc_match(am);
		
		    D(am) =
			std::max(D(am),
				 D(inner_am) + scoring_->arcmatch(am,true));
		}
	    }
	}
    }

    void 
    AlignerImpl::fill_D_entries_noLP(pos_type al, pos_type bl) {
	// get adj lists of arcs starting in al-1, bl-1
    
	for(ArcMatchIdxVec::const_iterator it=arc_matches_.common_left_end_list(al-1,bl-1).begin();
	    arc_matches_.common_left_end_list(al-1,bl-1).end() != it; ++it ) {
	
	    const ArcMatch &am = arc_matches_.arcmatch(*it);
	
	    pos_type ar = am.arcA().right()-1;
	    pos_type br = am.arcB().right()-1;

	    // only for arc matches which can occur paired with an inner arc match (no lonely pairs) 
	    // therefore check whether inner arc exists
	    // if stacking scores are used, the am has to be stackable to the inner arc,
	    // i.e. the joint probabilities have to be greater than 0
	    if (arc_matches_.exists_inner_arc_match(am)
		&&
		( ! scoring_->stacking() || scoring_->is_stackable_am(am) )
		) { 
		const ArcMatch& inner_am = arc_matches_.inner_arc_match(am);
	    
		infty_score_t m=Ms_[0](ar-1,br-1);
		if (params_->struct_local_) {
		    // if we align structure local
		    // we need to determine the maximum of entries in Ms_[k],
		    // where k is a non-open state
		    for (size_t k=1; k<4; k++) {
			m=std::max(m, Ms_[k](ar-1,br-1));
		    }
		}
	    
		//std::cout << scoring_->arcmatch(am,scoring_->stacking()) <<" "<< scoring_->arcmatch(inner_am) <<" "<< m << " ";
	    
		// score of outer arc match, where inner arc match is aligned too
		D(am) = 
		    std::max(m + scoring_->arcmatch(inner_am),
			     D(inner_am))
		    + scoring_->arcmatch(am,scoring_->stacking());
	    
		//std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << " ; ";
		//std::cout <<"["<< inner_am.arcA() << "," <<inner_am.arcB() <<"]:" <<D(inner_am)<< std::endl;

	    }
	}
    }



    // compute all entries D
    void
    AlignerImpl::align_D() {
	// ------------------------------------------------------------
	// General workflow:
	//
	// for all combinations of left arc ends al and ar
	// 1.) determine for which arc-pairs the D entries can be computed   
	// in one run, 2.) call align_in_arcmatch 3.) call fill_D_entries
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// traverse the left ends al,bl of arcs in descending order
	// (restrict by trace controller and r)
	//
	// for al in r_.get_endA() .. r_.get_startA
	for (pos_type al=r_.get_endA()+1; al>r_.get_startA(); ) { al--; 
	
	    pos_type max_bl = std::min(r_.get_endB(),params_->trace_controller_->max_col(al));
	    pos_type min_bl = std::max(r_.get_startB(),params_->trace_controller_->min_col(al));
	
	    // for bl in max_bl .. min_bl
	    for (pos_type bl=max_bl+1; bl > min_bl;) { bl--; 
	    
		if (! ( params_->constraints_->allowed_edge(al,bl)
			&& params_->trace_controller_->is_valid_match(al,bl) )
		    ) continue;
	    
		// ------------------------------------------------------------
		// get maximal right ends of arcs with left ends al,bl 
		// where max_diff_am conditions hold
		// and no_lonely_pairs condition holds
		//
	    
		pos_type max_ar=al;
		pos_type max_br=bl;
	    
		// get the maximal right ends of any arc match with left ends (al,bl)
		// in noLP mode, we don't consider cases without immediately enclosing arc match
		arc_matches_.get_max_right_ends(al,bl,
						&max_ar,&max_br,params_->no_lonely_pairs_);
	    
		// check whether there is an arc match at all
		if (al==max_ar || bl == max_br) continue;
	    
		// ------------------------------------------------------------
		// align under the maximal pair of arcs
		//
		align_in_arcmatch(al,max_ar,bl,max_br,params_->struct_local_);
	    
		//std::cout << al << ","<<bl<<":"<<std::endl
		//	      << Ms_[E_NO_NO] << std::endl;

		// ------------------------------------------------------------
		// fill D matrix entries
		//
		if (params_->no_lonely_pairs_) {
		    fill_D_entries_noLP(al,bl);
		} else {
		    fill_D_entries(al,bl);
		}
	    }
	}
    
	D_created_=true; // now the matrix D is built up
    }


    // align the top level in case of free end gaps
    //
    infty_score_t
    AlignerImpl::align_top_level_free_endgaps() {
    
	M_matrix_t &M=Ms_[E_NO_NO];
	infty_score_t max_score;
    
	init_state(E_NO_NO,
		   r_.get_startA()-1,r_.get_endA()+1,
		   r_.get_startB()-1,r_.get_endB()+1,
		   !free_endgaps_.allow_left_2(),false,
		   !free_endgaps_.allow_left_1(),false,
		   def_scoring_view_);
    
	// need to handle anchor constraints:
	// search maximum to the right of (or at) rightmost anchor constraint
	//
	AnchorConstraints::size_pair_t right_anchor = params_->constraints_->rightmost_anchor();
	//AnchorConstraints::size_pair_t left_anchor  = params_->constraints_->leftmost_anchor();
    
	for (pos_type i=r_.get_startA(); i<=r_.get_endA(); i++) {
	    Fs_[E_NO_NO]=infty_score_t::neg_infty;
	
	    // limit entries due to trace controller
	    pos_type min_col = std::max(r_.get_startB(),params_->trace_controller_->min_col(i));
	    pos_type max_col = std::min(r_.get_endB(),params_->trace_controller_->max_col(i));

	    for (pos_type j=min_col; j<=max_col; j++) {
		M(i,j) = align_noex( E_NO_NO, r_.get_startA()-1, r_.get_startB()-1, i, j,def_scoring_view_ );	      
	    }
	}
    
	//std::cout << "M-matrix:" <<std::endl << M << std::endl;

	max_score=M(r_.get_endA(),r_.get_endB());
	max_i_=r_.get_endA();
	max_j_=r_.get_endB();
    
    
	if (free_endgaps_.allow_right_2()) {
	    // search maximum in the rightmost row r_.get_endB()
	    // pay attention for anchor constraints AND trace controller
	
	    for (pos_type i=std::max(right_anchor.first+1,r_.get_startA()); i<=r_.get_endA(); i++) {
		if ( params_->trace_controller_->max_col(i)>=r_.get_endB() && M(i,r_.get_endB()) > max_score ) {
		    max_score = M(i,r_.get_endB());
		    max_i_=i; 
		    max_j_=r_.get_endB();
		}
	    }
	}
    
	if (free_endgaps_.allow_right_1()) {
	    // search maximum in the last column r_.get_endA()
	    // pay attention for anchor constraints AND trace controller
	
	    // limit entries due to trace controller
	    pos_type min_col = std::max(std::max(right_anchor.second+1,r_.get_startB()),params_->trace_controller_->min_col(r_.get_endA()));
	    pos_type max_col = std::min(r_.get_endB(),params_->trace_controller_->max_col(r_.get_endA()));

	
	    for (pos_type j=min_col; j<=max_col; j++) {
		if ( M(r_.get_endA(),j) > max_score ) {
		    max_score = M(r_.get_endA(),j);
		    max_i_=r_.get_endA();
		    max_j_=j;
		}
	    }
	}

    
    
	return max_score;
    }


    // align the top level in case of sequence local alignment
    //
    template<class ScoringView>
    infty_score_t
    AlignerImpl::align_top_level_locally(ScoringView sv) {
	//std::cout << r << std::endl;
    
	M_matrix_t &M=Ms_[E_NO_NO];
	infty_score_t max_score=infty_score_t::neg_infty;
	init_state(E_NO_NO,r_.get_startA()-1,r_.get_endA()+1,r_.get_startB()-1,r_.get_endB()+1,false,false,false,false,sv);
    
	// need to handle anchor constraints:
	// search maximum to the right of (or at) rightmost anchor constraint
	//
	AnchorConstraints::size_pair_t right_anchor = params_->constraints_->rightmost_anchor();
    
	AnchorConstraints::size_pair_t left_anchor = params_->constraints_->leftmost_anchor();

	//AnchorConstraints::size_pair_t right_anchor = AnchorConstraints::size_pair_t(r_.get_startA(),r_.get_startB());//dummy
    
	//std::cout << "right_anchor: "<<(right_anchor.first)<<","<<(right_anchor.second)<<std::endl;
    
	for (pos_type i=r_.get_startA(); i<=r_.get_endA(); i++) {
	    Fs_[E_NO_NO]=infty_score_t::neg_infty;

	    // limit entries due to trace controller
	    pos_type min_col = std::max(r_.get_startB(),params_->trace_controller_->min_col(i));
	    pos_type max_col = std::min(r_.get_endB(),params_->trace_controller_->max_col(i));

	    for (pos_type j=min_col; j<=max_col; j++) {
	    
		M(i,j) = align_noex( E_NO_NO,r_.get_startA()-1,r_.get_startB()-1,i,j, sv);	      
		//
		// score can be 0 (= drop prefix alignment) only if this is allowed due to constraints
		if ( i<left_anchor.first && j<left_anchor.second ) {
		    M(i,j) = std::max( (infty_score_t)0, M(i,j) );
		}
	    
		if (i>=right_anchor.first && j>=right_anchor.second && max_score < M(i,j)) {
		    max_score=M(i,j);
		    max_i_ = i;
		    max_j_ = j;
		}
	    }
	}
    
	// std::cout << "max: "<<max_i<<","<<max_j<<std::endl;
    
	// std::cout << M << std::endl;
    
	return max_score;
    }


    // special top level alignment for the scanning version
    // ATTENTION: no special anchor constraint handling done here (seems not very useful anyway)
    infty_score_t
    AlignerImpl::align_top_level_localB() {
	// std::cout <<"align local B " << r_.get_startA() << " " << r_.get_startB() << " "
	//           << r_.get_endA() << " " << r_.get_endB() << std::endl;
    
	M_matrix_t &M=Ms_[E_NO_NO];
	infty_score_t max_score=infty_score_t::neg_infty;
	init_state(E_NO_NO,r_.get_startA()-1,r_.get_endA()+1,r_.get_startB()-1,r_.get_endB()+1,true,false,false,false,def_scoring_view_);
    
	for (pos_type i=r_.get_startA(); i<=r_.get_endA(); i++) {

	    // limit entries due to trace controller
	    pos_type min_col = std::max(r_.get_startB(),params_->trace_controller_->min_col(i));
	    pos_type max_col = std::min(r_.get_endB(),params_->trace_controller_->max_col(i));

	    for (pos_type j=min_col; j<=max_col; j++) {
		M(i,j) = 
		    std::max( (infty_score_t)0,
			      align_noex( E_NO_NO, 
					  r_.get_startA()-1,r_.get_startB()-1,i,j,def_scoring_view_ ) );
		if (i==r_.get_endA() && max_score < M(i,j)) {
		    max_score=M(i,j);
		    max_i_ = i;
		    max_j_ = j;
		}
	    }
	}
	return max_score;
    }


    // compute the alignment score
    infty_score_t
    AlignerImpl::align() {
	// ------------------------------------------------------------
	// computes D matrix (if not already done) and then does the alignment on the top level
	// ------------------------------------------------------------
    
	if (!D_created_) align_D();

	if (params_->sequ_local_) {
	    return align_top_level_locally(def_scoring_view_);
	} else { // sequence global alignment
	
	    // align toplevel globally with potentially free endgaps (as given by description params_->free_endgaps)
	    return align_top_level_free_endgaps();
	
	    // Without free end gaps, we could call (much simpler):
	    // (we just call the more general method align_top_level_free_endgaps() since there is no performance penalty)
	    //
	    // align the subsequences that are specified by the restriction object
	    /*
	      align_in_arcmatch(r_.get_startA()-1,r_.get_endA()+1,r_.get_startB()-1,r_.get_endB()+1,false);
	      max_i=r_.get_endA();
	      max_j=r_.get_endB();
		
	      return Ms_[E_NO_NO](max_i,max_j);
	    */
	}
    }

    infty_score_t
    Aligner::align() {
	return pimpl_->align();
    }


    // ------------------------------------------------------------
    // Aligner: traceback

    void AlignerImpl::trace_arcmatch(const ArcMatch &am) {

	// std::cout << "trace_arcmatch " << am.arcA() << " " << am.arcB() <<std::endl;

	assert(params_->trace_controller_->is_valid_match(am.arcA().left(),am.arcB().left()));
	assert(params_->trace_controller_->is_valid_match(am.arcA().right(),am.arcB().right()));
    
	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();
    
	pos_type al=arcA.left();
	pos_type ar=arcA.right();
	pos_type bl=arcB.left();
	pos_type br=arcB.right();
    
	// --------------------
	// handle case of stacking
	if ( scoring_->stacking() ) {
	
	    if (arc_matches_.exists_inner_arc_match(am)) { 
		const ArcMatch &inner_am = arc_matches_.inner_arc_match(am);
	    
		if (D(am) == D(inner_am) + scoring_->arcmatch(am,true)) {
		    
		    const Arc & arcAI = inner_am.arcA();
		    const Arc & arcBI = inner_am.arcB();
		
		    alignment_.add_basepairA(arcAI.left(),arcAI.right());
		    alignment_.add_basepairB(arcBI.left(),arcBI.right());
		    alignment_.append(arcAI.left(),arcBI.left());
		
		    trace_arcmatch(inner_am);
		
		    alignment_.append(arcAI.right(),arcBI.right());
		
		    return;
		}
	    }
	}
    
	// --------------------
	// now handle the case that arc match is not stacked
    
	// first recompute M
	align_in_arcmatch(al,ar, bl,br,
			  params_->struct_local_);
    
	// then, trace in new M
	if (!params_->struct_local_) {
	    trace_in_arcmatch(0,al,ar-1,bl,br-1,false,def_scoring_view_);
	} else {
	    for (pos_type k=0; k<4; k++) {
		if (D(am) == Ms_[k](ar-1,br-1) + scoring_->arcmatch(am)) {
		    trace_in_arcmatch(k,al,ar-1,bl,br-1,false,def_scoring_view_);
		    break;
		}
	    }
	}
	return;
    }

    void AlignerImpl::trace_arcmatch_noLP(const ArcMatch &am) {
    
	//std::cout << "trace_arcmatch_noLP " << am.arcA() << " " << am.arcB() <<std::endl;
    
	assert(params_->trace_controller_->is_valid_match(am.arcA().left(),am.arcB().left()));
	assert(params_->trace_controller_->is_valid_match(am.arcA().right(),am.arcB().right()));
    
    
	assert(arc_matches_.exists_inner_arc_match(am));
    
	const ArcMatch &inner_am = arc_matches_.inner_arc_match(am);

	const Arc & arcAI = inner_am.arcA();
	const Arc & arcBI = inner_am.arcB();
    
	alignment_.add_basepairA(arcAI.left(),arcAI.right());
	alignment_.add_basepairB(arcBI.left(),arcBI.right());
	alignment_.append(arcAI.left(),arcBI.left());
    
	if (D(am) == D(inner_am) + scoring_->arcmatch(am,scoring_->stacking())) {
	    trace_arcmatch_noLP(inner_am);
	} else {
	    // first recompute M
	    align_in_arcmatch(arcAI.left(),arcAI.right(), arcBI.left(),arcBI.right(),
			      params_->struct_local_);
	
	    // then, trace in new M
	    if (!params_->struct_local_) {
		trace_in_arcmatch(0,arcAI.left(),arcAI.right()-1,arcBI.left(),arcBI.right()-1,false,def_scoring_view_);
	    } else {
		for (pos_type k=0; k<4; k++) {
		    if (D(am) == 
			Ms_[k](arcAI.right()-1,arcBI.right()-1)
			+ scoring_->arcmatch(am,scoring_->stacking()) + scoring_->arcmatch(inner_am))
			{
			    trace_in_arcmatch(k,arcAI.left(),arcAI.right()-1,arcBI.left(),arcBI.right()-1,false,def_scoring_view_);
			    break;
			}
		}
	    }
	}
	alignment_.append(arcAI.right(),arcBI.right());
    }

    // trace and handle all cases that do not involve exclusions
    template<class ScoringView>
    void
    AlignerImpl::trace_noex(int state,pos_type oal,pos_type i,
			    pos_type obl,pos_type j,
			    bool tl,
			    ScoringView sv) {
	M_matrix_t &M=Ms_[state];
    
	// determine where we get M(i,j) from
    
	// std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;
    
	// match
	if ( params_->constraints_->allowed_edge(i,j)
	     && params_->trace_controller_->is_valid(i-1,j-1)
	     && M(i,j) ==  M(i-1,j-1)+sv.scoring()->basematch(i,j) ) {
	    trace_in_arcmatch(state,oal,i-1,obl,j-1,tl,sv);
	    alignment_.append(i,j);
	    return;
	}

	if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
	    // del
	    if ( (!params_->constraints_->aligned_in_a(i))
		 && params_->trace_controller_->is_valid(i-1,j)
		 && M(i,j) == M(i-1,j)+sv.scoring()->gapA(i)) {
		trace_in_arcmatch(state,oal,i-1,obl,j,tl,sv);
		alignment_.append(i,-1);
		return;
	    }
	    // ins
	    if ( (!params_->constraints_->aligned_in_b(j))
		 && params_->trace_controller_->is_valid(i,j-1)
		 && M(i,j) == M(i,j-1)+sv.scoring()->gapB(j)) {
		trace_in_arcmatch(state,oal,i,obl,j-1,tl,sv);
		alignment_.append(-1,j);
		return;
	    }
	} else { // base del and ins, affine cost
	    // since, we didn't store all values in the E and F matrix
	    // we do the traceback in linear time per entry
	    // base del
	    score_t gap_cost=sv.scoring()->indel_opening();
	    for (pos_type k=1;
		 (i >= oal+k)
		     && (! params_->constraints_->aligned_in_a(i-k+1));
		 k++)
		{
		    // break if gap becomes invalid due to trace controller
		    // (it is safe to break, because of trace controller's monotonicity)
		    if (! params_->trace_controller_->is_valid(i-k,j)) break;

		    gap_cost += sv.scoring()->gapA(i-k+1);
	    
		    if ( M(i,j) == M(i-k,j) + gap_cost) {
			// gap in A of length k
			trace_in_arcmatch(state,oal,i-k,obl,j,tl,sv);
			for (pos_type l=k;l>0;l--) {
			    alignment_.append(i-l+1,-1);
			}
			return;
		    }
		}
	
	    // base ins
	    gap_cost=sv.scoring()->indel_opening();
	    for (pos_type k=1; 
		 (j >= obl+k)
		     && (! params_->constraints_->aligned_in_b(j-k+1));
		 k++)
		{
		    // break if gap becomes invalid due to trace controller
		    // (it is safe to break, because of trace controller's monotonicity)
		    if (! params_->trace_controller_->is_valid(i,j-k)) break;

		    gap_cost += sv.scoring()->gapB(j-k+1);
		    if (M(i,j) == M(i,j-k) + gap_cost) {
			// gap in B of length k
			trace_in_arcmatch(state,oal,i,obl,j-k,tl,sv);
			for (pos_type l=k;l>0;l--) {
			    alignment_.append(-1,j-l+1);
			}
		
			return;
		    }
		}
	}

	// only consider arc match cases if edge (i,j) is allowed and valid!
	if ( ! (params_->constraints_->allowed_edge(i,j) 
		&& params_->trace_controller_->is_valid(i-1,j-1))
	     ) return;
    
    
	// here (i,j) is allowed and valid
    
	//  arc match
	const pos_type &ar=i;
	const pos_type &br=j;
    
	for(ArcMatchIdxVec::const_iterator it=arc_matches_.common_right_end_list(ar,br).begin();
	    arc_matches_.common_right_end_list(ar,br).end() != it; ++it ) {
	
	    // NOTES: *it is the arc match index
	    //        we iterate only over valid arc matches, i.e.
	    //        constraints (including anchor c. and heuristic ones) are satisified
	
	    const ArcMatch &am = arc_matches_.arcmatch(*it);
	
	    const Arc &arcA=am.arcA();
	    const Arc &arcB=am.arcB();
	
	    if ( (arcA.left() <= oal) || (arcB.left() <= obl) ) continue;

	    const pos_type al=arcA.left();
	    const pos_type bl=arcB.left();
	    
	    if ( M(i,j) == M(al-1,bl-1) + sv.D(am)) {
		//
		// do the trace for alignment left of the arc match
		trace_in_arcmatch(state,oal,al-1,obl,bl-1,tl,sv);
	    
		//cout << "arcmatch "<<(al)<<","<<i<<";"<<(bl)<<","<<j<<" :: "
		//      <<(arcA->w)<<" + "<<(arcB->w)<< " + " << tau(al,bl,i,j)  <<endl;
	    
		alignment_.add_basepairA(al,ar);
		alignment_.add_basepairB(bl,br);
		alignment_.append(al,bl);
	    
		// do the trace below the arc match
	    
		if (params_->no_lonely_pairs_) {
		    trace_arcmatch_noLP(am);
		} else {
		    trace_arcmatch(am);
		}
	    
		alignment_.append(ar,br);
	    
		return;
	    }
	}
    }

    // do the trace within one arc match.
    // the cases without exclusions are delegated to trace_noex
    template <class ScoringView>
    void
    AlignerImpl::trace_in_arcmatch(int state,int al,int i,int bl,int j,bool tl,ScoringView sv) {
	//pre: M matrices for arc computed
	M_matrix_t &M=Ms_[state];

    
	// string state_text[]={"E_NO_NO", "E_X_NO", "E_NO_X", "E_X_X",
	// 			 "E_OP_NO", "E_NO_OP", "E_OP_X", "E_X_OP"};
	// cout << "trace_in_arcmatch "<<state_text[state]<<" al:"<<al<<" i:"<<i
	// 	 <<" bl:"<<bl<<" j:"<<j<<" :: "<< M(i,j) <<endl;
    
    
	assert(params_->trace_controller_->is_valid(i,j));
    
    
	// terminate traceback if
	// * trace on toplevel
	// * sequence local, and
	// * entry == 0
	if (tl && params_->sequ_local_ && M(i,j)==(infty_score_t)0) {
	    min_i_=i;
	    min_j_=j;
	    return;
	}

	if (i<=al) {
	    if (state==E_NO_NO
		|| state==E_OP_NO
		|| state==E_X_NO) {
	    	    
		// pad with gap edges, unless in special cases (local/semi-local alignment)
		if (!(tl && 
		      (params_->sequ_local_ || free_endgaps_.allow_left_1()))) {
		    for (int k=bl+1;k<=j;k++) {
			alignment_.append(-1,k);
		    }
		}
	    } else {
		//cout << "exclude B "<<state<<" "<< bl<<" - "<<j<<endl;
	    }
	    return;
	}

	if (j<=bl) {
	    if (state==E_NO_NO
		|| state==E_NO_OP
		|| state==E_NO_X) {
	    
		// pad with gap edges, unless in special cases (local/semi-local alignment)
		if (!(tl && 
		      ( params_->sequ_local_ || free_endgaps_.allow_left_2()))) {
		    for (int k=al+1;k<=i;k++) {
			alignment_.append(k,-1);
		    }
		}
	    } else {
		//cout << "exclude A "<<state<<" "<< al<<" - "<<i<<endl;
	    }
	    return;
	}


	switch(state) {
	case E_NO_NO:
	    trace_noex(state,al,i,bl,j,tl,sv);
	    break;
	case E_OP_NO:
	    if (M(i,j) == M(i-1,j)) {
		// cout << "exclude A "<<i<<endl;
		trace_in_arcmatch(state,al,i-1,bl,j,tl,sv);
	    } else if (M(i,j) == Ms_[E_NO_NO](i,j))
		trace_in_arcmatch(E_NO_NO,al,i,bl,j,tl,sv);
	    break;
	case E_NO_OP:
	    if (M(i,j) == M(i,j-1)) {
		// cout << "exclude B "<<j<<endl;
		trace_in_arcmatch(state,al,i,bl,j-1,tl,sv);
	    } else if (M(i,j) == Ms_[E_NO_NO](i,j))
		trace_in_arcmatch(E_NO_NO,al,i,bl,j,tl,sv);
	    break;
	case E_NO_X:
	    if (M(i,j) == Ms_[E_NO_OP](i,j)+sv.scoring()->exclusion()) {
		trace_in_arcmatch(E_NO_OP,al,i,bl,j,tl,sv);
	    } else
		trace_noex(state,al,i,bl,j,tl,sv);
	    break;
	case E_OP_X:
	    if (M(i,j) == M(i-1,j)) {
		//cout << "exclude A "<<i<<endl;
		trace_in_arcmatch(state,al,i-1,bl,j,tl,sv);
	    } else if (M(i,j) == Ms_[E_NO_X](i,j))
		trace_in_arcmatch(E_NO_X,al,i,bl,j,tl,sv);
	    break;
	case E_X_NO:
	    if (M(i,j) == Ms_[E_OP_NO](i,j)+sv.scoring()->exclusion()) {
		trace_in_arcmatch(E_OP_NO,al,i,bl,j,tl,sv);
	    } else
		trace_noex(state,al,i,bl,j,tl,sv);
	    break;
	case E_X_OP:
	    if (M(i,j) == M(i,j-1)) {
		// cout << "exclude B "<<j<<endl;
		trace_in_arcmatch(state,al,i,bl,j-1,tl,sv);
	    } else if (M(i,j) == Ms_[E_X_NO](i,j))
		trace_in_arcmatch(E_X_NO,al,i,bl,j,tl,sv);
	    break;
	case E_X_X:
	    if (M(i,j) == Ms_[E_OP_X](i,j)+sv.scoring()->exclusion()) {
		trace_in_arcmatch(E_OP_X,al,i,bl,j,tl,sv); 
	    } else if (M(i,j) == Ms_[E_X_OP](i,j)+sv.scoring()->exclusion()) {
		trace_in_arcmatch(E_X_OP,al,i,bl,j,tl,sv);
	    } else
		trace_noex(state,al,i,bl,j,tl,sv);
	    break;
	}
    }


    template<class ScoringView>
    void
    AlignerImpl::trace(ScoringView sv) {
	// pre: last call align_in_arcmatch(r_.get_startA()-1,r_.get_endA()+1,r_.get_startB()-1,r_.get_endB()+1);
	//      or align_top_level_locally for sequ_local_ alignent
    
	// reset the alignment strings (to empty strings)
	// such that they can be written again during the trace
	alignment_.clear();
    
	trace_in_arcmatch(E_NO_NO,r_.get_startA()-1,max_i_,r_.get_startB()-1,max_j_,true,sv);
    }

    void
    Aligner::trace() {pimpl_->trace(pimpl_->def_scoring_view_);}



    /* ------------------------------------------------------------
       Compute k-best alignments by intervall splitting
    */

    void 
    Aligner::set_restriction(const AlignerRestriction &r) {pimpl_->r_=r;}
    
    const AlignerRestriction &
    Aligner::get_restriction() const {return pimpl_->r_;}



    //! type of a task (used in computing k-best alignment)
    typedef std::pair<AlignerRestriction,infty_score_t> task_t;

    
    void Aligner::suboptimal(size_t k,
				 score_t threshold,
				 bool opt_normalized,
				 score_t normalized_L,
				 size_t output_width,
				 bool opt_verbose,
				 bool opt_local_output,
				 bool opt_pos_output,
				 bool opt_write_structure
				 ) {
	Aligner &a=*this;
	
	// compute alignment score for a
	infty_score_t a_score = (!opt_normalized)?a.align():a.normalized_align(normalized_L,false);
    
	// make priority queue tasks
	std::priority_queue<task_t, std::vector<task_t>, greater_second<task_t> > tasks;
	
	// put a into tasks
	tasks.push(task_t(a.get_restriction(),a_score));
	
	size_t i=1;
    
	while ( ( k<0 || i<=k ) ) {
	    //std::cout << "queue size: "<<tasks.size()<<std::endl;
	
	    // get best task from tasks
	    task_t task = tasks.top();
	
	    // pop topmost element
	    tasks.pop();
	
	    AlignerRestriction &task_r=task.first;
	    infty_score_t task_score=task.second;
	
	    if ( task_score < (infty_score_t)threshold+1 ) break;
	
	    a.set_restriction(task_r);
	
	    if (!opt_normalized) {
		a.align(); // alignment needs to be recomputed (this is fast)!
		// do backtrace and write tasks alignment
		a.trace();
	    } else {
		a.normalized_align(normalized_L,opt_verbose);
	    }
	
	    Alignment alignment = a.get_alignment();
	
	    // std::cout << "SCORE: " << task_score << std::endl;
	    
	    {
		// after major code changes, the following output was
		// stripped off the formerly present features like to
		// write structure and ?. TODO: check and
		// reimplement. SW - 2013 Jun 7
		if (opt_pos_output) {
		    std::cout << "HIT "<<task_score
			      <<alignment.local_startA()<<" "
			      <<alignment.local_startB()<<" "
			      <<alignment.local_endA()<<" "
			      <<alignment.local_endB()<<" "
			      <<std::endl;
		} else {
		    MultipleAlignment ma(alignment,true);
		    std::cout << "Score: "<<task_score<<std::endl;
		    ma.write(std::cout,120);
		}
	    }
	    
	    if (!opt_pos_output) std::cout << std::endl
					   << std::endl;
	
	    if (i==k) break; // break if enough solutions generated
	
	    // split the longer sequence according to local alignment
	    pos_type lenA=task_r.get_endA()-task_r.get_startA();
	    pos_type lenB=task_r.get_endB()-task_r.get_startB();
	
	    // make two clones of AlignerRestriction task_r
	    AlignerRestriction r1(task_r);
	    AlignerRestriction r2(task_r);
	
	    if (lenA>lenB) {
		// split A
		int splitA = (alignment.local_startA() + alignment.local_endA())/2;
		if (opt_verbose) std::cout <<"Split A at "<<splitA<<std::endl;
		r1.set_endA(splitA);
		r2.set_startA(splitA);
	    } else {
		int splitB = (alignment.local_startB() + alignment.local_endB())/2;
		if (opt_verbose) std::cout <<"Split B at "<<splitB<<std::endl;
		r1.set_endB(splitB);
		r2.set_startB(splitB);
	    }
	
	    // compute alignment scores for both splits
	
	    a.set_restriction(r1);
	    infty_score_t a1_score = (!opt_normalized)?a.align():a.normalized_align(normalized_L,false);
	    a.set_restriction(r2);
	    infty_score_t a2_score = (!opt_normalized)?a.align():a.normalized_align(normalized_L,false);
	
	    // std::cout <<"a1_score: " << a1_score << std::endl;
	    // std::cout <<"a2_score: " << a2_score << std::endl;
		
	    // put both splits into <tasks>
	
	    tasks.push(task_t(r1,a1_score));
	    tasks.push(task_t(r2,a2_score));
		
	    ++i; // count enumerated alignments
	}
    }
    

    // ------------------------------------------------------------
    // Normalized local alignment using Dinkelbach's algorithm for
    // fractional programming
    //

    infty_score_t
    Aligner::normalized_align(score_t L, bool opt_verbose) {
    
	// The D matrix is filled as in non-normalized alignment. Because
	// alignments of the subsequences enclosed by arcs are essentially
	// global, their scores can be optimized in the same way as for
	// non-normalized alignments.
	if (!pimpl_->D_created_) pimpl_->align_D();

	if (pimpl_->mod_scoring_) delete pimpl_->mod_scoring_;
	pimpl_->mod_scoring_=new Scoring(*pimpl_->scoring_); // make mod_scoring point to a copy of scoring
    
	// Apply Dinkelbach's algorithm
    
	score_t new_lambda=0;
	score_t lambda=-1;

	size_t iteration=0;
    
	// iterate until convergence
	while ( lambda != new_lambda )
	    {
		++iteration;
		if (opt_verbose) std::cout << "Perform Dinkelbach iteration "<<iteration<<std::endl;
	
		lambda=new_lambda;
		
		// modify the scoring by lambda
		pimpl_->mod_scoring_->modify_by_parameter(lambda);
		pimpl_->mod_scoring_view_.set_lambda(lambda);
	
		infty_score_t score = pimpl_->align_top_level_locally(pimpl_->mod_scoring_view_);
	
		pimpl_->alignment_.clear();
	
		// perform a traceback for normalized alignment
		pimpl_->trace(pimpl_->mod_scoring_view_);
	
		// compute length (of alignment) as sum of lengths of
		// aligned subsequences from the trace
		pos_type length=pimpl_->max_i_-pimpl_->min_i_+1+pimpl_->max_j_-pimpl_->min_j_+1;
	
		// get score for the best alignment in the modified problem
		// but for unmodified scoring. Because for each position,
		// lambda was subtracted, we can just add length*lambda!
		score += length*lambda;
	
		new_lambda = score.finite_value()/(length+L);
	
		if (opt_verbose) std::cout << "Score: "<<score<<" Length: "<<length<<" Normalized Score: "<<new_lambda<<std::endl;

		if (opt_verbose) {
		    MultipleAlignment ma(pimpl_->alignment_,true);
		    std::cout << "Score: "<<(infty_score_t)new_lambda<<std::endl;
		    ma.write(std::cout,120);
		}
		
		if (opt_verbose) std::cout<<std::endl;
	    }
	return (infty_score_t)new_lambda;
    }

} //end namespace LocARNA
