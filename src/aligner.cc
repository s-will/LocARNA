#include "aligner.hh"
#include "anchor_constraints.hh"
#include "trace_controller.hh"
// #include "d_matrix.hh"

#include <math.h>
#include <assert.h>

#include <queue>

#include <iostream>


using namespace std;


/*
  NEW SEMANTIC OF CONSTRAINTS
  
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

Aligner::Aligner(const Aligner &a)
    :scoring(a.scoring),
     mod_scoring(0),
     params(a.params),
     seqA(a.seqA),
     seqB(a.seqB),
     arc_matches(a.arc_matches),
     bpsA(a.bpsA),
     bpsB(a.bpsB),
     r(a.r),
     Dmat(a.Dmat),
     Ms(a.Ms),
     Es(a.Es),
     Fs(a.Fs),
     min_i(a.min_i),
     min_j(a.min_j),
     max_i(a.max_i),
     max_j(a.max_j),
     D_created(a.D_created),
     alignment(a.alignment),
     def_scoring_view(this),
     mod_scoring_view(this)
{}

Aligner::Aligner(const Sequence &seqA_, 
		 const Sequence &seqB_,
		 const ArcMatches &arc_matches_,
		 const AlignerParams *ap_,
		 const Scoring *s_
		 )
    : scoring(s_),
      mod_scoring(0),
      params(ap_),
      seqA(seqA_), seqB(seqB_),
      arc_matches(arc_matches_),
      bpsA(arc_matches_.get_base_pairsA()),
      bpsB(arc_matches_.get_base_pairsB()),
      r(1,1,seqA_.length(),seqB_.length()),
      D_created(false),
      alignment(seqA_,seqB_),
      def_scoring_view(this),
      mod_scoring_view(this)
{
    Ms.resize(params->STRUCT_LOCAL?8:1);
    Es.resize(params->STRUCT_LOCAL?4:1);
    Fs.resize(params->STRUCT_LOCAL?4:1);
    
    Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
    Dmat.fill(infty_score_t::neg_infty);
    
    for (size_type k=0; k<(params->STRUCT_LOCAL?8:1); k++) {
	Ms[k].resize(seqA.length()+1,seqB.length()+1);
    }
    for (size_type k=0; k<(params->STRUCT_LOCAL?4:1); k++) {
	Es[k].resize(seqB.length()+1);
    }
}



Aligner::~Aligner() {
    if (mod_scoring!=0) delete mod_scoring;
}

// standard cases in alignment: base match, base in/del, arc match
// (for structure local alignment this is extended by exclusion handling)
// 
// if lonely basepairs are disallowed, there is special treatment
//
// align_noex has a side effect: it computes entry in Es[state] and Fs[state]
template<class ScoringView>
infty_score_t
Aligner::align_noex(int state, size_type al, size_type bl, size_type i, size_type j,ScoringView sv) {
    
    assert(0<=state && state<4);
    
    assert(params->trace_controller.is_valid(i,j));
    
    M_matrix_t &M = Ms[state];
    ScoreVector &E = Es[state];
    infty_score_t &F = Fs[state];
    
    // compute E entry
    if ( (! params->constraints.aligned_in_a(i)) ) {
      // due to constraints, i can be deleted
      E[j] = 
	std::max( E[j] + sv.scoring()->gapA(i,j),
		  M(i-1,j) + sv.scoring()->gapA(i,j) + sv.scoring()->indel_opening() );
      E[j] = E[j].normalized_neg();
    } else {
      // due to constraints, i cannot be deleted
      E[j] = infty_score_t::neg_infty;
    }
    
    // compute F entry
    if ( (! params->constraints.aligned_in_b(j)) ) {
      // due to constraints, j can be inserted
      F=std::max( F + sv.scoring()->gapB(i,j),
		  M(i,j-1) + sv.scoring()->gapB(i,j) + sv.scoring()->indel_opening() );
      F = F.normalized_neg();
    } else {
      // due to constraints, j cannot be inserted
      F = infty_score_t::neg_infty;
    }
    

    infty_score_t max_score = infty_score_t::neg_infty;
    
    // base match
    if (params->constraints.allowed_edge(i,j)) {
	max_score = M(i-1,j-1) + sv.scoring()->basematch(i,j);
    }
    
    // base del
    max_score=max(max_score, E[j]);
    
    // base ins
    max_score=max(max_score, F);

    // arc match
    
    // standard case for arc match (without restriction to lonely pairs)
    //    
    
    if ( params->constraints.allowed_edge(i,j) ) {
	const BasePairs::RightAdjList &adjlA = bpsA.right_adjlist(i);
	const BasePairs::RightAdjList &adjlB = bpsB.right_adjlist(j);
	
	// for all pairs of arcs in A and B that have right ends i and j, respectively
	//
	for (BasePairs::RightAdjList::const_iterator arcA=adjlA.begin();
	     arcA!=adjlA.end() && arcA->left() > al  ; ++arcA) {
	    for (BasePairs::RightAdjList::const_iterator arcB=adjlB.begin();
		 arcB!=adjlB.end() && arcB->left() > bl ; ++arcB) {
		
		// no need to check (params->constraints.allowed_edge(arcA->left(),arcB->left()))
		// or other "constraints"
		// because for these arc matches holds that sv.D(*arcA,*arcB)==neg_infty
		
		
		infty_score_t new_score =
		  M(arcA->left()-1,arcB->left()-1)
		  + sv.D(*arcA,*arcB);
		
		if (new_score > max_score) {
		  //std::cout << *arcA << "-"<< *arcB << ": "<<M(arcA->left()-1,arcB->left()-1)<<"+"<<D(*arcA,*arcB)<<"="<<new_score<<std::endl;
		  max_score=new_score;
		}
	    }
	}
    }
    
    return max_score.normalized_neg();

    
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
	
// 	max_score = max( max_score, M(arcA.left()-1,arcB.left()-1) + D[am.idx()] );
	
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
Aligner::init_state(int state, size_type al, size_type ar, size_type bl, size_type br, 
		    bool globalA, bool exclA,
		    bool globalB, bool exclB, 
		    ScoringView sv) {
    
    // alignments that have empty subsequence in A (i=al) and
    // end with gap in alistr of B do not exist ==> -infty 
    if (state<4) {
	ScoreVector   &E = Es[state];
	for (size_type j=bl; j<br; j++) {
	    E[j]=infty_score_t::neg_infty;
	}
    }
    
    M_matrix_t &M = Ms[state];
    
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
    
    size_type i;
    for (i=al+1; i<ar; i++) {

	if (params->trace_controller.min_col(i)>bl) break; // fill only as long as column bl is accessible

	if (!indel_score.is_neg_infty()) {
	    if (params->constraints.aligned_in_a(i)) {
		indel_score=infty_score_t::neg_infty;
	    }
	    else if (!exclA && globalA) {
		indel_score += sv.scoring()->gapA(i,bl);
	    }
	}
	M(i,bl) = (infty_score_t)indel_score;
    }

    // fill entries left of valid entries 
    for ( ; i<ar; i++) {
	assert(params->trace_controller.min_col(i)>bl);
	M(i,params->trace_controller.min_col(i)-1) = infty_score_t::neg_infty; 
    }
    
    // init first row al
    //
    indel_score=(infty_score_t)sv.scoring()->indel_opening();
    if (exclB) {
      indel_score = (infty_score_t)sv.scoring()->exclusion();
    } else if (!globalB) {
      indel_score = (infty_score_t)0;
    }
    
    size_type j;
    for (j=bl+1 ; j < min(br, params->trace_controller.max_col(al)+1) ; j++) {
      if (!indel_score.is_neg_infty()) {
	if (params->constraints.aligned_in_b(j)) {
	  indel_score=infty_score_t::neg_infty;
	}
	else if (!exclB && globalB && !indel_score.is_neg_infty()) {
	  indel_score += sv.scoring()->gapB(al,j);
	}
      }
      M(al,j) = (infty_score_t)indel_score;
    }

    // fill entries above valid entries 
    // here j points to one position right of the last initialized entry in row al
    for (i=al+1; i<ar; i++) {
	for (;
	     j<min(br,params->trace_controller.max_col(i)+1); ++j) {
	    M(i-1,j)=infty_score_t::neg_infty;
	}
    }

}

// ----------------------------------------
// recomputes M matrix/matrices
// after the call the matrix is filled in the range [al..ar-1] x [bl..br-1]
void
Aligner::align_in_arcmatch(size_type al,size_type ar,size_type bl,size_type br,
				bool allow_exclusion) {

    assert(br>0); // if br<=0 we run into trouble below when computing br-1

    // cout << al << " " << ar <<" " << bl << " " << br <<endl;

    // When using RMatrix as type of the M matrices,
    // the M matrix/matrices can be restricted to the range of the arc match
    // for (size_type state=0; state < ((allow_exclusion)?8:1); state++) {
	//std::cout <<state<<" "<<Ms[state].sizes().first<<" "<<Ms[state].sizes().second<<" "<<al<<" "<<ar<<" "<<bl<<" "<<br<<std::endl; 
	//Ms[state].restrict(al,ar-1,bl,br-1);
    // }
    
    
    // if in a sequence the state is not open than gaps with cost scoring->gap() have to be introduced.
    // In open cases the maximal similarity is 0 in the first row/column,
    // since an exclusion can be introduced.
    
    init_state(E_NO_NO,al,ar,bl,br,true ,false,true ,false,def_scoring_view);
    
    if (allow_exclusion) {
	init_state(E_X_NO, al,ar,bl,br, true , true , true , false,def_scoring_view);
	init_state(E_NO_X, al,ar,bl,br, true , false, true , true ,def_scoring_view);
	init_state(E_X_X,  al,ar,bl,br, true , true , true , true ,def_scoring_view);
	
	// open states
	init_state(E_OP_NO,al,ar,bl,br, false, false, true , false,def_scoring_view);
	init_state(E_NO_OP,al,ar,bl,br, true , false, false, false,def_scoring_view);
	init_state(E_X_OP, al,ar,bl,br, true , true , false, false,def_scoring_view);
	init_state(E_OP_X, al,ar,bl,br, false, false, true , true ,def_scoring_view);
    }

    // ----------------------------------------
    // alignment for state E_NO_NO
    //
    
    for (size_type i=al+1; i<ar; i++) {
	Fs[E_NO_NO]=infty_score_t::neg_infty;
	
	// limit entries due to trace controller
	size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	
	for (size_type j=min_col; j<=max_col; j++) {
	    Ms[E_NO_NO](i,j)=align_noex(E_NO_NO,al,bl,i,j,def_scoring_view);
	}
    }
    
    //
    // end state E_NO_NO
    // ----------------------------------------
    
    if (allow_exclusion) {
	int state;

	state=E_OP_NO;
	for (size_type i=al+1; i<ar; i++) { 
	
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(
				     params->constraints.aligned_in_a(i)?infty_score_t::neg_infty:Ms[state](i-1,j),
				     Ms[E_NO_NO](i,j)
				     );
	    }
	}
	    
	state=E_NO_OP;
	for (size_type i=al+1; i<ar; i++) {
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(
				     params->constraints.aligned_in_b(j)?infty_score_t::neg_infty:Ms[state](i,j-1),
				     Ms[E_NO_NO](i,j)
				     );
	    }
	}

	state=E_NO_X;
	for (size_type i=al+1; i<ar; i++) {
	    Fs[state]=infty_score_t::neg_infty;
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(align_noex(state,al,bl,i,j,def_scoring_view),
				     Ms[E_NO_OP](i,j)+scoring->exclusion());
	    }
	}
	
	state=E_OP_X;
	for (size_type i=al+1; i<ar; i++) {
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(
				     params->constraints.aligned_in_a(i)?infty_score_t::neg_infty:Ms[state](i-1,j),
				     Ms[E_NO_X](i,j)
				     );
	    }
	}
	
	state=E_X_NO;
	for (size_type i=al+1; i<ar; i++) {
	    Fs[state]=infty_score_t::neg_infty;
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(align_noex(state,al,bl,i,j,def_scoring_view),
				     Ms[E_OP_NO](i,j)+scoring->exclusion());
	    }
	}

	state=E_X_OP;
	for (size_type i=al+1; i<ar; i++) {
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = max(
				     params->constraints.aligned_in_b(j)?infty_score_t::neg_infty:Ms[state](i,j-1),
				     Ms[E_X_NO](i,j)
				     );
	    }
	}

	state=E_X_X;
	for (size_type i=al+1; i<ar; i++) {
	    Fs[state]=infty_score_t::neg_infty;
	    // limit entries due to trace controller
	    size_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    size_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    
	    for (size_type j=min_col; j<=max_col; j++) {
		Ms[state](i,j) = 
		    max(align_noex(state,al,bl,i,j,def_scoring_view),
			max(Ms[E_OP_X](i,j)+scoring->exclusion(),
			    Ms[E_X_OP](i,j)+scoring->exclusion())); }
	}
    }
}

// compute the entries in the D matrix that
// can be computed from the matrix/matrices M
// for the subproblem al,bl,max_ar,max_br
//
// pre: M matrices are computed by a call to 
//      align_in_arcmatch(al,max_ar,bl,max_br,params->STRUCT_LOCAL)
//
void 
Aligner::fill_D_entries(size_type al, size_type bl)
{
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();
	arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();

	
	size_type ar = arcA.right();
	size_type br = arcB.right();
	
	infty_score_t m=Ms[0](ar-1,br-1);
	if (params->STRUCT_LOCAL) {
	    // if we align structure local
	    // we need to determine the maximum of entries in Ms[k],
	    // where k is a non-open state
	    for (size_type k=1; k<4; k++) {
		m=max(m, Ms[k](ar-1,br-1));
	    }
	}
	
	D(am) = m + scoring->arcmatch(am);
	
	//std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << std::endl;

	if (scoring->stacking()) {
	    if (arc_matches.is_stackable(am)) {
		const ArcMatch &inner_am = arc_matches.inner_arc_match(am);
		
		D(am) =
		    max(D(am),
			D(inner_am) + scoring->arcmatch_stacked(am));
	    }
	}
    }
}

void 
Aligner::fill_D_entries_noLP(size_type al, size_type bl) {
    // get adj lists of arcs starting in al-1, bl-1
    
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al-1,bl-1).begin();
	arc_matches.common_left_end_list(al-1,bl-1).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	size_type ar = am.arcA().right()-1;
	size_type br = am.arcB().right()-1;

	// only for arc matches which can occur paired with an inner arc match (no lonely pairs) 
	//
	// when using stacking scores, the am has to be stackable to the inner arc,
	// without using stacking scores, it suffices that an inner arc exists!
	if ((scoring->stacking() && arc_matches.is_stackable(am))
	    || (!scoring->stacking() && arc_matches.exists_inner_arc_match(am))) { 
	    const ArcMatch& inner_am = arc_matches.inner_arc_match(am);
	    
	    infty_score_t m=Ms[0](ar-1,br-1);
	    if (params->STRUCT_LOCAL) {
		// if we align structure local
		// we need to determine the maximum of entries in Ms[k],
		// where k is a non-open state
		for (size_type k=1; k<4; k++) {
		    m=max(m, Ms[k](ar-1,br-1));
		}
	    }
	    
	    //std::cout << scoring->arcmatch(am,scoring->stacking()) <<" "<< scoring->arcmatch(inner_am) <<" "<< m << " ";
	    
	    // score of outer arc match, where inner arc match is aligned too
	    D(am) = 
		max(m + scoring->arcmatch(inner_am),
		    D(inner_am))
		+ scoring->arcmatch(am,scoring->stacking());

	    //std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << " ; ";
	    //std::cout <<"["<< inner_am.arcA() << "," <<inner_am.arcB() <<"]:" <<D(inner_am)<< std::endl;

	}
    }
}



// compute all entries D
void
Aligner::align_D() {
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
    // for al in r.get_endA() .. r.get_startA
    for (size_type al=r.get_endA()+1; al>r.get_startA(); ) { al--; 
	
	size_type max_bl = min(r.get_endB(),params->trace_controller.max_col(al));
	size_type min_bl = max(r.get_startB(),params->trace_controller.min_col(al));
	
	// for bl in max_bl .. min_bl
	for (size_type bl=max_bl+1; bl > min_bl;) { bl--; 
	    
	    if (! ( params->constraints.allowed_edge(al,bl)
		    && params->trace_controller.is_valid_match(al,bl) )
		) continue;
	    
	    // ------------------------------------------------------------
	    // get maximal right ends of arcs with left ends al,bl 
	    // where max_diff_am conditions hold
	    // and no_lonely_pairs condition holds
	    //
	    
	    size_type max_ar=al;
	    size_type max_br=bl;
	    
	    // get the maximal right ends of any arc match with left ends (al,bl)
	    arc_matches.get_max_right_ends(al,bl,&max_ar,&max_br,params->no_lonely_pairs);
	    
	     // check whether there is an arc match at all
	    if (al==max_ar || bl == max_br) continue;
	    
	    // ------------------------------------------------------------
	    // align under the maximal pair of arcs
	    //
	    align_in_arcmatch(al,max_ar,bl,max_br,params->STRUCT_LOCAL);
	    
	    //std::cout << al << ","<<bl<<":"<<std::endl
	    //	      << Ms[E_NO_NO] << std::endl;

	    // ------------------------------------------------------------
	    // fill D matrix entries
	    //
	    if (params->no_lonely_pairs) {
		fill_D_entries_noLP(al,bl);
	    } else {
		fill_D_entries(al,bl);
	    }
	}
    }
    
    D_created=true; // now the matrix D is built up
}


// align the top level in case of free end gaps
//
infty_score_t
Aligner::align_top_level_free_endgaps() {
    
    M_matrix_t &M=Ms[E_NO_NO];
    infty_score_t max_score;
    
    init_state(E_NO_NO,
	       r.get_startA()-1,r.get_endA()+1,
	       r.get_startB()-1,r.get_endB()+1,
	       !params->free_endgaps.allow_left_2(),false,
	       !params->free_endgaps.allow_left_1(),false,
	       def_scoring_view);
    
    // need to handle anchor constraints:
    // search maximum to the right of (or at) rightmost anchor constraint
    //
    AnchorConstraints::size_pair_t right_anchor = params->constraints.rightmost_anchor();
    AnchorConstraints::size_pair_t left_anchor  = params->constraints.leftmost_anchor();
    
    for (size_type i=r.get_startA(); i<=r.get_endA(); i++) {
	Fs[E_NO_NO]=infty_score_t::neg_infty;
	
	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(i));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(i));

	for (size_type j=min_col; j<=max_col; j++) {
	    M(i,j) = align_noex( E_NO_NO, r.get_startA()-1, r.get_startB()-1, i, j,def_scoring_view );	      
	}
    }
    
    //std::cout << "M-matrix:" <<std::endl << M << std::endl;

    max_score=M(r.get_endA(),r.get_endB());
    max_i=r.get_endA();
    max_j=r.get_endB();
    
    
    if (params->free_endgaps.allow_right_2()) {
	// search maximum in the rightmost row r.get_endB()
	// pay attention for anchor constraints
	
	for (size_type i=std::max(right_anchor.first+1,r.get_startA()); i<=r.get_endA(); i++) {
	    if ( M(i,r.get_endB()) > max_score ) {
		max_score = M(i,r.get_endB());
		max_i=i; 
		max_j=r.get_endB();
	    }
	}
    }
    
    if (params->free_endgaps.allow_right_1()) {
	// search maximum in the last column r.get_endA()
	// pay attention for anchor constraints
	
	for (size_type j=std::max(right_anchor.second+1,r.get_startB()); j<=r.get_endB(); j++) {
	    if ( M(r.get_endA(),j) > max_score ) {
		max_score = M(r.get_endA(),j);
		max_i=r.get_endA();
		max_j=j;
	    }
	}
    }

    
    
    return max_score;
}


// align the top level in case of sequence local alignment
//
template<class ScoringView>
infty_score_t
Aligner::align_top_level_locally(ScoringView sv) {
    //std::cout << r << std::endl;
    
    M_matrix_t &M=Ms[E_NO_NO];
    infty_score_t max_score=infty_score_t::neg_infty;
    init_state(E_NO_NO,r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1,false,false,false,false,sv);
    
    // need to handle anchor constraints:
    // search maximum to the right of (or at) rightmost anchor constraint
    //
    AnchorConstraints::size_pair_t right_anchor = params->constraints.rightmost_anchor();
    
    AnchorConstraints::size_pair_t left_anchor = params->constraints.leftmost_anchor();

    //AnchorConstraints::size_pair_t right_anchor = AnchorConstraints::size_pair_t(r.get_startA(),r.get_startB());//dummy
    
    //std::cout << "right_anchor: "<<(right_anchor.first)<<","<<(right_anchor.second)<<std::endl;
    
    for (size_type i=r.get_startA(); i<=r.get_endA(); i++) {
	Fs[E_NO_NO]=infty_score_t::neg_infty;

	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(i));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(i));

	for (size_type j=min_col; j<=max_col; j++) {
	    
	    M(i,j) = align_noex( E_NO_NO,r.get_startA()-1,r.get_startB()-1,i,j, sv);	      
	    //
	    // score can be 0 (= drop prefix alignment) only if this is allowed due to constraints
	    if ( i<left_anchor.first && j<left_anchor.second ) {
	      M(i,j) = std::max( (infty_score_t)0, M(i,j) );
	    }
	    
	    if (i>=right_anchor.first && j>=right_anchor.second && max_score < M(i,j)) {
		max_score=M(i,j);
		max_i = i;
		max_j = j;
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
Aligner::align_top_level_localB() {
    // std::cout <<"align local B " << r.get_startA() << " " << r.get_startB() << " "
    //           << r.get_endA() << " " << r.get_endB() << std::endl;
    
    M_matrix_t &M=Ms[E_NO_NO];
    infty_score_t max_score=infty_score_t::neg_infty;
    init_state(E_NO_NO,r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1,true,false,false,false,def_scoring_view);
    
    for (size_type i=r.get_startA(); i<=r.get_endA(); i++) {

	// limit entries due to trace controller
	size_type min_col = std::max(r.get_startB(),params->trace_controller.min_col(i));
	size_type max_col = std::min(r.get_endB(),params->trace_controller.max_col(i));

	for (size_type j=min_col; j<=max_col; j++) {
	    M(i,j) = 
		max( (infty_score_t)0,
		     align_noex( E_NO_NO, 
				 r.get_startA()-1,r.get_startB()-1,i,j,def_scoring_view ) );
	    if (i==r.get_endA() && max_score < M(i,j)) {
		max_score=M(i,j);
		max_i = i;
		max_j = j;
	    }
	}
    }
    return max_score;
}


// compute the alignment score
infty_score_t
Aligner::align() {
    // ------------------------------------------------------------
    // computes D matrix (if not already done) and then does the alignment on the top level
    // ------------------------------------------------------------
    
    if (!D_created) align_D();

    if (params->SEQU_LOCAL) {
	return align_top_level_locally(def_scoring_view);
    } else { // sequence global alignment
	
	// align toplevel globally with potentially free endgaps (as given by description params->free_endgaps)
	return align_top_level_free_endgaps();
	
	// Without free end gaps, we could call (much simpler):
	// (we just call the more general method align_top_level_free_endgaps() since there is no performance penalty)
	//
	// align the subsequences that are specified by the restriction object
	/*
	align_in_arcmatch(r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1,false);
	max_i=r.get_endA();
	max_j=r.get_endB();
		
	return Ms[E_NO_NO](max_i,max_j);
	*/
    }
}

// ------------------------------------------------------------
// Aligner: traceback

void Aligner::trace_arcmatch(const ArcMatch &am) {
    
    const Arc &arcA=am.arcA();
    const Arc &arcB=am.arcB();
    
    size_type al=arcA.left();
    size_type ar=arcA.right();
    size_type bl=arcB.left();
    size_type br=arcB.right();
    
    // --------------------
    // handle case of stacking
    if ( scoring->stacking() ) {
	
	if (arc_matches.exists_inner_arc_match(am)) { 
	    const ArcMatch &inner_am = arc_matches.inner_arc_match(am);
	    
	    if (D(am) == D(inner_am) + scoring->arcmatch_stacked(am)) {
		
		const Arc & arcAI = inner_am.arcA();
		const Arc & arcBI = inner_am.arcB();
		
		alignment.add_basepairA(arcAI.left(),arcAI.right());
		alignment.add_basepairB(arcBI.left(),arcBI.right());
		alignment.append(arcAI.left(),arcBI.left());
		
		trace_arcmatch(inner_am);
		
		alignment.append(arcAI.right(),arcBI.right());
		
		return;
	    }
	}
    }
    
    // --------------------
    // now handle the case that arc match is not stacked
    
    // first recompute M
    align_in_arcmatch(al,ar, bl,br,
		      params->STRUCT_LOCAL);
    
    // then, trace in new M
    if (!params->STRUCT_LOCAL) {
	trace_in_arcmatch(0,al,ar-1,bl,br-1,false,def_scoring_view);
    } else {
	for (size_type k=0; k<4; k++) {
	    if (D(am) == Ms[k](ar-1,br-1) + scoring->arcmatch(am)) {
		trace_in_arcmatch(k,al,ar-1,bl,br-1,false,def_scoring_view);
		break;
	    }
	}
    }
    return;
}

void Aligner::trace_arcmatch_noLP(const ArcMatch &am) {
    
    assert(arc_matches.exists_inner_arc_match(am));
    
    const ArcMatch &inner_am = arc_matches.inner_arc_match(am);

    const Arc & arcAI = inner_am.arcA();
    const Arc & arcBI = inner_am.arcB();
    
    alignment.add_basepairA(arcAI.left(),arcAI.right());
    alignment.add_basepairB(arcBI.left(),arcBI.right());
    alignment.append(arcAI.left(),arcBI.left());
    
    if (D(am) == D(inner_am) + scoring->arcmatch(am,scoring->stacking())) {
	trace_arcmatch_noLP(inner_am);
    } else {
	// first recompute M
	align_in_arcmatch(arcAI.left(),arcAI.right(), arcBI.left(),arcBI.right(),
			  params->STRUCT_LOCAL);
	
	// then, trace in new M
	if (!params->STRUCT_LOCAL) {
	    trace_in_arcmatch(0,arcAI.left(),arcAI.right()-1,arcBI.left(),arcBI.right()-1,false,def_scoring_view);
	} else {
	    for (size_type k=0; k<4; k++) {
		if (D(am) == 
		    Ms[k](arcAI.right()-1,arcBI.right()-1)
		    + scoring->arcmatch(am,scoring->stacking()) + scoring->arcmatch(inner_am))
		{
		    trace_in_arcmatch(k,arcAI.left(),arcAI.right()-1,arcBI.left(),arcBI.right()-1,false,def_scoring_view);
		    break;
		}
	    }
	}
    }
    alignment.append(arcAI.right(),arcBI.right());
}

// trace and handle all cases that do not involve exclusions
template<class ScoringView>
void
Aligner::trace_noex(int state,size_type oal,size_type i,
		    size_type obl,size_type j,
		    bool tl,
		    ScoringView sv) {
    M_matrix_t &M=Ms[state];
    
    // determine where we get M(i,j) from
    
    // std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;
    
    // match
    if ( params->constraints.allowed_edge(i,j)
	 && M(i,j) ==  M(i-1,j-1)+sv.scoring()->basematch(i,j) ) {
	trace_in_arcmatch(state,oal,i-1,obl,j-1,tl,sv);
	alignment.append(i,j);
	return;
    }

    if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
	// del
	if ( (!params->constraints.aligned_in_a(i))
	     && M(i,j) == M(i-1,j)+sv.scoring()->gapA(i,j)) {
	    trace_in_arcmatch(state,oal,i-1,obl,j,tl,sv);
	    alignment.append(i,-1);
	    return;
	}
	// ins
	if ( (!params->constraints.aligned_in_b(j))
	     && M(i,j) == M(i,j-1)+sv.scoring()->gapB(i,j)) {
	    trace_in_arcmatch(state,oal,i,obl,j-1,tl,sv);
	    alignment.append(-1,j);
	    return;
	}
    } else { // base del and ins, affine cost
	// since, we didn't store all values in the E and F matrix
	// we do the traceback in linear time per entry
	// base del
	score_t gap_cost=sv.scoring()->indel_opening();
	for (size_type k=1;
	     (i >= oal+k)
		 && (! params->constraints.aligned_in_a(i-k+1));
	     k++)
	{
	    gap_cost += sv.scoring()->gapA(i-k+1,j);
	    
	    if ( M(i,j) == M(i-k,j) + gap_cost) {
		// gap in A of length k
		trace_in_arcmatch(state,oal,i-k,obl,j,tl,sv);
		for (size_type l=k;l>0;l--) {
		    alignment.append(i-l+1,-1);
		}
		return;
	    }
	}
	
	// base ins
	gap_cost=sv.scoring()->indel_opening();
	for (size_type k=1; 
	     (j >= obl+k)
		 && (! params->constraints.aligned_in_b(j-k+1));
	     k++)
	{
	    gap_cost += sv.scoring()->gapB(i,j-k+1);
	    if (M(i,j) == M(i,j-k) + gap_cost) {
		// gap in B of length k
		trace_in_arcmatch(state,oal,i,obl,j-k,tl,sv);
		for (size_type l=k;l>0;l--) {
		    alignment.append(-1,j-l+1);
		}
		
		return;
	    }
	}
    }

    // only consider arc match cases if edge (i,j) is allowed!
    if ( ! params->constraints.allowed_edge(i,j) ) return;
    
    // here (i,j) is allowed
    
    //  arc match
    const size_type &ar=i;
    const size_type &br=j;
    
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(ar,br).begin();
	arc_matches.common_right_end_list(ar,br).end() != it; ++it ) {
	
	// NOTES: *it is the arc match index
	//        we iterate only over valid arc matches, i.e.
	//        constraints (including anchor c. and heuristic ones) are satisified
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();
	
	if ( (arcA.left() <= oal) || (arcB.left() <= obl) ) continue;

        const size_type al=arcA.left();
	const size_type bl=arcB.left();
	    
	if ( M(i,j) == M(al-1,bl-1) + sv.D(am)) {
	    //
	    // do the trace for alignment left of the arc match
	    trace_in_arcmatch(state,oal,al-1,obl,bl-1,tl,sv);
	    
	    //cout << "arcmatch "<<(al)<<","<<i<<";"<<(bl)<<","<<j<<" :: "
	    //      <<(arcA->w)<<" + "<<(arcB->w)<< " + " << tau(al,bl,i,j)  <<endl;
	    
	    alignment.add_basepairA(al,ar);
	    alignment.add_basepairB(bl,br);
	    alignment.append(al,bl);
	    
	    // do the trace below the arc match
	    
	    if (params->no_lonely_pairs) {
		trace_arcmatch_noLP(am);
	    } else {
		trace_arcmatch(am);
	    }
	    
	    alignment.append(ar,br);
	    
	    return;
	}
    }
}

// do the trace within one arc match.
// the cases without exclusions are delegated to trace_noex
template <class ScoringView>
void
Aligner::trace_in_arcmatch(int state,int al,int i,int bl,int j,bool tl,ScoringView sv) {
    //pre: M matrices for arc computed
    M_matrix_t &M=Ms[state];

    /*
    string state_text[]={"E_NO_NO", "E_X_NO", "E_NO_X", "E_X_X",
      "E_OP_NO", "E_NO_OP", "E_OP_X", "E_X_OP"};
      cout << "trace_in_arcmatch "<<state_text[state]<<" al:"<<al<<" i:"<<i
      <<" bl:"<<bl<<" j:"<<j<<" :: "<< M(i,j) <<endl;
    */
    
    // terminate traceback if
    // * trace on toplevel
    // * sequence local, and
    // * entry == 0
    if (tl && params->SEQU_LOCAL && M(i,j)==(infty_score_t)0) {
	min_i=i;
	min_j=j;
	return;
    }

    if (i<=al) {
	if (state==E_NO_NO
	    || state==E_OP_NO
	    || state==E_X_NO) {
	    	    
	    // pad with gap edges, unless in special cases (local/semi-local alignment)
	    if (!(tl && (params->SEQU_LOCAL || params->free_endgaps.allow_left_1()))) {
		for (int k=bl+1;k<=j;k++) {
		    alignment.append(-1,k);
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
	    if (!(tl && ( params->SEQU_LOCAL || params->free_endgaps.allow_left_2()))) {
		for (int k=al+1;k<=i;k++) {
		    alignment.append(k,-1);
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
	} else if (M(i,j) == Ms[E_NO_NO](i,j))
	    trace_in_arcmatch(E_NO_NO,al,i,bl,j,tl,sv);
	break;
    case E_NO_OP:
	if (M(i,j) == M(i,j-1)) {
	    // cout << "exclude B "<<j<<endl;
	    trace_in_arcmatch(state,al,i,bl,j-1,tl,sv);
	} else if (M(i,j) == Ms[E_NO_NO](i,j))
	    trace_in_arcmatch(E_NO_NO,al,i,bl,j,tl,sv);
	break;
    case E_NO_X:
	if (M(i,j) == Ms[E_NO_OP](i,j)+sv.scoring()->exclusion()) {
	    trace_in_arcmatch(E_NO_OP,al,i,bl,j,tl,sv);
	} else
	    trace_noex(state,al,i,bl,j,tl,sv);
	break;
    case E_OP_X:
	if (M(i,j) == M(i-1,j)) {
	    //cout << "exclude A "<<i<<endl;
	    trace_in_arcmatch(state,al,i-1,bl,j,tl,sv);
	} else if (M(i,j) == Ms[E_NO_X](i,j))
	    trace_in_arcmatch(E_NO_X,al,i,bl,j,tl,sv);
	break;
    case E_X_NO:
	if (M(i,j) == Ms[E_OP_NO](i,j)+sv.scoring()->exclusion()) {
	    trace_in_arcmatch(E_OP_NO,al,i,bl,j,tl,sv);
	} else
	    trace_noex(state,al,i,bl,j,tl,sv);
	break;
    case E_X_OP:
	if (M(i,j) == M(i,j-1)) {
	    // cout << "exclude B "<<j<<endl;
	    trace_in_arcmatch(state,al,i,bl,j-1,tl,sv);
	} else if (M(i,j) == Ms[E_X_NO](i,j))
	    trace_in_arcmatch(E_X_NO,al,i,bl,j,tl,sv);
	break;
    case E_X_X:
	if (M(i,j) == Ms[E_OP_X](i,j)+sv.scoring()->exclusion()) {
	    trace_in_arcmatch(E_OP_X,al,i,bl,j,tl,sv); 
	} else if (M(i,j) == Ms[E_X_OP](i,j)+sv.scoring()->exclusion()) {
	    trace_in_arcmatch(E_X_OP,al,i,bl,j,tl,sv);
	} else
	    trace_noex(state,al,i,bl,j,tl,sv);
	break;
    }
}


template<class ScoringView>
void
Aligner::trace(ScoringView sv) {
    // pre: last call align_in_arcmatch(r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1);
    //      or align_top_level_locally for SEQU_LOCAL alignent
    
    // reset the alignment strings (to empty strings)
    // such that they can be written again during the trace
    alignment.clear();
    
    trace_in_arcmatch(E_NO_NO,r.get_startA()-1,max_i,r.get_startB()-1,max_j,true,sv);
}

void
Aligner::trace() {trace(def_scoring_view);}

/* ------------------------------------------------------------
   Compute k-best alignments by intervall splitting
*/

template <class T>
class greater_second {
public:
    bool operator() (const T& a, const T&b) {
	return a.second < b.second;
    }
};

typedef std::pair<AlignerRestriction,infty_score_t> task_t;

//! enumerate suboptimal alignments local alignments
void Aligner::suboptimal(size_type k,
			 score_t threshold,
			 bool opt_normalized,
			 score_t normalized_L,
			 size_type output_width,
			 bool opt_verbose,
			 bool opt_local_output,
			 bool opt_pos_output,
			 bool opt_write_structure
			 ) {
    Aligner &a=*this;
    
    // compute alignment score for a
    infty_score_t a_score = (!opt_normalized)?a.align():a.normalized_align(normalized_L,false);
    
    // make priority queue tasks
    std::priority_queue<task_t, vector<task_t>, greater_second<task_t> > tasks;
    
    // put a into tasks
    tasks.push(task_t(a.get_restriction(),a_score));

    size_type i=1;
    
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
	
	alignment.write(cout, output_width, task_score, 
			opt_local_output, opt_pos_output, opt_write_structure);
	
	if (!opt_pos_output) std::cout << std::endl
				       << std::endl;
	
	if (i==k) break; // break if enough solutions generated
	
	// split the longer sequence according to local alignment
	size_type lenA=task_r.get_endA()-task_r.get_startA();
	size_type lenB=task_r.get_endB()-task_r.get_startB();
	
	// make two clones of AlignerRestriction task_r
	AlignerRestriction r1(task_r);
	AlignerRestriction r2(task_r);
	
	if (lenA>lenB) {
	    // split A
	    int splitA = (alignment.get_local_startA() + alignment.get_local_endA())/2;
	    if (opt_verbose) std::cout <<"Split A at "<<splitA<<std::endl;
	    r1.set_endA(splitA);
	    r2.set_startA(splitA);
	} else {
	    int splitB = (alignment.get_local_startB() + alignment.get_local_endB())/2;
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
    if (!D_created) align_D();
    
    mod_scoring=new Scoring(*scoring); // make mod_scoring point to a copy of scoring
    
    // Apply Dinkelbach's algorithm
    
    score_t new_lambda=0;
    score_t lambda=-1;

    size_type iteration=0;
    
    // iterate until convergence
    while ( lambda != new_lambda )
    {
	++iteration;
	if (opt_verbose) std::cout << "Perform Dinkelbach iteration "<<iteration<<std::endl;
	
	lambda=new_lambda;
		
	// modify the scoring by lambda
	mod_scoring->modify_by_parameter(lambda);
	mod_scoring_view.set_lambda(lambda);
	
	infty_score_t score = align_top_level_locally(mod_scoring_view);
	
	alignment.clear();
	
	// perform a traceback for normalized alignment
	trace(mod_scoring_view);
	
	// compute length (of alignment) as sum of lengths of
	// aligned subsequences from the trace
	size_type length=max_i-min_i+1+max_j-min_j+1;
	
	// get score for the best alignment in the modified problem
	// but for unmodified scoring. Because for each position,
	// lambda was subtracted, we can just add length*lambda!
	score += length*lambda;
	
	new_lambda = score.finite_value()/(length+L);
	
	if (opt_verbose) std::cout << "Score: "<<score<<" Length: "<<length<<" Normalized Score: "<<new_lambda<<std::endl;

	if (opt_verbose) alignment.write(std::cout,
					 120,
					 (infty_score_t)new_lambda,
					 true,
					 false,
					 true
					 );
	if (opt_verbose) std::cout<<endl;
    }
    return (infty_score_t)new_lambda;
}
