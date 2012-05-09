#include "aligner_n.hh"
#include "anchor_constraints.hh"
#include "trace_controller.hh"
// #include "d_matrix.hh"

#include <math.h>
#include <assert.h>

#include <queue>

#include <iostream>


using namespace std;

namespace LocARNA {






// ------------------------------------------------------------
// AlignerN: align / compute similarity
//

AlignerN::AlignerN(const AlignerN &a)
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
 IAmat(a.IAmat),
 IBmat(a.IBmat),
// IA(a.IA),
// IB(a.IB),
 Ms(a.Ms),
 min_i(a.min_i),
 min_j(a.min_j),
 max_i(a.max_i),
 max_j(a.max_j),
 D_created(a.D_created),
 alignment(a.alignment),
 def_scoring_view(this),
 mod_scoring_view(this)
{}

AlignerN::AlignerN(const Sequence &seqA_,
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

	Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
	Dmat.fill(infty_score_t::neg_infty);

//	IA.resize(seqA.length()+1);
//	IB.resize(seqB.length()+1);
	IAmat.resize(seqA.length()+1, bpsB.num_bps());
	IBmat.resize(bpsA.num_bps(), seqB.length()+1);

	for (size_t k=0; k<(params->STRUCT_LOCAL?8:1); k++) {
		Ms[k].resize(seqA.length()+1,seqB.length()+1);
	}

}



AlignerN::~AlignerN() {
	if (mod_scoring!=0) delete mod_scoring;
}

template<class ScoringView>
infty_score_t
AlignerN::compute_IA(pos_type al, Arc arcB, pos_type i, ScoringView sv) {

	// compute IA entry
	infty_score_t max_score = infty_score_t::neg_infty;
	if ( (! params->constraints.aligned_in_a(i)) ) {
		// due to constraints, i can be deleted
		max_score =
				std::max(	max_score,
						IAmat(i-1, arcB.idx()) + sv.scoring()->gapA(i, arcB.right() )); //TODO: Review: Check for matrix functionality!

		max_score = max_score.normalized_neg(); //TODO: normalization required?
	}

	const BasePairs::RightAdjList &adjlA = bpsA.right_adjlist(i); //TODO: Ask: right_adjlist(i) returns a sorted list??

	//---------------------cost computation for left sided gap blocks ----------------------------
	std::vector<infty_score_t> blockGapCostsA;
	blockGapCostsA.resize( i - al + 1, infty_score_t::neg_infty); //one additional element for al
	//
	infty_score_t gap_score = (infty_score_t)0;

	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type lastPos;
	for (lastPos = al + 1; lastPos <= i; lastPos++) { //TODO: to be optimized. can be integrated in the arc loop

		//todo: !Review: Shall we check this trace_controller for locarna_n?
		if (params->trace_controller.min_col(lastPos) > arcB.left()) break; // fill only as long as column bl is accessible, remaining elements have been initialized with neg_infinity

		if (!gap_score.is_neg_infty()) {
			if (params->constraints.aligned_in_a(lastPos)) {
				gap_score = infty_score_t::neg_infty;
			}
			else {
				gap_score += sv.scoring()->gapA( lastPos, arcB.left());
			}
		}
		blockGapCostsA[lastPos - al ] = gap_score;
	}
	//-------------------------------------------------------------------------------------

	for (BasePairs::RightAdjList::const_iterator arcA = adjlA.begin();
			arcA != adjlA.end() && arcA->left() > al  ; ++arcA) {
		infty_score_t new_score = blockGapCostsA[arcA->left() - al ]
		                                         + sv.scoring()->gapA (i, arcB.right())
		                                         + sv.D(*arcA, arcB ) //todo: ugly code!
		                                         + sv.scoring()->arcgap (*arcA);

		if (new_score > max_score) {
			max_score = new_score;
		}

	}
	//TODO: Review: Check for matrix functionality!
	return (max_score.normalized_neg()); //TODO: is normalized_neg? required ?
	//TODO: IB computation, To be done similar to  IA

}

template<class ScoringView>
infty_score_t
AlignerN::compute_IB(Arc arcA, pos_type bl, pos_type k, ScoringView sv) {

	// compute IA entry
	infty_score_t max_score = infty_score_t::neg_infty;
	if ( (! params->constraints.aligned_in_b(k)) ) {
		// due to constraints, i can be deleted
		max_score =
				std::max(	max_score,
						IBmat( arcA.idx(), k-1) + sv.scoring()->gapB(arcA.right(), k)); //TODO: Review: Check for matrix functionality!

		max_score = max_score.normalized_neg(); //TODO: normalization required?
	}

	const BasePairs::RightAdjList &adjlB = bpsB.right_adjlist(k); //TODO: Ask: right_adjlist(i) returns a sorted list??

	//---------------------cost computation for left sided gap blocks ----------------------------
	std::vector<infty_score_t> blockGapCostsB;
	blockGapCostsB.resize( k - bl + 1, infty_score_t::neg_infty); //one additional element for al
	//
	infty_score_t gap_score = (infty_score_t)0;

	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type lastPos;
	for (lastPos = bl + 1; lastPos <= k; lastPos++) { //TODO: to be optimized. can be integrated in the arc loop

		//todo: !Review: Shall we check this trace_controller for locarna_n?
		if (params->trace_controller.min_col(lastPos) > arcA.left()) break; // fill only as long as column bl is accessible, remaining elements have been initialized with neg_infinity

		if (!gap_score.is_neg_infty()) {
			if (params->constraints.aligned_in_b(lastPos)) {
				gap_score = infty_score_t::neg_infty;
			}
			else {
				gap_score += sv.scoring()->gapB( arcA.left(), lastPos);
			}
		}
		blockGapCostsB[lastPos - bl ] = gap_score;
	}
	//-------------------------------------------------------------------------------------

	for (BasePairs::RightAdjList::const_iterator arcB = adjlB.begin();
			arcB != adjlB.end() && arcB->left() > bl  ; ++arcB) {
		infty_score_t new_score = blockGapCostsB[arcB->left() - bl ]
		                                         + sv.scoring()->gapB (arcA.right(), k)
		                                         + sv.D(arcA, *arcB ) //todo: ugly code!
		                                         + sv.scoring()->arcgap (*arcB);

		if (new_score > max_score) {
			max_score = new_score;
		}

	}
	//TODO: Review: Check for matrix functionality!
	return (max_score.normalized_neg()); //TODO: is normalized_neg? required ?
	//TODO: IB computation, To be done similar to  IA

}

// standard cases in alignment: base match, base in/del, arc match
// (for structure local alignment this is extended by exclusion handling)
//
// if lonely basepairs are disallowed, there is special treatment
//
template<class ScoringView>
infty_score_t
AlignerN::compute_M_entry(int state, pos_type al, pos_type bl, pos_type i, pos_type j,ScoringView sv) {

	assert(0<=state && state<4);

	assert(params->trace_controller.is_valid(i,j)); //todo: trace_controller for locarna_n?

	M_matrix_t &M = Ms[state];

	infty_score_t max_score = infty_score_t::neg_infty;

	// base match
	if (params->constraints.allowed_edge(i,j)) { //todo: constraints for locarna_n
		max_score = M(i-1,j-1) + sv.scoring()->basematch(i,j);
	}


	// base del
	max_score = std::max(max_score, M(i-1,j) + sv.scoring()->gapA(i,j));

	// base ins
	max_score = std::max(max_score, M(i,j-1) + sv.scoring()->gapB(i,j));

	// arc match

	// standard case for arc match (without restriction to lonely pairs)
	//
	//TODO: Important: Instead of common_right_end_list, two iterative loops of right_adjlist could be used. That might be more efficient but arc match scores are not available for an arbitrary pair of arcs
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list( i, j ).begin();
			arc_matches.common_right_end_list(i,j).end() != it; ++it ) {

		const ArcMatch &am = arc_matches.arcmatch(*it);

		const Arc &arcX=am.arcA();
		const Arc &arcY=am.arcB();

		pos_type xl = arcX.left();
		pos_type yl = arcY.left();

		if ( xl > al && yl > bl)
		{
			infty_score_t new_score =  M(arcX.left(),arcY.left()-1) + sv.D(arcX,arcY) + sv.scoring()->arcmatch(am);
			if (new_score > max_score) {
				max_score=new_score;
			}
		}
	}

	/*
	if ( params->constraints.allowed_edge(i,j) ) { //todo: constraints for locarna_n
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
				infty_score_t new_score =  M(arcA->left()-1,arcB->left()-1) + sv.D(*arcA,*arcB) + scoring->arcmatch(am);

				//				infty_score_t new_score = M(arcA->left()-1,arcB->left()-1) + sv.scoring()->arcmatch(arcA) + sv.scoring()->arcmatch(arcB);

				if (new_score > max_score) {
					//std::cout << *arcA << "-"<< *arcB << ": "<<M(arcA->left()-1,arcB->left()-1)<<"+"<<D(*arcA,*arcB)<<"="<<new_score<<std::endl;
					max_score=new_score;
				}
			}
		}
	}
	 */


	return max_score.normalized_neg();



}


// generic initalization method.
//

//TODO: Implement AlignerN::init_state method!
// Initialization in case of local alignment 
template <class ScoringView>
void
AlignerN::init_state(int state, pos_type al, pos_type ar, pos_type bl, pos_type br,
		bool globalA, bool exclA,
		bool globalB, bool exclB,
		ScoringView sv) {

    // alignments that have empty subsequence in A (i=al) and
	// end with gap in alistr of B do not exist ==> -infty

	M_matrix_t &M = Ms[state];

	// al,bl can only be reached in states, where this is legal with cost 0 for empty alignment
	M(al,bl) = (infty_score_t)0;

	std::cout << "COL "<<bl<<" AL: "<<al<<" AR: "<<ar<<std::endl;

	// init first col bl
	//
	infty_score_t indel_score = (infty_score_t)0; //TODO: indel_opening score for locarna_n
	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type i;
	for (i=al+1; i<ar; i++) {

		if (params->trace_controller.min_col(i)>bl) break; // fill only as long as column bl is accessible //todo: trace_controller for locarna_n

		if (!indel_score.is_neg_infty()) {
			if (params->constraints.aligned_in_a(i)) {
				indel_score=infty_score_t::neg_infty;
			}
			else {
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
	indel_score=(infty_score_t)0;  //TODO: indel_opening score for locarna_n

	pos_type j;
	for (j=bl+1 ; j < std::min(br, params->trace_controller.max_col(al)+1) ; j++) { //todo: trace_controller for locarna_n
		if (!indel_score.is_neg_infty()) {
			if (params->constraints.aligned_in_b(j)) {
				indel_score=infty_score_t::neg_infty;
			}
			else {
				indel_score += sv.scoring()->gapB(al,j);
			}
		}
		M(al,j) = (infty_score_t)indel_score;
	}

}

void AlignerN::fill_IA_entries ( pos_type al, Arc arcB, pos_type x)
{

	IAmat(al, arcB.idx()) = (infty_score_t)0; //TODO: Ask reinitialize IA each time?

	for (pos_type i = al+1; i <= x; i++) { //todo: Evaluation order of IA & IB elements? (i <= x)?

		// limit entries due to trace controller
		//todo: do we need to limit IA entries due to trace controller?
		IAmat(i, arcB.idx()) = compute_IA(al, arcB, i, def_scoring_view); //TODO: Review: Check for matrix functionality!
	}

}


void AlignerN::fill_IB_entries ( Arc arcA, pos_type bl, pos_type x)
{

	IBmat(arcA.idx(), bl) = (infty_score_t)0; //TODO: Ask reinitialize IA each time?

	for (pos_type k = bl+1; k <= x; k++) { //todo: Evaluation order of IA & IB elements? (i <= x)?

		// limit entries due to trace controller
		//todo: do we need to limit IA entries due to trace controller?
		IBmat(arcA.idx(), k) = compute_IB(arcA, bl, k, def_scoring_view); //TODO: Review: Check for matrix functionality!
	}

}
void
AlignerN::align_in_arcmatch(pos_type al,pos_type ar,pos_type bl,pos_type br,
		bool allow_exclusion) {

	assert(br>0); // if br<=0 we run into trouble below when computing br-1




	init_state(E_NO_NO,al,ar,bl,br,true ,false,true ,false,def_scoring_view);




	for (pos_type i=al+1; i<ar; i++) {

		// limit entries due to trace controller
		//todo: do we need to limit M entries due to trace controller?
		pos_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
		pos_type max_col = std::min(br-1,params->trace_controller.max_col(i));

		for (pos_type j=min_col; j<=max_col; j++) {
			Ms[E_NO_NO](i,j) = compute_M_entry(E_NO_NO,al,bl,i,j,def_scoring_view); //TODO: not using normalization has been hard coded?!
		}
	}


	if (allow_exclusion) { //TODO support exclusions
		std::cerr << "ERROR Exclusions are not supported!" << std::endl;
		assert ( ! allow_exclusion );
	}
}

// compute the entries in the D matrix that
// can be computed from the matrix/matrices M, IA, IB
// for the subproblem al,bl,max_ar,max_br
//
// pre: M,IA,IB matrices are computed by a call to
void 
AlignerN::fill_D_entries(pos_type al, pos_type bl)
{
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin(); //TODO: MUST be updated! which positions of matrices D,IA,IB should be computed?
			arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {

		const ArcMatch &am = arc_matches.arcmatch(*it);

		const Arc &arcA=am.arcA();
		const Arc &arcB=am.arcB();


		pos_type ar = arcA.right();
		pos_type br = arcB.right();

		infty_score_t m= Ms[0](ar,br);
		infty_score_t ia= IAmat(ar,br);
		infty_score_t ib= IBmat(ar,br);

		if (params->STRUCT_LOCAL) {
			std::cerr << "ERROR Structural Locality not implemented!" << std::endl;
			assert (! params->STRUCT_LOCAL);
			}

		D(am) = std::max(m, ia);
		D(am) = std::max(D(am), ib );

		//std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << std::endl;

		if (scoring->stacking()) {
			std::cout << "Warning! stacking not implemented!" << std::endl;
		}
	}
}



// compute all entries D
void
AlignerN::align_D() {
	// for al in r.get_endA() .. r.get_startA
	for (pos_type al=r.get_endA()+1; al>r.get_startA(); ) { al--;

	const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(al);
	if ( adjlA.empty() )
		continue;

//	pos_type max_bl = std::min(r.get_endB(),params->trace_controller.max_col(al)); //TODO: check trace_controller or not? It seems not!
//	pos_type min_bl = std::max(r.get_startB(),params->trace_controller.min_col(al));

	pos_type max_bl = r.get_endB();
	pos_type min_bl = r.get_startB();

	// for bl in max_bl .. min_bl
	for (pos_type bl=max_bl+1; bl > min_bl;) { bl--;

//TODO: Take care of max ar, br!!
	pos_type max_ar=r.get_endA();	//todo: any  restriction (e.g. trace_controller) required?
	pos_type max_br=r.get_endB();

	align_in_arcmatch(al,max_ar,bl,max_br,params->STRUCT_LOCAL);

	const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(bl);
	const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(al);

	for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
			arcB != adjlB.end(); arcB++)
	{
		fill_IA_entries(al, *arcB, max_ar ); //todo: max_ar?
	}

	//TODO: Fill IB
	for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
			arcA != adjlA.end(); arcA++)
	{
		fill_IB_entries(*arcA, bl, max_br ); //todo: max_ar?
	}





	//std::cout << al << ","<<bl<<":"<<std::endl
	//	      << Ms[E_NO_NO] << std::endl;

	// ------------------------------------------------------------
	// fill D matrix entries
	//
	if (params->no_lonely_pairs) {
		std:cerr << "no_lonely_pairs not implemented!" << std::endl;
		assert(! params->no_lonely_pairs);
	} else {
		fill_D_entries(al,bl);
	}
	}
	}

	D_created=true; // now the matrix D is built up
}


// compute the alignment score
infty_score_t
AlignerN::align() {
	// ------------------------------------------------------------
	// computes D matrix (if not already done) and then does the alignment on the top level
	// ------------------------------------------------------------

	if (!D_created) align_D();

	if (params->SEQU_LOCAL) {
		std::cerr << "SEQU_LOCAL is not supported by locarna_n\n" << std::endl;
		assert (! params->SEQU_LOCAL );
	} else { // sequence global alignment

		// align toplevel globally with potentially free endgaps (as given by description params->free_endgaps)
		//return align_top_level_free_endgaps(); //TODO: implement align_top_level_free_endgaps()


		//TODO: Which one to be used at outermost? align M or align D? What has been used in the Sankoff?
		align_in_arcmatch(r.get_startA()-1, r.get_endA(), r.get_startB()-1, r.get_endB(), false);
		return Ms[E_NO_NO]( r.get_startA()-1, r.get_startB()-1); //no free end gaps
	}
}



// ------------------------------------------------------------

template <class ScoringView>
	void AlignerN::trace_IA (pos_type i, const Arc &arcB, ScoringView sv)
{

}
// AlignerN: traceback
//TODO: Implement traceback!
template<class ScoringView>
void AlignerN::trace_arcmatch(const ArcMatch &am, ScoringView sv) {

	// std::cout << "trace_arcmatch " << am.arcA() << " " << am.arcB() <<std::endl;
	assert(! params->STRUCT_LOCAL); //TODO: support structure local


	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();

	pos_type al=arcA.left();
	pos_type ar=arcA.right();
	pos_type bl=arcB.left();
	pos_type br=arcB.right();

	// --------------------
	// case of stacking: not supported
	if ( scoring->stacking() ) { //TODO: support stacking
		std::cerr << "aligner_n: stacking is not supported!" << std::endl;
	}

	// --------------------
	// now handle the case that arc match is not stacked

	//first compute IA
	fill_IA_entries(al, arcB, ar-1);  //todo: If a one-dimensional IAmat were used, re-fill IA_entries for traces
	if ( IA(ar-1, arcB ) == sv.D(am) )
	{
		trace_IA(ar-1, arcB, sv);
		return;
	}
	//TODO: trace_IB in the similar way

	// first recompute M
	align_in_arcmatch(al,ar, bl,br,	params->STRUCT_LOCAL);

	trace_in_arcmatch(0,al,ar-1,bl,br-1,false,def_scoring_view);

	return;
}
// trace and handle all cases that do not involve exclusions
template<class ScoringView>
void AlignerN::trace_noex(int state,pos_type oal,pos_type i, pos_type obl,pos_type j, bool tl, ScoringView sv) {

	assert (state == E_NO_NO);
	M_matrix_t &M=Ms[state];

	// determine where we get M(i,j) from

	// std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;

	// match
	if ( params->constraints.allowed_edge(i,j)
			&& params->trace_controller.is_valid(i-1,j-1) // todo: check trace_controller?
			&& M(i,j) ==  M(i-1,j-1) + sv.scoring()->basematch(i,j) ) {
		trace_in_arcmatch(state, oal, i-1, obl, j-1, tl, sv);
		alignment.append(i,j);
		return;
	}

	if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
		// del
		if ( (!params->constraints.aligned_in_a(i))
				&& params->trace_controller.is_valid(i-1,j) //todo: check trace_controller?
				&& M(i,j) == M(i-1,j) + sv.scoring()->gapA(i,j)) {
			trace_in_arcmatch(state, oal, i-1, obl, j, tl, sv);
			alignment.append(i,-1);
			return;
		}
		// ins
		if ( (!params->constraints.aligned_in_b(j))
				&& params->trace_controller.is_valid(i,j-1) //todo: check trace_controller?
				&& M(i,j) == M(i,j-1)+sv.scoring()->gapB(i,j)) {
			trace_in_arcmatch(state, oal, i, obl, j-1, tl, sv);
			alignment.append(-1,j);
			return;
		}
	} else {
		// base del and ins, affine cost
		std::cout << "aligner_n: affine gap cost is not supported" << std::endl; //TODO: Support affine gap cost
		assert ( sv.scoring()->indel_opening() == 0 );
	}

	// only consider arc match cases if edge (i,j) is allowed and valid!
	if ( ! (params->constraints.allowed_edge(i,j)
			&& params->trace_controller.is_valid(i-1,j-1)) //todo: check trace_controller?
	) return;


	// here (i,j) is allowed and valid

	//  arc match
	const pos_type &ar = i;
	const pos_type &br = j;

	for(ArcMatchIdxVec::const_iterator it = arc_matches.common_right_end_list(ar,br).begin();
			arc_matches.common_right_end_list(ar,br).end() != it; ++it ) {

		// NOTES: *it is the arc match index
		//        we iterate only over valid arc matches, i.e.
		//        constraints (including anchor c. and heuristic ones) are satisfied

		const ArcMatch &am = arc_matches.arcmatch(*it);

		const Arc &arcA = am.arcA();
		const Arc &arcB = am.arcB();

		if ( (arcA.left() <= oal) || (arcB.left() <= obl) ) continue;

		const pos_type al = arcA.left();
		const pos_type bl = arcB.left();

		if ( M(i,j) == M(al-1, bl-1) + sv.D(am) + sv.scoring()->arcmatch(am) ) {
			//
			// do the trace for alignment left of the arc match
			trace_in_arcmatch(state, oal, al-1, obl, bl-1, tl, sv);

			//cout << "arcmatch "<<(al)<<","<<i<<";"<<(bl)<<","<<j<<" :: "
			//      <<(arcA->w)<<" + "<<(arcB->w)<< " + " << tau(al,bl,i,j)  <<endl;

			alignment.add_basepairA(al,ar);
			alignment.add_basepairB(bl,br);
			alignment.append(al,bl);

			// do the trace below the arc match

			if (params->no_lonely_pairs) {
				std:cerr << "no_lonely_pairs not implemented!" << std::endl;
				assert(! params->no_lonely_pairs);
			} else {
				trace_arcmatch(am, sv);
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
AlignerN::trace_in_arcmatch(int state,int al,int i,int bl,int j,bool tl,ScoringView sv) {
	//pre: M matrices for arc computed
	M_matrix_t &M=Ms[state];
	assert (state == E_NO_NO);
	assert(! params->SEQU_LOCAL); //Local seq alignment not implemented yet.
	// string state_text[]={"E_NO_NO", "E_X_NO", "E_NO_X", "E_X_X",
	// 			 "E_OP_NO", "E_NO_OP", "E_OP_X", "E_X_OP"};
	// cout << "trace_in_arcmatch "<<state_text[state]<<" al:"<<al<<" i:"<<i
	// 	 <<" bl:"<<bl<<" j:"<<j<<" :: "<< M(i,j) <<endl;

	if (i<=al) {
		for (int k=bl+1;k<=j;k++) { //TODO: end gaps cost is not free
			alignment.append(-1,k);
		}
		return;
	}

	if (j<=bl) {
		for (int k=al+1;k<=i;k++) {
			alignment.append(k,-1); //TODO: end gaps cost is not free
		}
	}


	switch(state) {
	case E_NO_NO:
		trace_noex(state, al, i, bl, j, tl, sv);
		break;
	}
}


template<class ScoringView>
void
AlignerN::trace(ScoringView sv) {
	// pre: last call align_in_arcmatch(r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1);
	//      or align_top_level_locally for SEQU_LOCAL alignent

	// reset the alignment strings (to empty strings)
	// such that they can be written again during the trace
	alignment.clear();

	// free end gap version: trace_in_arcmatch(E_NO_NO,r.get_startA()-1,max_i,r.get_startB()-1,max_j,true,sv);
	trace_in_arcmatch(E_NO_NO, r.get_startA()-1, r.get_endA(), r.get_startB()-1, r.get_endB(), true, sv);

}

void
AlignerN::trace() {trace(def_scoring_view);}

/* ------------------------------------------------------------
   Compute k-best alignments by intervall splitting
 */


//! type of a task (used in computing k-best alignment)
typedef std::pair<AlignerRestriction,infty_score_t> task_t;

} //end namespace LocARNA
