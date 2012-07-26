#include "aligner_n.hh"
#include "anchor_constraints.hh"
#include "trace_controller.hh"
// #include "d_matrix.hh"

#include <math.h>
#include <assert.h>

#include <queue>

#include <iostream>


//using namespace std;

namespace LocARNA {

bool trace_debugging_output=false;


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

	IAmat.resize(seqA.length()+1, bpsB.num_bps());
	IBmat.resize(bpsA.num_bps(), seqB.length()+1);

	for (size_t k=0; k<(params->STRUCT_LOCAL?8:1); k++) {
		Ms[k].resize(seqA.length()+1,seqB.length()+1);
	}

}


//destructor
AlignerN::~AlignerN() {
	if (mod_scoring!=0) delete mod_scoring;
}

template<class ScoringView>
infty_score_t
AlignerN::compute_IX(pos_type xl, const Arc& arcY, pos_type i, bool isA, ScoringView sv) {

	bool constraints_aligned_pos;
	const BasePairs &bpsX = isA? bpsA : bpsB;
	const BasePairs::RightAdjList &adjlX = bpsX.right_adjlist(i);//TODO: Ask: right_adjlist(i) returns a sorted list??

	if ( isA )
	{
	    constraints_aligned_pos = params->constraints.aligned_in_a(i);
	}
	else
	{
	    constraints_aligned_pos = params->constraints.aligned_in_b(i);
	}
	// compute IA entry
	infty_score_t max_score = infty_score_t::neg_infty;
	if ( !constraints_aligned_pos  ) {
		// due to constraints, i can be deleted
		max_score =
				std::max(	max_score,
						IX(i-1, arcY, isA) + sv.scoring()->gapX(i, arcY.right(), isA )); //TODO: Be Careful: in gapX first and second variables does NOT indicate seqA or seqB

	}


	//---------------------cost computation for left sided gap blocks ----------------------------
	//TODO: blockGapCosts can be computed just once and reused in method calls, variable can be transformed to an upper level attribute
	std::vector<infty_score_t> blockGapCostsX;
	blockGapCostsX.resize( i - xl + 1, infty_score_t::neg_infty); //one additional element for al
	//
	infty_score_t gap_score = (infty_score_t)0;

	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type lastPos;
	for (lastPos = xl + 1; lastPos <= i; lastPos++) { //TODO: to be optimized. can be integrated in the arc loop

		if (isA && params->trace_controller.min_col(lastPos) > arcY.left()) break; // fill only as long as column bl is accessible, remaining elements have been initialized with neg_infinity
		if ( !isA && params->trace_controller.max_col(arcY.left()) < lastPos ) break;


		if (!gap_score.is_neg_infty()) {
			if ( (isA && params->constraints.aligned_in_a(lastPos))
					|| ( !isA && params->constraints.aligned_in_b(lastPos)) ) {
				gap_score = infty_score_t::neg_infty;
			}
			else {
				gap_score += sv.scoring()->gapX( lastPos, arcY.left(), isA);
			}
		}
		else
		{
			break; //no more block of deletion/insertion is possible;
		}

		blockGapCostsX[lastPos - xl ] = gap_score;
	}
	//-------------------------------------------------------------------------------------

	for (BasePairs::RightAdjList::const_iterator arcX = adjlX.begin();
			arcX != adjlX.end() && arcX->left() > xl  ; ++arcX) {
		infty_score_t new_score = blockGapCostsX[arcX->left() - xl ]
		                                         + sv.scoring()->gapX (i, arcY.right(), isA)
		                                         + sv.D(*arcX, arcY, isA ) //todo: ugly code!
		                                         + sv.scoring()->arcDel (*arcX, isA);

		if (new_score > max_score) {
			max_score = new_score;
		}

	}
	return max_score;

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
//	std::cout << "compute_M_entry al: " << al << " bl: " << bl << " i: " << i << " j: " << j << std::endl;
	M_matrix_t &M = Ms[state];

	infty_score_t max_score = infty_score_t::neg_infty;

	// base match
	if (params->constraints.allowed_edge(i,j)) {
		max_score = M(i-1,j-1) + sv.scoring()->basematch(i,j);
	}


	// base del
	max_score = std::max(max_score, M(i-1,j) + sv.scoring()->gapA(i,j));

	// base ins
	max_score = std::max(max_score, M(i,j-1) + sv.scoring()->gapB(i,j));

	// arc match

	// standard case for arc match (without restriction to lonely pairs)
	//
	//TODO: Instead of common_right_end_list, two iterative loops of right_adjlist could be used. That might be more efficient but arc match scores are not available for an arbitrary pair of arcs
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list( i, j ).begin();
			arc_matches.common_right_end_list(i,j).end() != it; ++it ) {

		const ArcMatch &am = arc_matches.arcmatch(*it);

		const Arc &arcX=am.arcA();
		const Arc &arcY=am.arcB();

		pos_type xl = arcX.left();
		pos_type yl = arcY.left();

		if ( xl > al && yl > bl)
		{
			infty_score_t new_score =  M(arcX.left()-1,arcY.left()-1) + sv.D(arcX,arcY) + sv.scoring()->arcmatch(am);
			if (new_score > max_score) {
				max_score=new_score;
//				std::cout << "compute_M_entry arcs " << arcX << " , " << arcY << " new score: " << new_score << "arc match score: " << sv.scoring()->arcmatch(am) << std::endl;
			}
		}
	}

	return max_score;



}


// generic initalization method.
//
template <class ScoringView>
void
AlignerN::init_state(int state, pos_type al, pos_type ar, pos_type bl, pos_type br, ScoringView sv) {

    // alignments that have empty subsequence in A (i=al) and
	// end with gap in alistr of B do not exist ==> -infty

	if (trace_debugging_output)
	    std::cout << "init_state al: " << al << " bl: " << bl << " ar: " << ar << " br: " << br << std::endl;
	M_matrix_t &M = Ms[state];

	// al,bl can only be reached in states, where this is legal with cost 0 for empty alignment
	M(al,bl) = (infty_score_t)0;

//	std::cout << "COL is"<<bl<<" AL: "<<al<<" AR: "<<ar<<std::endl;

	// init first col bl
	//
	infty_score_t indel_score = (infty_score_t)0; //TODO: indel_opening score for locarna_n
	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type i;
	for (i=al+1; i<ar; i++) {

		if (params->trace_controller.min_col(i)>bl) break; // fill only as long as column bl is accessible

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
	indel_score=(infty_score_t)0;

	pos_type j;
	for (j=bl+1 ; j < std::min(br, params->trace_controller.max_col(al)+1) ; j++) { // trace_controller for locarna_n
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

void AlignerN::fill_IA_entries ( pos_type al, Arc arcB, pos_type max_ar)
{

	IAmat(al, arcB.idx()) = infty_score_t::neg_infty;
	for (pos_type i = al+1; i < max_ar; i++) {

		// limit entries due to trace controller
		//IAmat(i, arcB.idx()) = compute_IA(al, arcB, i, def_scoring_view);

		IAmat(i, arcB.idx()) = compute_IX(al, arcB, i, true, def_scoring_view);
	}
//	cout << "fill_IA_entries al: "<< al << " arcB.idx: " << arcB.idx() << " arcB.left: " << arcB.left() << " arcB.right: " << arcB.right() << " IAmat: " << std::endl << IAmat << std::endl;


}



void AlignerN::fill_IB_entries ( Arc arcA, pos_type bl, pos_type max_br)
{

	IBmat(arcA.idx(), bl) = infty_score_t::neg_infty;

	for (pos_type k = bl+1; k < max_br; k++) {

		// limit entries due to trace controller
//		IBmat(arcA.idx(), k) = compute_IB(arcA, bl, k, def_scoring_view);
		IBmat(arcA.idx(), k) = compute_IX(bl, arcA, k, false, def_scoring_view);
	}
//	cout << "fill_IB_entries arcA: " << arcA << " bl: "<< bl <<  " IBmat: " << std::endl << IBmat << std::endl;

}
void
AlignerN::align_M(pos_type al,pos_type ar,pos_type bl,pos_type br,
		bool allow_exclusion) {

	assert(br>0); // if br<=0 we run into trouble below when computing br-1




	init_state(E_NO_NO, al, ar, bl, br, def_scoring_view);




	for (pos_type i=al+1; i<ar; i++) {

		// limit entries due to trace controller
		//need to limit M entries due to trace controller?
		pos_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
		pos_type max_col = std::min(br-1,params->trace_controller.max_col(i));

		for (pos_type j=min_col; j<=max_col; j++) {
			Ms[E_NO_NO](i,j) = compute_M_entry(E_NO_NO,al,bl,i,j,def_scoring_view);
		}
	}

	assert ( ! allow_exclusion );

//	std::cout << "align_M aligned M is :" << std::endl << Ms[E_NO_NO] << std::endl;
}

// compute the entries in the D matrix that
// can be computed from the matrix/matrices M, IA, IB
// for the subproblem al,bl,max_ar,max_br
//
// pre: M,IA,IB matrices are computed by a call to
void 
AlignerN::fill_D_entries(pos_type al, pos_type bl)
{
	if (trace_debugging_output) std::cout << "fill_D_entries al: " << al << " bl: " << bl << std::endl;
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();
			arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {

		const ArcMatch &am = arc_matches.arcmatch(*it);

		const Arc &arcA=am.arcA();
		const Arc &arcB=am.arcB();


		pos_type ar = arcA.right();
		pos_type br = arcB.right();

		infty_score_t m= Ms[0](ar-1, br-1);
		infty_score_t ia= IAmat(ar-1,arcB.idx());
		infty_score_t ib= IBmat(arcA.idx(),br-1);

		assert (! params->STRUCT_LOCAL);

		D(am) = std::max(m, ia);
		D(am) = std::max(D(am), ib );

//		std::cout <<"["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << std::endl;

		assert(! scoring->stacking());
	}
}



// compute all entries D
void
AlignerN::align_D() {
	// for al in r.get_endA() .. r.get_startA
	for (pos_type al=r.get_endA()+1; al>r.get_startA(); ) { al--;
//	std::cout << "align_D al: " << al << std:.endl;
	const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(al);
	if ( adjlA.empty() )
	{
//		std::cout << "empty left_adjlist(al)" << endl;
		continue;
	}
//	pos_type max_bl = std::min(r.get_endB(),params->trace_controller.max_col(al)); //TODO: check trace_controller or not? It seems not!
//	pos_type min_bl = std::max(r.get_startB(),params->trace_controller.min_col(al));

	pos_type max_bl = r.get_endB(); //tracecontroller not considered
	pos_type min_bl = r.get_startB(); //tracecontroller not considered

	// for bl in max_bl .. min_bl
	for (pos_type bl=max_bl+1; bl > min_bl;) { bl--;

	const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(bl);

	if ( adjlB.empty() )
	{
//		std::cout << "empty left_adjlist(bl)" << endl;
		continue;
	}
//	std::cout << "size of left_adjlist(bl): " << adjlB.size() << std::endl;
	pos_type max_ar=adjlA.begin()->right();	//tracecontroller not considered
	pos_type max_br=adjlB.begin()->right();

	//find the rightmost possible basepair for left base al
	for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
			arcA != adjlA.end(); arcA++)
	{
		if (max_ar < arcA->right() )
			max_ar = arcA->right();
	}
	//find the rightmost possible basepair for left base bl
	for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
			arcB != adjlB.end(); arcB++)
	{
		if (max_br < arcB->right() )
			max_br = arcB->right();
	}
	align_M(al,max_ar,bl,max_br,params->STRUCT_LOCAL);

	for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
			arcB != adjlB.end(); arcB++)
	{
		fill_IA_entries(al, *arcB, max_ar );
	}

	for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
			arcA != adjlA.end(); arcA++)
	{
		fill_IB_entries(*arcA, bl, max_br );
	}





	//std::cout << al << ","<<bl<<":"<<std::endl
	//	      << Ms[E_NO_NO] << std::endl;

	// ------------------------------------------------------------
	// fill D matrix entries
	//
		assert(! params->no_lonely_pairs);

		fill_D_entries(al,bl);

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
		exit(-1);
	} else { // sequence global alignment

		// align toplevel globally with potentially free endgaps (as given by description params->free_endgaps)
		//return align_top_level_free_endgaps(); //TODO: implement align_top_level_free_endgaps()

		align_M(r.get_startA()-1, r.get_endA()+1, r.get_startB()-1, r.get_endB()+1, false);
		return Ms[E_NO_NO]( r.get_endA(), r.get_endB()); //no free end gaps
	}
}

// ------------------------------------------------------------

template <class ScoringView>
	void AlignerN::trace_IX (pos_type xl, pos_type i, const Arc &arcY, bool isA, ScoringView sv)
{
	if (trace_debugging_output) std::cout << "****trace_IX****" << (isA?"A ":"B ") << " (" << xl << ","<< i << "] , " << arcY << std::endl;

	if ( i <= xl )
	{
		for (size_type k = arcY.left()+1; k < arcY.right(); k++) { //TODO: implement free end gaps cost
			if (isA)
				alignment.append(-1,k);
			else
				alignment.append(k,-1);
		}
		return;
	}

	if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
		if (isA)
		{
			if( !params->constraints.aligned_in_a(i)  //todo: is this check necessary? does scoring()->gapA take care of constraints?
					//todo: valid trace check for IA?
					&& IA(i, arcY) == IA(i-1, arcY) + sv.scoring()->gapA(i, arcY.right() ) )
			{
				trace_IX( xl, i-1, arcY, isA, sv);
				alignment.append(i, -1);
				return;
			}
		}
		else //isB
			if ( !params->constraints.aligned_in_b(i) && IB(arcY, i) == IB(arcY, i-1) + sv.scoring()->gapB(arcY.right(), i ) )
			{
				trace_IX( xl, i-1, arcY, isA, sv);
				alignment.append(-1, i);
				return;
			}

	}else {
		// base del and ins, affine cost
		assert ( sv.scoring()->indel_opening() == 0 );
	}


	const BasePairs::RightAdjList &adjlX = isA? bpsA.right_adjlist(i): bpsB.right_adjlist(i); //TODO: Ask: right_adjlist(i) returns a sorted list?

	//---------------------cost computation for left sided gap blocks ----------------------------
	std::vector<infty_score_t> blockGapCostsX;
	blockGapCostsX.resize( i - xl + 1, infty_score_t::neg_infty); //one additional element for al
	//
	infty_score_t gap_score = (infty_score_t)0;

	// handling of anchor constraints:
	// anchored positions must not be excluded,
	// nor deleted

	pos_type lastPos;
	for (lastPos = xl + 1; lastPos <= i; lastPos++) { //TODO: to be optimized. can be integrated in the arc loop

			if (isA && params->trace_controller.min_col(lastPos) > arcY.left()) break; // fill only as long as column bl is accessible, remaining elements have been initialized with neg_infinity
			if ( !isA && params->trace_controller.max_col(arcY.left()) < lastPos ) break; // todo: validate the code


			if (!gap_score.is_neg_infty()) {
				if ( (isA && params->constraints.aligned_in_a(lastPos))
						|| ( !isA && params->constraints.aligned_in_b(lastPos)) ) {
					gap_score = infty_score_t::neg_infty;
				}
				else {
					gap_score += sv.scoring()->gapX( lastPos, arcY.left(), isA);
				}
			}
			else
			{
				break; //no more block of deletion/insertion is possible;
			}

			blockGapCostsX[lastPos - xl ] = gap_score;
	}
	//-------------------------------------------------------------------------------------

	for (BasePairs::RightAdjList::const_iterator arcX = adjlX.begin();
			arcX != adjlX.end() && arcX->left() > xl  ; ++arcX) {

		if ( IX(i, arcY, isA) == blockGapCostsX[arcX->left() - xl ]
		                                        + sv.scoring()->gapX (i, arcY.right(), isA)
		                                        + sv.D(*arcX, arcY, isA ) //todo: ugly code!
		                                        + sv.scoring()->arcDel (*arcX, isA))
		{
			if (trace_debugging_output) std::cout << "Arc Deletion for X " << (isA?"A ":"B ") << std::endl;
			if (isA)
			{
				alignment.add_basepairA(arcX->left(), arcX->right());
				for (size_type k = xl+1; k <= arcX->left(); k++) { //TODO: end gaps cost is not free
					alignment.append(k, -1);
				}

				trace_D(*arcX, arcY, sv);

				alignment.append(arcX->right(), -1);
			}
			else
			{
				alignment.add_basepairB(arcX->left(), arcX->right());
				for (size_type k = xl+1; k <= arcX->left(); k++) { //TODO: end gaps cost is not free
					alignment.append(-1, k);
				}

				trace_D(arcY, *arcX, sv);

				alignment.append(-1, arcX->right());

			}
		}

	}


}


// AlignerN: traceback
template<class ScoringView>
void AlignerN::trace_D(const Arc &arcA, const Arc &arcB, ScoringView sv) {

	if (trace_debugging_output) std::cout << "****trace_D****" << arcA << " " << arcB <<std::endl;
	assert(! params->STRUCT_LOCAL);



	pos_type al=arcA.left();
	pos_type ar=arcA.right();
	pos_type bl=arcB.left();
	pos_type br=arcB.right();

	// --------------------
	// case of stacking: not supported
	assert(! scoring->stacking());

	// --------------------
	// now handle the case that arc match is not stacked

	//first compute IA
	fill_IA_entries(al, arcB, ar);
	if ( IA(ar-1, arcB ) == sv.D(arcA, arcB) )
	{
		trace_IX(al, ar-1, arcB, true, sv);
		return;
	}

	fill_IB_entries(arcA, bl, br);
	if ( IB(arcA, br-1 ) == sv.D(arcA, arcB) )
	{
		trace_IX(bl, br-1, arcA, false, sv);

		return;
	}

	// first recompute M
	align_M(al,ar, bl,br,	params->STRUCT_LOCAL);

	trace_M(0,al,ar-1,bl,br-1,false,def_scoring_view);

	return;
}

template<class ScoringView>
void AlignerN::trace_D(const ArcMatch &am, ScoringView sv) {

	trace_D(am.arcA(), am.arcB(), sv);
}

// trace and handle all cases that do not involve exclusions
template<class ScoringView>
void AlignerN::trace_M_noex(int state,pos_type oal,pos_type i, pos_type obl,pos_type j, bool tl, ScoringView sv) {

	assert (state == E_NO_NO);
	M_matrix_t &M=Ms[state];

	// determine where we get M(i,j) from

	// std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;

	// match
	if ( params->constraints.allowed_edge(i,j)
			&& params->trace_controller.is_valid(i-1,j-1) // todo: check trace_controller?
			&& M(i,j) ==  M(i-1,j-1) + sv.scoring()->basematch(i,j) ) {
		if (trace_debugging_output) std::cout << "base match " << i << " , " << j << std::endl;
		trace_M(state, oal, i-1, obl, j-1, tl, sv);
		alignment.append(i,j);
		return;
	}

	if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
		// del
		if ( (!params->constraints.aligned_in_a(i))
				&& params->trace_controller.is_valid(i-1,j) //todo: check trace_controller?
				&& M(i,j) == M(i-1,j) + sv.scoring()->gapA(i,j)) {
			trace_M(state, oal, i-1, obl, j, tl, sv);
			alignment.append(i,-1);
			return;
		}
		// ins
		if ( (!params->constraints.aligned_in_b(j))
				&& params->trace_controller.is_valid(i,j-1) //todo: check trace_controller?
				&& M(i,j) == M(i,j-1)+sv.scoring()->gapB(i,j)) {
			trace_M(state, oal, i, obl, j-1, tl, sv);
			alignment.append(-1,j);
			return;
		}
	} else {
		// base del and ins, affine cost
		assert ( sv.scoring()->indel_opening() == 0 );
	}

	// only consider arc match cases if edge (i,j) is allowed and valid!
	if ( ! (params->constraints.allowed_edge(i,j)
			&& params->trace_controller.is_valid(i-1,j-1)) //todo: check trace_controller?
	)
		{
			std::cerr << "WARNING: unallowed edge or trace is invalid" << std::endl;
			return;
		}


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
			trace_M(state, oal, al-1, obl, bl-1, tl, sv);

			if (trace_debugging_output) std::cout << "arcmatch "<<(al)<<","<<i<<";"<<(bl)<<","<<j<<" :: "   << std::endl;

			alignment.add_basepairA(al,ar);
			alignment.add_basepairB(bl,br);
			alignment.append(al,bl);

			// do the trace below the arc match

			assert(! params->no_lonely_pairs);
			trace_D(am, sv);
			alignment.append(ar,br);



			return;
		}
	}
	if (trace_debugging_output) std::cout << "WARNING: No trace found!" << std::endl;
}

// do the trace within one arc match.
// the cases without exclusions are delegated to trace_noex
template <class ScoringView>
void
AlignerN::trace_M(int state,int al,int i,int bl,int j,bool tl,ScoringView sv) {
	//pre: M matrices for arc computed
	M_matrix_t &M=Ms[state];
	assert (state == E_NO_NO);
	assert(! params->SEQU_LOCAL); //Local seq alignment not implemented yet.
	if (trace_debugging_output) std::cout << "******trace_M***** " << " al:" << al << " i:" << i <<" bl:"<< bl << " j:" << j << " :: " <<  M(i,j) << std::endl;

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
		return;
	}


	switch(state) {
	case E_NO_NO:
		trace_M_noex(state, al, i, bl, j, tl, sv);
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

	// free end gap version: trace_M(E_NO_NO,r.get_startA()-1,max_i,r.get_startB()-1,max_j,true,sv);
	trace_M(E_NO_NO, r.get_startA()-1, r.get_endA(), r.get_startB()-1, r.get_endB(), true, sv); //TODO: right side for trace_M differs with align_M

}

void
AlignerN::trace() {trace(def_scoring_view);}


//! type of a task (used in computing k-best alignment)
typedef std::pair<AlignerRestriction,infty_score_t> task_t;

} //end namespace LocARNA
