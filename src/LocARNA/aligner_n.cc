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
 mapperA(a.mapperA),
 mapperB(a.mapperB),
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
	const SparsificationMapper &mapperA_,
	const SparsificationMapper &mapperB_,
	const ArcMatches &arc_matches_,
	const AlignerParams *ap_,
	const Scoring *s_
)
: scoring(s_),
  mod_scoring(0),
  params(ap_),
  seqA(seqA_), seqB(seqB_),
  mapperA(mapperA_),
  mapperB(mapperB_),
  arc_matches(arc_matches_),
  bpsA(arc_matches_.get_base_pairsA()),
  bpsB(arc_matches_.get_base_pairsB()),
  r(1,1,seqA_.length(),seqB_.length()),
  D_created(false),
  alignment(seqA_,seqB_),
  def_scoring_view(this),
  mod_scoring_view(this)
{
    assert(!params->STRUCT_LOCAL);
    Ms.resize(params->STRUCT_LOCAL?8:1);

    Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
    Dmat.fill(infty_score_t::neg_infty);

    IAmat.resize(mapperA.get_max_info_vec_size()+1, bpsB.num_bps());
    IBmat.resize(bpsA.num_bps(), mapperB.get_max_info_vec_size()+1);
        for (size_t k=0; k<(params->STRUCT_LOCAL?8:1); k++) {
	Ms[k].resize(mapperA.get_max_info_vec_size()+1,mapperB.get_max_info_vec_size()+1);
    }

}


//destructor
AlignerN::~AlignerN() {
    if (mod_scoring!=0) delete mod_scoring;
}

// Compute/Returns aligning to the gap, the sequence range between leftSide & rightSide, not including right/left side
template<class ScoringView>
infty_score_t AlignerN::getGapCostBetween( pos_type leftSide, pos_type rightSide, bool isA, ScoringView sv) //todo: Precompute the matrix?!
{
    if (trace_debugging_output)
	cout << "getGapCostBetween: leftSide:" <<  leftSide << " rightSide:" << rightSide << "isA:" << isA << endl;
    assert(leftSide < rightSide);
    if (leftSide >= rightSide)
	return infty_score_t::neg_infty;

    infty_score_t gap_score = (infty_score_t)0;
    for (pos_type lastPos = leftSide+1; lastPos < rightSide; lastPos++)
    {
	if ( (isA && params->constraints.aligned_in_a(lastPos))
		|| ( !isA && params->constraints.aligned_in_b(lastPos)) ) {
	    return infty_score_t::neg_infty;
	}
	else {
	    gap_score += sv.scoring()->gapX( lastPos, isA);
	}

    }
    return gap_score;
}

/*
template<class ScoringView>
void AlignerN::compute_gap_costs( pos_type xl, pos_type xr, const Arc& arcY, std::vector<infty_score_t> &blockGapCostsX, bool isA, ScoringView sv )
{
    //---------------------cost computation for left sided gap blocks ----------------------------
    blockGapCostsX.resize( xr - xl , infty_score_t::neg_infty); //one additional element for al
    //
    infty_score_t gap_score = (infty_score_t)0;

    // handling of anchor constraints:
    // anchored positions must not be excluded,
    // nor deleted

    pos_type lastPos;
    for (lastPos = xl + 1; lastPos < xr; lastPos++) { //TODO: to be optimized. can be integrated in the arc loop

	//		if (isA && params->trace_controller.min_col(lastPos) > arcY.left()) break; // fill only as long as column bl is accessible, remaining elements have been initialized with neg_infinity
	//		if ( !isA && params->trace_controller.max_col(arcY.left()) < lastPos ) break;


	if (!gap_score.is_neg_infty()) {
	    if ( (isA && params->constraints.aligned_in_a(lastPos))
		    || ( !isA && params->constraints.aligned_in_b(lastPos)) ) {
		gap_score = infty_score_t::neg_infty;
	    }
	    else {
		gap_score += sv.scoring()->gapX( lastPos, isA);
	    }
	}
	else
	{
	    break; //no more block of deletion/insertion is possible;
	}

	blockGapCostsX[lastPos - xl ] = gap_score;
    }
    //-------------------------------------------------------------------------------------
}
*/

// Compute an element of the matrix IA/IB
template<class ScoringView>
infty_score_t AlignerN::compute_IX(index_t xl, const Arc& arcY, matidx_t i_index, bool isA, ScoringView sv) {

    bool constraints_aligned_pos = false; //constraints are ignored, params->constraints.aligned_in_a(i_seq_pos);
    const BasePairs &bpsX = isA? bpsA : bpsB;
    const SparsificationMapper &mapperX = isA ? mapperA : mapperB;

    seq_pos_t i_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index);
    seq_pos_t i_prev_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index-1); //TODO: Check border i_index==1,0

    const ArcIdxVec &arcIdxVecX = mapperX.valid_arcs_right_adj(xl, i_index);

    infty_score_t max_score = infty_score_t::neg_infty;

    //base insertion/deletion
    if ( !constraints_aligned_pos  ) {
	infty_score_t ins_del_score = IX(i_index-1, arcY, isA) + sv.scoring()->gapX(i_seq_pos, isA) + getGapCostBetween(i_prev_seq_pos, i_seq_pos, isA, sv);

	max_score = std::max( max_score, ins_del_score);
    }

    //arc deletion + align left side of the arc to gap
    for (ArcIdxVec::const_iterator arcIdx = arcIdxVecX.begin(); arcIdx != arcIdxVecX.end(); ++arcIdx)
    {
	const Arc& arcX = bpsX.arc(*arcIdx);
	infty_score_t new_score =  sv.D(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA)
						    + getGapCostBetween(xl, arcX.left(), isA, sv);

	if (new_score > max_score) {
	    max_score = new_score;
	}
    }

    return max_score;
}

//Compute an entry of matrix M
template<class ScoringView>
infty_score_t
AlignerN::compute_M_entry(int state, index_t al, index_t bl, matidx_t i_index, matidx_t j_index, ScoringView sv) {

    assert(state == 0);
    if (trace_debugging_output)	std::cout << "compute_M_entry al: " << al << " bl: " << bl << " i: " << i_index << " j: " << j_index << std::endl;
    M_matrix_t &M = Ms[state];

    bool constraints_alowed_edge = true;// constraints are ignored,  params->constraints.allowed_edge(i_seq_pos, j_seq_pos)
    infty_score_t max_score = infty_score_t::neg_infty;

    //define variables for sequence positions
    seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
    seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);
    seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1); //TODO: Check border i_index==1,0
    seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1); //TODO: Check border j_index==1,0

    // base match
    if ( constraints_alowed_edge &&
	    mapperA.pos_unpaired(al, i_index) && mapperB.pos_unpaired(bl, j_index) ) {
	max_score = M(i_index-1, j_index-1) + sv.scoring()->basematch(i_seq_pos, j_seq_pos)
							+ getGapCostBetween(i_prev_seq_pos, i_seq_pos, true, sv) + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false, sv);//todo: precompute range gapcosts with memoru complexity O(n^2)?
    }

    // base del
    infty_score_t del_score = M(i_index-1, j_index) + sv.scoring()->gapA(i_seq_pos) + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true, sv);
    max_score = std::max(max_score,  del_score );

    // base ins
    infty_score_t ins_score = M(i_index,j_index-1) + sv.scoring()->gapB(j_seq_pos) + getGapCostBetween( j_prev_seq_pos, j_seq_pos, false, sv);
    max_score = std::max(max_score, ins_score);

    //list of valid arcs ending at i/j
    const ArcIdxVec& arcsA = mapperA.valid_arcs_right_adj(al, i_index);
    const ArcIdxVec& arcsB = mapperB.valid_arcs_right_adj(bl, j_index);

    // arc match
    for (ArcIdxVec::const_iterator arcAIdx = arcsA.begin(); arcAIdx != arcsA.end(); ++arcAIdx)
    {
	const Arc& arcA = bpsA.arc(*arcAIdx);

	matidx_t  arcA_left_index_before   = mapperA.first_valid_mat_pos_before(al, arcA.left());
	seq_pos_t arcA_left_seq_pos_before = mapperA.get_pos_in_seq_new(al, arcA_left_index_before);

	for (ArcIdxVec::const_iterator arcBIdx = arcsB.begin(); arcBIdx != arcsB.end(); ++arcBIdx)
	{
	    const Arc& arcB = bpsB.arc(*arcBIdx);
	    matidx_t arcB_left_index_before = mapperB.first_valid_mat_pos_before(bl, arcB.left());
	    seq_pos_t arcB_left_seq_pos_before = mapperB.get_pos_in_seq_new(bl, arcB_left_index_before);

	    infty_score_t arc_match_score =  M(arcA_left_index_before, arcB_left_index_before) + sv.D( arcA, arcB ) +  sv.scoring()->arcmatch(arcA, arcB) //toask: Should I also care about scoring scheme for IA,IB?
								    + getGapCostBetween( arcA_left_seq_pos_before, arcA.left(), true, sv)
								    + getGapCostBetween( arcB_left_seq_pos_before, arcB.left(), false, sv);
	    if (arc_match_score > max_score) {
		max_score = arc_match_score;
//		if (trace_debugging_output)	std::cout << "compute_M_entry arcs " << arcX << " , " << arcY << " new score: " << new_score << "arc match score: " << sv.scoring()->arcmatch(am) << std::endl;
	    }
	}
    }

    return max_score;
}


// initializing matrix M
//
template <class ScoringView>
void
AlignerN::init_M(int state, pos_type al, pos_type ar, pos_type bl, pos_type br, ScoringView sv) {

    // alignments that have empty subsequence in A (i=al) and
    // end with gap in alistr of B do not exist ==> -infty

    if (trace_debugging_output)
	std::cout << "init_state al: " << al << " bl: " << bl << " ar: " << ar << " br: " << br << std::endl;

    M_matrix_t &M = Ms[state];

    //empty sequences A,B
    M(0,0) = (infty_score_t)0;

    // init first column
    //
    infty_score_t indel_score = (infty_score_t)0; //tomark: indel_opening score for locarna_n
    for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {

	seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al,i_index);
	if (trace_debugging_output)
	    std::cout << "i_index:" << i_index << " i_seq_pos:" << i_seq_pos << std::endl;
	//		if (params->trace_controller.min_col(i)>bl) break; // no trace controller in this version
	//tocheck:toask: check alignment constraints in the invalid positions between valid gaps
	if (!indel_score.is_neg_infty()) { //checked for optimization
/*	    if (params->constraints.aligned_in_a(i_seq_pos) )
	    {
		indel_score=infty_score_t::neg_infty;
	    }
	    else */ {
		seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al,i_index-1);
		indel_score = indel_score + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true, sv) + sv.scoring()->gapA(i_seq_pos); //toask: infty_score_t operator+ overloading
	    }
	}
	M(i_index,0) = indel_score;
    }

    // init first row
    //
    indel_score=(infty_score_t)0;
    for (matidx_t j_index=1 ; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {
	seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl,j_index);
	if (!indel_score.is_neg_infty()) { //checked for optimization
	   /* if (params->constraints.aligned_in_b(j_seq_pos)) {
		indel_score=infty_score_t::neg_infty;
	    }
	    else*/ {
		seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl,j_index-1);

		indel_score = indel_score + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false, sv) + sv.scoring()->gapB(j_seq_pos); //toask: infty_score_t operator+ overloading
	    }
	}
	M(0,j_index) = indel_score;
    }

}

//fill IA entries for a column with fixed al, arcB
void AlignerN::fill_IA_entries ( pos_type al, Arc arcB, pos_type max_ar)
{
    IAmat(0, arcB.idx()) = infty_score_t::neg_infty;
    for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {

	IAmat(i_index, arcB.idx()) = compute_IX(al, arcB, i_index, true, def_scoring_view);
    }
    //	cout << "fill_IA_entries al: "<< al << " arcB.idx: " << arcB.idx() << " arcB.left: " << arcB.left() << " arcB.right: " << arcB.right() << " IAmat: " << std::endl << IAmat << std::endl;
}

//fill IB entries for a row with fixed arcA, bl
void AlignerN::fill_IB_entries ( Arc arcA, pos_type bl, pos_type max_br)
{
    IBmat(arcA.idx(), 0) = infty_score_t::neg_infty;

    for (pos_type j_index = 1; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {		// limit entries due to trace controll
	IBmat(arcA.idx(), j_index) = compute_IX(bl, arcA, j_index, false, def_scoring_view);
    }
    //	cout << "fill_IB_entries arcA: " << arcA << " bl: "<< bl <<  " IBmat: " << std::endl << IBmat << std::endl;
}

//compute/align matrix M
void AlignerN::align_M(pos_type al,pos_type ar,pos_type bl,pos_type br, bool allow_exclusion) {

    assert(br>0); //todo: adding appropriate assertions

    //initialize M
    init_M(E_NO_NO, al, ar, bl, br, def_scoring_view);

    if (trace_debugging_output)	cout << "init_M finished" << endl;

    //iterate through valid entries
    for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {
	/*
	    //tomark: constraints
	    // limit entries due to trace controller
		pos_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
		pos_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	 */
	for (matidx_t j_index = 1; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {
	    Ms[E_NO_NO](i_index,j_index) = compute_M_entry(E_NO_NO,al,bl,i_index,j_index,def_scoring_view); //toask: where should we care about non_default scoring views
	}
    }

    assert ( ! allow_exclusion );

    if (trace_debugging_output)
	std::cout << "align_M aligned M is :" << std::endl << Ms[E_NO_NO] << std::endl;
}

// compute the entries in the D matrix that
// can be computed from the matrix/matrices M, IA, IB
// for the subproblem al,bl,max_ar,max_br
// pre: M,IA,IB matrices are computed by a call to
void AlignerN::fill_D_entries(pos_type al, pos_type bl)
{
    if (trace_debugging_output)
	std::cout << "fill_D_entries al: " << al << " bl: " << bl << std::endl;

    UnmodifiedScoringViewN sv = def_scoring_view; //toask: where should we care about non_default scoring views

    //iterate through arcs begining at al,bl
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();  //tocheck:toask:todo: IMPORTANT! can we use arc_matches to get the common endlist?? arcA,arcB may not be matched!
	    arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {

	const ArcMatch &am = arc_matches.arcmatch(*it);

	const Arc &arcA=am.arcA();
	const Arc &arcB=am.arcB();


	//define variables for sequence positions & sparsed indices
	seq_pos_t ar_seq_pos = arcA.right();
	seq_pos_t br_seq_pos = arcB.right();

	if (trace_debugging_output)
	    cout << "arcA:" << arcA << " arcB:" << arcB << "ar_seq_pos" << ar_seq_pos << "br_seq_pos" << br_seq_pos << endl;


	matidx_t ar_prev_mat_idx_pos = mapperA.first_valid_mat_pos_before(al, ar_seq_pos);
	matidx_t br_prev_mat_idx_pos = mapperB.first_valid_mat_pos_before(bl, br_seq_pos);

	if (trace_debugging_output)
	    cout << "arcA:" << arcA << " arcB:" << arcB << " ar_prev_mat_idx_pos:" << ar_prev_mat_idx_pos << " br_prev_mat_idx_pos:" << br_prev_mat_idx_pos << endl;

	seq_pos_t ar_prev_seq_pos = mapperA.get_pos_in_seq_new(al, ar_prev_mat_idx_pos);
	infty_score_t jumpGapCostA = getGapCostBetween(ar_prev_seq_pos, ar_seq_pos, true, sv);

	seq_pos_t br_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, br_prev_mat_idx_pos);
	infty_score_t jumpGapCostB = getGapCostBetween(br_prev_seq_pos, br_seq_pos, false, sv);

	//M,IA,IB scores
	infty_score_t m= Ms[0](ar_prev_mat_idx_pos, br_prev_mat_idx_pos) + jumpGapCostA + jumpGapCostB;
	infty_score_t ia= IAmat(ar_prev_mat_idx_pos,arcB.idx()) + jumpGapCostA;
	infty_score_t ib= IBmat(arcA.idx(),br_prev_mat_idx_pos) + jumpGapCostB;

	if (trace_debugging_output)	cout << "m=" << m << " ia=" << ia << " ib=" << ib << endl;

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
    for (pos_type al=r.get_endA()+1; al>r.get_startA(); ) {
	al--;
	if (trace_debugging_output) std::cout << "align_D al: " << al << std::endl;

	const BasePairs::LeftAdjList &adjlA = bpsA.left_adjlist(al);
	if ( adjlA.empty() )
	{
	    if (trace_debugging_output)	std::cout << "empty left_adjlist(al=)" << al << endl;
	    continue;
	}
	//	pos_type max_bl = std::min(r.get_endB(),params->trace_controller.max_col(al)); //tomark: trace_controller
	//	pos_type min_bl = std::max(r.get_startB(),params->trace_controller.min_col(al));

	pos_type max_bl = r.get_endB();
	pos_type min_bl = r.get_startB();

	// for bl in max_bl .. min_bl
	for (pos_type bl=max_bl+1; bl > min_bl;) {
	    bl--;

	    const BasePairs::LeftAdjList &adjlB = bpsB.left_adjlist(bl);

	    if ( adjlB.empty() )
	    {
		if (trace_debugging_output)	std::cout << "empty left_adjlist(bl=)" << bl << endl;
		continue;
	    }

	    // ------------------------------------------------------------
	    // old code for finding maximum arc ends:
	    
	    // pos_type max_ar=adjlA.begin()->right();	//tracecontroller not considered
	    // pos_type max_br=adjlB.begin()->right();

	    // //find the rightmost possible basepair for left base al
	    // for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
	    // 	 arcA != adjlA.end(); arcA++)
	    // 	{
	    // 	    if (max_ar < arcA->right() )
	    // 		max_ar = arcA->right();
	    // 	}
	    // //find the rightmost possible basepair for left base bl
	    // for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
	    // 	 arcB != adjlB.end(); arcB++)
	    // 	{
	    // 	    if (max_br < arcB->right() )
	    // 		max_br = arcB->right();
	    // 	}
	    
	    // ------------------------------------------------------------
	    // from aligner.cc: find maximum arc ends
	    pos_type max_ar=al;
	    pos_type max_br=bl;
	    
	    // get the maximal right ends of any arc match with left ends (al,bl)
	    // in noLP mode, we don't consider cases without immediately enclosing arc match
	    assert(params->no_lonely_pairs==false);
	    arc_matches.get_max_right_ends(al,bl,&max_ar,&max_br,params->no_lonely_pairs);
	    
	     // check whether there is an arc match at all
	    if (al==max_ar || bl == max_br) continue;


	    //compute matrix M
	    align_M(al,max_ar,bl,max_br,params->STRUCT_LOCAL);

	    //compute IA
	    for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
		    arcB != adjlB.end(); arcB++)
	    {
		fill_IA_entries(al, *arcB, max_ar );
	    }

	    //comput IB
	    for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
		    arcA != adjlA.end(); arcA++)
	    {
		fill_IB_entries(*arcA, bl, max_br );
	    }

	    // ------------------------------------------------------------
	    // now fill matrix D entries
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
	//return align_top_level_free_endgaps(); //tomark: implement align_top_level_free_endgaps()
	seq_pos_t ps_al = r.get_startA()-1;
	seq_pos_t ps_ar = r.get_endA()+1;
	seq_pos_t ps_bl = r.get_startB()-1;
	seq_pos_t ps_br = r.get_endB()+1;
	matidx_t last_index_A = mapperA.number_of_valid_mat_pos(ps_al)-1;
	seq_pos_t last_valid_seq_pos_A = mapperA.get_pos_in_seq_new(ps_al, last_index_A);
	matidx_t last_index_B = mapperB.number_of_valid_mat_pos(ps_bl)-1;
	seq_pos_t last_valid_seq_pos_B = mapperB.get_pos_in_seq_new(ps_bl, last_index_B);
	if(trace_debugging_output)
	    cout << "Align top level with ps_al, last_index_A, ps_bl, last_index_B," << ps_al << last_index_A << ps_bl << last_index_B << endl;
	align_M(ps_al, last_index_A, ps_bl, last_index_B, false); //tocheck: always use get_startA-1 (not zero) in sparsification_mapper and other parts
	return Ms[E_NO_NO]( last_index_A, last_index_B)
						+ getGapCostBetween( last_valid_seq_pos_A, ps_ar, true, def_scoring_view)  //toask: where should we care about non_default scoring views
						+ getGapCostBetween( last_valid_seq_pos_B, ps_br, false, def_scoring_view) ; //no free end gaps
    }
}

// ------------------------------------------------------------

template <class ScoringView>
void AlignerN::trace_IX (pos_type xl, matidx_t i_index, const Arc &arcY, bool isA, ScoringView sv)
{
    if (trace_debugging_output) std::cout << "****trace_IX****" << (isA?"A ":"B ") << " (" << xl << ","<< i_index << "] , " << arcY << std::endl;
    const BasePairs &bpsX = isA? bpsA : bpsB;
    const SparsificationMapper &mapperX = isA ? mapperA : mapperB;
    bool constraints_aligned_pos = false;

    seq_pos_t i_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index);
    seq_pos_t i_prev_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index-1);


    if ( i_seq_pos <= xl )
    {
	for (size_type k = arcY.left()+1; k < arcY.right(); k++) {
	    if (isA)
		alignment.append(-2,k);
	    else
		alignment.append(k,-2);
	}
	return;
    }

    if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost

	if( !constraints_aligned_pos
		&&
		IX(i_index, arcY, isA) == IX(i_index-1, arcY, isA) + sv.scoring()->gapX(i_seq_pos, isA) + getGapCostBetween(i_prev_seq_pos, i_seq_pos, isA, sv) )
	{
	    trace_IX( xl, i_index-1, arcY, isA, sv);
	    for ( size_type k = i_prev_seq_pos + 1; k <= i_seq_pos; k++)
	    {
		if (isA)
		    alignment.append(k, -2);
		else
		    alignment.append(-2, k);
	    }
	    return;
	}
    }else {
	// base del and ins, affine cost tbd
	assert ( sv.scoring()->indel_opening() == 0 );
    }

    const ArcIdxVec &arcIdxVecX = mapperX.valid_arcs_right_adj(xl, i_index);

    for (ArcIdxVec::const_iterator arcIdx = arcIdxVecX.begin(); arcIdx != arcIdxVecX.end(); ++arcIdx)
    {
	const Arc& arcX = bpsX.arc(*arcIdx);
	infty_score_t current_score =  sv.D(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA) + getGapCostBetween(xl, arcX.left(), isA, sv);

	if ( IX(i_index, arcY, isA) == current_score) {

	    if (trace_debugging_output) std::cout << "Arc Deletion for X " << (isA?"A ":"B ") << std::endl;
	    if (isA)
	    {
		alignment.add_basepairA(arcX.left(), arcX.right());
		for (size_type k = xl+1; k <= arcX.left(); k++) {
		    alignment.append(k, -2);
		}

		trace_D(arcX, arcY, sv);

		alignment.append(arcX.right(), -2);
	    }
	    else
	    {
		alignment.add_basepairB(arcX.left(), arcX.right());
		for (size_type k = xl+1; k <= arcX.left(); k++) {
		    alignment.append(-2, k);
		}

		trace_D(arcY, arcX, sv);

		alignment.append(-2, arcX.right());

	    }
	}

    }

}


// AlignerN: traceback
template<class ScoringView>
void AlignerN::trace_D(const Arc &arcA, const Arc &arcB, ScoringView sv) {

    if (trace_debugging_output) std::cout << "****trace_D****" << arcA << " " << arcB <<std::endl;
    assert(! params->STRUCT_LOCAL);

    seq_pos_t al = arcA.left();
    seq_pos_t ar_seq_pos = arcA.right();
    seq_pos_t bl = arcB.left();
    seq_pos_t br_seq_pos = arcB.right();
    seq_pos_t ar_prev_mat_idx_pos = mapperA.first_valid_mat_pos_before(al, ar_seq_pos);
    seq_pos_t ar_prev_seq_pos = mapperA.get_pos_in_seq_new(al, ar_prev_mat_idx_pos);
    infty_score_t jumpGapCostA = getGapCostBetween(ar_prev_seq_pos, ar_seq_pos, true, sv);

    matidx_t br_prev_mat_idx_pos = mapperB.first_valid_mat_pos_before(bl, br_seq_pos); //tocheck: ar or ar-1?
    seq_pos_t br_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, br_prev_mat_idx_pos);
    infty_score_t jumpGapCostB = getGapCostBetween(br_prev_seq_pos, br_seq_pos, false, sv);

    // --------------------
    // case of stacking: not supported
    assert(! scoring->stacking());

    // --------------------
    // now handle the case that arc match is not stacked

    //first compute IA
    fill_IA_entries(al, arcB, ar_seq_pos);
    if ( sv.D(arcA, arcB) == IA( ar_prev_mat_idx_pos, arcB ) + jumpGapCostA )
    {
	trace_IX(al, ar_prev_mat_idx_pos, arcB, true, sv);
	for ( size_type k = ar_prev_seq_pos + 1; k < ar_seq_pos; k++)
	{
	    alignment.append(k, -1);
	}
	return;
    }

    fill_IB_entries(arcA, bl, br_seq_pos);
    if (sv.D(arcA, arcB) ==  IB(arcA, br_prev_mat_idx_pos ) + jumpGapCostB )
    {
	trace_IX(bl, br_prev_mat_idx_pos, arcA, false, sv);
	for ( size_type k = br_prev_seq_pos + 1; k < br_seq_pos; k++)
	{
	    alignment.append(-1, k);
	}
	return;
    }

    // first recompute M
    align_M(al,ar_seq_pos, bl, br_seq_pos,	params->STRUCT_LOCAL);

    trace_M(0, al, ar_prev_mat_idx_pos, bl, br_prev_mat_idx_pos, false, def_scoring_view);
    for ( size_type k = ar_prev_seq_pos + 1; k < ar_seq_pos; k++)
    {
	alignment.append(k, -1);
    }
    for ( size_type k = br_prev_seq_pos + 1; k < br_seq_pos; k++)
    {
	alignment.append(-1, k);
    }

    return;
}

template<class ScoringView>
void AlignerN::trace_D(const ArcMatch &am, ScoringView sv) {

    trace_D(am.arcA(), am.arcB(), sv);
}

// trace and handle all cases that do not involve exclusions
template<class ScoringView>
void AlignerN::trace_M_noex(int state,pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool tl, ScoringView sv) {

    assert (state == E_NO_NO);
    M_matrix_t &M=Ms[state];


    seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
    seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);
    if ( i_seq_pos == al && j_seq_pos == bl )
    	return;


    seq_pos_t i_prev_seq_pos = 0;
    if ( i_seq_pos > al )
	i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1); //TODO: Check border i_index==1,0
    seq_pos_t j_prev_seq_pos = 0;
    if (j_seq_pos > bl )
	j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1); //TODO: Check border j_index==1,0
    bool constraints_alowed_edge = true; // constraints are not considered,  params->constraints.allowed_edge(i_seq_pos, j_seq_pos)
    bool constraints_aligned_pos_A = false; // TOcheck: Probably unnecessary, constraints are not considered
    bool constraints_aligned_pos_B = false; // TOcheck: Probably unnecessaryconstraints are not considered
    // determine where we get M(i,j) from



    // std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;

    // match

    if (  i_seq_pos > al && j_seq_pos > bl &&
	    constraints_alowed_edge
	    && mapperA.pos_unpaired(al, i_index) && mapperB.pos_unpaired(bl, j_index)
	    && M(i_index,j_index) ==  M(i_index-1, j_index-1) + sv.scoring()->basematch(i_seq_pos, j_seq_pos)
	    + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true, sv) + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false, sv) ) {
	if (trace_debugging_output) std::cout << "base match " << i_index << " , " << j_index << std::endl;

	trace_M(state, al, i_index-1, bl, j_index-1, tl, sv);

/*	for ( size_type k = i_prev_seq_pos + 1; k < i_seq_pos; k++)
	{
	    alignment.append(k, -1);
	}
	for ( size_type k = j_prev_seq_pos + 1; k < j_seq_pos; k++)
	{
	    alignment.append(-1, k);
	}
	*/
	alignment.append(i_seq_pos,j_seq_pos);
	return;
    }

    if ( sv.scoring()->indel_opening() == 0 ) { // base del and ins, linear cost
	// del
	if (  i_seq_pos > al &&
		!constraints_aligned_pos_A
		&& M(i_index,j_index) ==
			M(i_index-1, j_index) + sv.scoring()->gapA(i_seq_pos) + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true, sv) )
	{
	    trace_M(state, al, i_index-1, bl, j_index, tl, sv);
	    alignment.append(i_seq_pos, -1 );
	    /* for ( size_type k = i_prev_seq_pos + 1; k <= i_seq_pos; k++)
	    {
		alignment.append(k, -1);
	    }*/
	    return;
	}

	// ins
	if (  j_seq_pos > bl &&
		!constraints_aligned_pos_B
		&& M(i_index,j_index) ==
			M(i_index,j_index-1) + sv.scoring()->gapB(j_seq_pos) + getGapCostBetween( j_prev_seq_pos, j_seq_pos, false, sv))
	{
	    trace_M(state, al, i_index, bl, j_index-1, tl, sv);
	    alignment.append(-1, j_seq_pos);
	    /*for ( size_type k = j_prev_seq_pos + 1; k <= j_seq_pos; k++)
	    {
		alignment.append(-1, k);
	    }*/
	    return;
	}

    }
    else {
	// base del and ins, affine cost
	assert ( sv.scoring()->indel_opening() == 0 );
    }

    // only consider arc match cases if edge (i,j) is allowed and valid! (assumed valid)
    if ( ! constraints_alowed_edge  )
    {
	std::cerr << "WARNING: unallowed edge" << std::endl;
	return;
    }


    // here (i,j) is allowed and valid,

    //  arc match

    const ArcIdxVec& arcsA = mapperA.valid_arcs_right_adj(al, i_index);
        const ArcIdxVec& arcsB = mapperB.valid_arcs_right_adj(bl, j_index);
        for (ArcIdxVec::const_iterator arcAIdx = arcsA.begin(); arcAIdx != arcsA.end(); ++arcAIdx)
        {

            const Arc& arcA = bpsA.arc(*arcAIdx);

            matidx_t arcA_left_index_before = mapperA.first_valid_mat_pos_before(al, arcA.left());
            seq_pos_t arcA_left_seq_pos_before = mapperA.get_pos_in_seq_new(al, arcA_left_index_before);
            for (ArcIdxVec::const_iterator arcBIdx = arcsB.begin(); arcBIdx != arcsB.end(); ++arcBIdx)
            {
        	const Arc& arcB = bpsB.arc(*arcBIdx);
        	matidx_t arcB_left_index_before = mapperB.first_valid_mat_pos_before(bl, arcB.left());
        	seq_pos_t arcB_left_seq_pos_before = mapperB.get_pos_in_seq_new(bl, arcB_left_index_before);


        	infty_score_t loop_match_score =  M(arcA_left_index_before, arcB_left_index_before) + sv.D( arcA, arcB ) +  sv.scoring()->arcmatch(arcA, arcB) //toask: Should I also care about scoring scheme for IA,IB?
            									    + getGapCostBetween( arcA_left_seq_pos_before, arcA.left(), true, sv)
            									    + getGapCostBetween( arcB_left_seq_pos_before, arcB.left(), false, sv);
        	if ( M(i_index, j_index) == loop_match_score ) {

        	    if (trace_debugging_output) std::cout << "arcmatch "<< arcA <<";"<< arcB << " :: "   << std::endl;

        	    trace_M(state, al, arcA_left_index_before, bl, arcB_left_index_before, tl, sv);

        	    /*for ( size_type k = arcA_left_seq_pos_before + 1; k < arcA.left(); k++)
        	    {
        		alignment.append(k, -1);
        	    }
        	    for ( size_type k = arcB_left_seq_pos_before + 1; k < arcB.left(); k++)
        	    {
        		alignment.append(-1, k);
        	    }
        	    */
        	    alignment.add_basepairA(arcA.left(), arcA.right());
        	    alignment.add_basepairB(arcB.left(), arcB.right());
        	    alignment.append(arcA.left(),arcB.left());

        	    // do the trace below the arc match

        	    assert(! params->no_lonely_pairs);
        	    trace_D(arcA, arcB, sv);
        	    alignment.append(arcA.right(),arcB.right());



        	    return;

        	}
            }
        }



    if (trace_debugging_output) std::cout << "WARNING: No trace found!" << std::endl;
}

// do the trace within one arc match.
// the cases without exclusions are delegated to trace_noex
template <class ScoringView>
void
AlignerN::trace_M(int state,pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool tl, ScoringView sv) {
    //pre: M matrices for arc computed
    M_matrix_t &M=Ms[state];
    assert (state == E_NO_NO);
    assert(! params->SEQU_LOCAL); //Local seq alignment not implemented yet.

    seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
    seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);
    if (trace_debugging_output) std::cout << "******trace_M***** " << " al:" << al << " i:" << i_seq_pos <<" bl:"<< bl << " j:" << j_seq_pos << " :: " <<  M(i_index,j_index) << std::endl;

//    if ( i_seq_pos <= al ) {
//	for (int k = bl+1; k <= j_seq_pos; k++) { //TODO: end gaps cost is not free
//
//	    if ( ((al == r.get_startA()-1) && mapperB.is_valid_pos_external(k))
//		    || ( (al != r.get_startA()-1) && mapperB.is_valid_pos((k)) )
//		alignment.append(-1,k);
//
//
//	}
//    }

//    if (j_seq_pos <= bl) {
//	for (int k = al+1;k <= i_seq_pos; k++) {
//
//	    if ( ((bl == r.get_startB()-1) && mapperA.is_valid_pos_external(k))
//		    || ( (bl != r.get_startB()-1) && mapperA.is_valid_pos(k)) )
//		alignment.append(k,-1);
//
//	}
//	return;
//    }


    switch(state) {
    case E_NO_NO:
	trace_M_noex(state, al, i_index, bl, j_index, tl, sv);
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
    seq_pos_t ps_al = r.get_startA() - 1;
    matidx_t last_mat_idx_pos_A = mapperA.number_of_valid_mat_pos(ps_al) -1;//tocheck: check the correctness
    seq_pos_t last_seq_pos_A = mapperA.get_pos_in_seq_new(ps_al, last_mat_idx_pos_A);

    seq_pos_t ps_bl = r.get_startB() - 1;
    matidx_t last_mat_idx_pos_B = mapperB.number_of_valid_mat_pos(ps_bl) -1;//tocheck: check the correctness
    seq_pos_t last_seq_pos_B = mapperB.get_pos_in_seq_new(ps_bl, last_mat_idx_pos_B);

    trace_M(E_NO_NO, ps_al, last_mat_idx_pos_A, ps_bl, last_mat_idx_pos_B, true, sv); //TODO: right side for trace_M differs with align_M
/*    for ( size_type k = last_seq_pos_A + 1; k <= r.get_endA(); k++)//tocheck: check the correctness
    {
		alignment.append(k, -1);
    }
    for ( size_type k = last_seq_pos_B + 1; k <= r.get_endB(); k++)//tocheck: check the correctness
    {
	    alignment.append(-1, k);
    }
*/
}

void
AlignerN::trace() {
    stopwatch.start("trace");

    trace(def_scoring_view);

    stopwatch.stop("trace");
}


//! type of a task (used in computing k-best alignment)
typedef std::pair<AlignerRestriction,infty_score_t> task_t;

} //end namespace LocARNA
