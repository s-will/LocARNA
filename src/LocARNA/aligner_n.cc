#include "aux.hh"
#include "global_stopwatch.hh"

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

    bool trace_debugging_output=false; //!< a static switch to enable generating debugging logs


    // ------------------------------------------------------------
    // AlignerN: align / compute similarity
    //

    AlignerN::AlignerN(const AlignerN &a)
	: params(a.params),
	  scoring(a.scoring),
	  mod_scoring(0),
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
	  IADmat(a.IADmat),
	  IBDmat(a.IBDmat),
	  Emat(a.Emat),
	  Fmat(a.Fmat),
	  Ms(a.Ms),
	  gapCostAmat(a.gapCostAmat),
	  gapCostBmat(a.gapCostBmat),
	  min_i(a.min_i),
	  min_j(a.min_j),
	  max_i(a.max_i),
	max_j(a.max_j),
	D_created(a.D_created),
	alignment(a.alignment),
	def_scoring_view(this),
	mod_scoring_view(this) 
    {}

    AlignerN::AlignerN(const AlignerParams &ap_)
	: params(dynamic_cast<const AlignerNParams *>(&ap_)),
	  scoring(params->scoring_),
	  mod_scoring(0),
	  seqA(*params->seqA_), seqB(*params->seqB_),
	  mapperA(*params->sparsification_mapperA_),
	  mapperB(*params->sparsification_mapperB_),
	  arc_matches(*params->arc_matches_),
	  bpsA(params->arc_matches_->get_base_pairsA()),
	  bpsB(params->arc_matches_->get_base_pairsB()),
	  r(1,1,params->seqA_->length(),params->seqB_->length()),
	  D_created(false),
	  alignment(*params->seqA_,*params->seqB_),
	def_scoring_view(this),
	mod_scoring_view(this)
    {
	assert(!params->struct_local_);
	Ms.resize(params->struct_local_?8:1);

	Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
	Dmat.fill(infty_score_t::neg_infty);

	IAmat.resize(mapperA.get_max_info_vec_size()+1, bpsB.num_bps());
	IBmat.resize(bpsA.num_bps(), mapperB.get_max_info_vec_size()+1);

	IADmat.resize(bpsA.num_bps(),bpsB.num_bps());
	IADmat.fill(infty_score_t::neg_infty);

	IBDmat.resize(bpsA.num_bps(),bpsB.num_bps());
	IBDmat.fill(infty_score_t::neg_infty);

	for (size_t k=0; k<(params->struct_local_?8:1); k++) {
	    Ms[k].resize(mapperA.get_max_info_vec_size()+1,mapperB.get_max_info_vec_size()+1);
	}
	Emat.resize(mapperA.get_max_info_vec_size()+1, mapperB.get_max_info_vec_size()+1);
	Fmat.resize(mapperA.get_max_info_vec_size()+1, mapperB.get_max_info_vec_size()+1);


	gapCostAmat.resize(seqA.length()+3, seqA.length()+3);
	gapCostBmat.resize(seqB.length()+3, seqB.length()+3);


    }


    //destructor
    AlignerN::~AlignerN() {
	if (mod_scoring!=0) delete mod_scoring;
    }

    // Computes and stores score of aligning a the subsequence between different possible leftSides & rightSides to the gap
    template <class ScoringView>
    void AlignerN::computeGapCosts(bool isA, ScoringView sv)
    {
	if (trace_debugging_output)	cout << "computeGapCosts " << (isA?'A':'B') << std::endl;
	const Sequence& seqX = isA?seqA:seqB;
	ScoreMatrix& gapCostXmat = isA?gapCostAmat:gapCostBmat;
	for( pos_type leftSide = 0;  leftSide <= seqX.length(); leftSide++)
	    {
		infty_score_t gap_score = (infty_score_t)0;
		gapCostXmat(leftSide, leftSide) = gap_score;
		for( pos_type lastPos = leftSide+1;  lastPos <= seqX.length(); lastPos++)
		    {
			if ( (isA && params->constraints_->aligned_in_a(lastPos))
			     || ( !isA && params->constraints_->aligned_in_b(lastPos)) ) {
			    gap_score = infty_score_t::neg_infty;
			    //		break;
			}
			else {
			    gap_score += sv.scoring()->gapX( lastPos, isA);
			}

			gapCostXmat(leftSide,lastPos+1) = gap_score;
		    }
		/*for (;lastPos <= seqX.length(); lastPos++)
		  {
		  gapCostXmat(leftSide,lastPos+1) = gap_score;
		  }*/
	    }
	if (trace_debugging_output)
	    cout << "computed computeGapCosts " << (isA?'A':'B') << std::endl;

    }
    // Returns score of aligning a the subsequence between leftSide & rightSide to the gap, not including right/left side
    inline
    infty_score_t AlignerN::getGapCostBetween( pos_type leftSide, pos_type rightSide, bool isA) //todo: Precompute the matrix?!
    {
	//    if (trace_debugging_output)	cout << "getGapCostBetween: leftSide:" <<  leftSide << " rightSide:" << rightSide << "isA:" << isA << endl;
	assert(leftSide < rightSide);

	return (isA?gapCostAmat(leftSide,rightSide):gapCostBmat(leftSide, rightSide));
    }


    // Compute an element of the matrix IA/IB
    template<class ScoringView>
    infty_score_t AlignerN::compute_IX(index_t xl, const Arc& arcY, matidx_t i_index, bool isA, ScoringView sv) {

	bool constraints_aligned_pos = false; //constraints are ignored, params->constraints_->aligned_in_a(i_seq_pos);
	const BasePairs &bpsX = isA? bpsA : bpsB;
	const SparsificationMapper &mapperX = isA ? mapperA : mapperB;

	seq_pos_t i_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index);
	seq_pos_t i_prev_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index-1); //TODO: Check border i_index==1,0

	const ArcIdxVec &arcIdxVecX = mapperX.valid_arcs_right_adj(xl, i_index);

	infty_score_t max_score = infty_score_t::neg_infty;

	//base insertion/deletion
	if ( !constraints_aligned_pos  ) {
	    infty_score_t gap_score =  getGapCostBetween(i_prev_seq_pos, i_seq_pos, isA)  + sv.scoring()->gapX(i_seq_pos, isA) ;
	    if (gap_score.is_finite())
		{    	// convert the base gap score to the loop gap score
		    gap_score  =	(infty_score_t)(sv.scoring()->loop_indel_score( gap_score.finite_value())); // todo: unclean interface and casting
		}
	    infty_score_t base_indel_score = IX(i_index-1, arcY, isA) + gap_score;
	    max_score = std::max( max_score, base_indel_score);
	}

	//arc deletion + align left side of the arc to gap
	for (ArcIdxVec::const_iterator arcIdx = arcIdxVecX.begin(); arcIdx != arcIdxVecX.end(); ++arcIdx)
	    {
		const Arc& arcX = bpsX.arc(*arcIdx);
		infty_score_t gap_score =  getGapCostBetween(xl, arcX.left(), isA);
		if (gap_score.is_finite())
		    {    // convert the base gap score to the loop gap score
			gap_score = (infty_score_t)(sv.scoring()->loop_indel_score( gap_score.finite_value()));
		    }
		infty_score_t arc_indel_score_extend = IXD(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA) + gap_score ;
		if (arc_indel_score_extend > max_score) {
		    max_score = arc_indel_score_extend;
		}

		infty_score_t arc_indel_score_open = sv.D(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA) + gap_score + sv.scoring()->indel_opening_loop();
		if (arc_indel_score_open > max_score) {
		    max_score = arc_indel_score_open;
		}
	    }

	return max_score;
    }


    //Compute an entry of matrix E
    template<class ScoringView>
    infty_score_t
    AlignerN::compute_E_entry(int state, index_t al, matidx_t i_index, matidx_t j_index, seq_pos_t i_seq_pos, seq_pos_t i_prev_seq_pos, ScoringView sv)
    {
	assert(state == 0);

	bool constraints_aligned_pos_A = false; // TOcheck: Probably unnecessary, constraints are not considered
	if (i_seq_pos <= al || constraints_aligned_pos_A) //check possibility of base deletion
	    return infty_score_t::neg_infty;

	// base del
	infty_score_t gap_cost = getGapCostBetween(i_prev_seq_pos, i_seq_pos, true) + sv.scoring()->gapA(i_seq_pos) ;
	infty_score_t extend_score = gap_cost + Emat(i_index-1,j_index) ;
	infty_score_t open_score = 	Ms[state](i_index-1, j_index) + gap_cost + sv.scoring()->indel_opening();
	return  (std::max(extend_score, open_score ));
    }


    //Compute an entry of matrix F
    template<class ScoringView>
    infty_score_t
    AlignerN::compute_F_entry(int state, index_t bl, matidx_t i_index, matidx_t j_index, seq_pos_t j_seq_pos, seq_pos_t j_prev_seq_pos, ScoringView sv)
    {
	assert(state == 0);

	bool constraints_aligned_pos_B = false; // TOcheck: Probably unnecessary, constraints are not considered
	if (j_seq_pos <= bl || constraints_aligned_pos_B) //check possibility of base deletion
	    return infty_score_t::neg_infty;

	// base ins
	infty_score_t gap_cost = getGapCostBetween(j_prev_seq_pos, j_seq_pos, false) + sv.scoring()->gapB(j_seq_pos);
	infty_score_t extend_score = Fmat(i_index,j_index-1) + gap_cost;
	infty_score_t open_score = 	Ms[state](i_index, j_index-1) + gap_cost  + sv.scoring()->indel_opening() ;
	return  (std::max(extend_score, open_score ));
    }

    //Compute an entry of matrix M
    template<class ScoringView>
    infty_score_t
    AlignerN::compute_M_entry(int state, index_t al, index_t bl, matidx_t i_index, matidx_t j_index, ScoringView sv) {

	assert(state == 0);
	M_matrix_t &M = Ms[state];

	bool constraints_alowed_edge = true;// constraints are ignored,  params->constraints_->allowed_edge(i_seq_pos, j_seq_pos)
	infty_score_t max_score = infty_score_t::neg_infty;

	//define variables for sequence positions
	seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
	seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);
	seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1); //TODO: Check border i_index==1,0
	seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1); //TODO: Check border j_index==1,0

	//    if (trace_debugging_output)	std::cout << "compute_M_entry al: " << al << " bl: " << bl << "i:" << i_seq_pos << "/i_index: " << i_index << " j: " << j_seq_pos << "/j_index:"<< j_index << std::endl;


	infty_score_t opening_cost_A;
	if (i_prev_seq_pos < (i_seq_pos - 1)) //implicit base deletion because of sparsification
	    opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
	else
	    opening_cost_A = (infty_score_t)0;

	infty_score_t opening_cost_B;
	if (j_prev_seq_pos < (j_seq_pos - 1)) //implicit base insertion because of sparsification
	    opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
	else
	    opening_cost_B  = (infty_score_t)0;


	// base match
	if ( constraints_alowed_edge &&
	     mapperA.pos_unpaired(al, i_index) && mapperB.pos_unpaired(bl, j_index) ) {
	    infty_score_t gap_match_score = getGapCostBetween(i_prev_seq_pos, i_seq_pos, true) + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false) + (sv.scoring()->basematch(i_seq_pos, j_seq_pos));
	    max_score = std::max( max_score,
				  (infty_score_t)(gap_match_score + opening_cost_B + Emat(i_index-1, j_index-1)) );
	    max_score = std::max( max_score,
				  (infty_score_t)(gap_match_score + opening_cost_A + Fmat(i_index-1, j_index-1)) );
	    max_score = std::max( max_score,
				  (infty_score_t)(gap_match_score + opening_cost_A + opening_cost_B + M(i_index-1, j_index-1)) );

	}

	// base del, for efficiency compute_E/F entry invoked within compute_M_entry
	Emat(i_index, j_index) = compute_E_entry(state, al, i_index, j_index, i_seq_pos, i_prev_seq_pos, sv);
	max_score = std::max(max_score,  Emat(i_index, j_index));

	// base ins
	Fmat(i_index, j_index) = compute_F_entry(state, bl, i_index, j_index, j_seq_pos, j_prev_seq_pos, sv);
	max_score = std::max(max_score,  Fmat(i_index, j_index));

	//list of valid arcs ending at i/j
	const ArcIdxVec& arcsA = mapperA.valid_arcs_right_adj(al, i_index);
	const ArcIdxVec& arcsB = mapperB.valid_arcs_right_adj(bl, j_index);

	// arc match
	for (ArcIdxVec::const_iterator arcAIdx = arcsA.begin(); arcAIdx != arcsA.end(); ++arcAIdx)
	    {
		const Arc& arcA = bpsA.arc(*arcAIdx);

		matidx_t  arcA_left_index_before   = mapperA.first_valid_mat_pos_before(al, arcA.left());
		seq_pos_t arcA_left_seq_pos_before = mapperA.get_pos_in_seq_new(al, arcA_left_index_before);
		if (arcA_left_seq_pos_before < (arcA.left() - 1)) //implicit base deletion because of sparsification
		    opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
		else
		    opening_cost_A = (infty_score_t)0;

		for (ArcIdxVec::const_iterator arcBIdx = arcsB.begin(); arcBIdx != arcsB.end(); ++arcBIdx)
		    {
			const Arc& arcB = bpsB.arc(*arcBIdx);


			matidx_t arcB_left_index_before = mapperB.first_valid_mat_pos_before(bl, arcB.left());
			seq_pos_t arcB_left_seq_pos_before = mapperB.get_pos_in_seq_new(bl, arcB_left_index_before);

			if (arcB_left_seq_pos_before < (arcB.left() - 1)) //implicit base insertion because of sparsification
			    opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
			else
			    opening_cost_B = (infty_score_t)0;
			if (trace_debugging_output)	std::cout << "\tmatching arcs: arcA" << arcA << "arcB:" << arcB << " D(arcA,arcB)=" << sv.D( arcA, arcB ) << " sv.scoring()->arcmatch(arcA, arcB)="<<sv.scoring()->arcmatch(arcA, arcB) << "M(" << arcA_left_index_before <<"," << arcB_left_index_before << ")=" << M(arcA_left_index_before, arcB_left_index_before) << std::endl;

			infty_score_t gap_match_score = getGapCostBetween( arcA_left_seq_pos_before, arcA.left(), true)	+ getGapCostBetween( arcB_left_seq_pos_before, arcB.left(), false)
			    + sv.D( arcA, arcB ) + sv.scoring()->arcmatch(arcA, arcB);
			infty_score_t arc_match_score =  gap_match_score  + opening_cost_A + opening_cost_B + M(arcA_left_index_before, arcB_left_index_before);
			//	    if (trace_debugging_output) std::cout << "gap_match_score:" << gap_match_score << std::endl;
			//	    if (trace_debugging_output) std::cout << "arc_match_score:" << arc_match_score << std::endl;

			arc_match_score = std::max( arc_match_score,
						    (infty_score_t)(gap_match_score + opening_cost_B + Emat (arcA_left_index_before, arcB_left_index_before)) );
			arc_match_score = std::max( arc_match_score,
						    (infty_score_t)(gap_match_score + opening_cost_A + Fmat (arcA_left_index_before, arcB_left_index_before)) ) ;
			if (arc_match_score > max_score) {
			    max_score = arc_match_score;
			    //		if (trace_debugging_output)	std::cout << "compute_M_entry arcs " << arcA << " , " << arcB << "arc match score: " << arc_match_score << std::endl;
			}
		    }
	    }

	return max_score;
    }


    // initializing matrix M
    //
    template <class ScoringView>
    void
    AlignerN::init_M_E_F(int state, pos_type al, pos_type ar, pos_type bl, pos_type br, ScoringView sv) {

	// alignments that have empty subsequence in A (i=al) and
	// end with gap in alistr of B do not exist ==> -infty

	if (trace_debugging_output)
	    std::cout << "init_state al: " << al << " bl: " << bl << " ar: " << ar << " br: " << br << std::endl;

	M_matrix_t &M = Ms[state];

	//empty sequences A,B
	M(0,0) = (infty_score_t)0;

	Emat(0,0) = infty_score_t::neg_infty;//tocheck:validity
	Fmat(0,0) = infty_score_t::neg_infty;//tocheck:validity


	// init first column
	//
	infty_score_t indel_score = (infty_score_t)(sv.scoring()->indel_opening());
	for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {

	    seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al,i_index);
	    if (trace_debugging_output)
		std::cout << "i_index:" << i_index << " i_seq_pos:" << i_seq_pos << std::endl;
	    //		if (params->trace_controller.min_col(i)>bl) break; // no trace controller in this version
	    //tocheck:toask: check alignment constraints in the invalid positions between valid gaps
	    if (!indel_score.is_neg_infty()) { //checked for optimization
		/*	    if (params->constraints_->aligned_in_a(i_seq_pos) )
			    {
			    indel_score=infty_score_t::neg_infty;
			    }
			    else */ {
		    seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al,i_index-1);
		    indel_score = indel_score + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true) + sv.scoring()->gapA(i_seq_pos);
		}
	    }
	    Emat(i_index, 0) = indel_score;
	    Fmat(i_index, 0) = infty_score_t::neg_infty;
	    M(i_index,0) = indel_score;//same as Emat(i_index, 0);

	}

	// init first row
	//
	indel_score = (infty_score_t)(sv.scoring()->indel_opening());
	for (matidx_t j_index=1 ; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {
	    seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl,j_index);
	    if (!indel_score.is_neg_infty()) { //checked for optimization
		/* if (params->constraints_->aligned_in_b(j_seq_pos)) {
		   indel_score=infty_score_t::neg_infty;
		   }
		   else*/ {
		    seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl,j_index-1);

		    indel_score = indel_score + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false) + sv.scoring()->gapB(j_seq_pos); //toask: infty_score_t operator+ overloading
		}
	    }
	    Emat(0,j_index) = infty_score_t::neg_infty;
	    Fmat(0,j_index) = indel_score;
	    M(0,j_index) = indel_score; // same as Fmat(0,j_index);

	}

    }

    //fill IA entries for a column with fixed al, arcB
    void AlignerN::fill_IA_entries ( pos_type al, Arc arcB, pos_type max_ar)
    {
	if (trace_debugging_output)
	    cout << "fill_IA_entries: " <<  "al=" << al << "max_ar=" << max_ar << ", arcB=" << arcB << endl;

	IAmat(0, arcB.idx()) = infty_score_t::neg_infty;
	for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {

	    IAmat(i_index, arcB.idx()) = compute_IX(al, arcB, i_index, true, def_scoring_view);

	    //fill IAD matrix entries //tocheck: verify
	    seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
	    seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1);
	    if (bpsA.exists_arc(al,i_seq_pos))
		{
		    const Arc& arcA = bpsA.arc(al, i_seq_pos);
		    IADmat(arcA.idx(), arcB.idx()) = ((IAmat(i_index-1,arcB.idx()))) + getGapCostBetween(i_prev_seq_pos, i_seq_pos, true);
		}
	}
	//	cout << "fill_IA_entries al: "<< al << " arcB.idx: " << arcB.idx() << " arcB.left: " << arcB.left() << " arcB.right: " << arcB.right() << " IAmat: " << std::endl << IAmat << std::endl;
    }

    //fill IB entries for a row with fixed arcA, bl
    void AlignerN::fill_IB_entries ( Arc arcA, pos_type bl, pos_type max_br)
    {
	if (trace_debugging_output)
	    cout << "fill_IB_entries: " << "arcA=" << arcA<< ", bl=" << bl << "max_br=" << max_br << endl;
	IBmat(arcA.idx(), 0) = infty_score_t::neg_infty;

	for (pos_type j_index = 1; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {		// limit entries due to trace control


	    IBmat(arcA.idx(), j_index) = compute_IX(bl, arcA, j_index, false, def_scoring_view);
	    //fill IBD matrix entries
	    seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);
	    //	cout << "j_seq_pos=" << j_seq_pos  << endl;
	    seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1);
	    //	cout <<" j_prev_seq_pos=" << j_prev_seq_pos << endl;
	    if (bpsB.exists_arc(bl,j_seq_pos))
		{
		    const Arc& arcB = bpsB.arc(bl,j_seq_pos);
		    if (trace_debugging_output)
			cout << "exists arcB" << arcB << "  current IBDmat(" << arcA.idx() << "," << arcB.idx()<<")=" << IBDmat(arcA.idx(), arcB.idx()) << endl;

		    IBDmat(arcA.idx(), arcB.idx()) = IBmat(arcA.idx(), j_index-1) + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false);
		    if (trace_debugging_output)
			cout << "IBDmat(" << arcA.idx() << "," << arcB.idx()<<")=" << IBDmat(arcA.idx(), arcB.idx()) << endl;
		}
	}
	//	cout << "fill_IB_entries arcA: " << arcA << " bl: "<< bl <<  " IBmat: " << std::endl << IBmat << std::endl;
    }

    //compute/align matrix M
    void AlignerN::fill_M_entries(pos_type al,pos_type ar,pos_type bl,pos_type br, bool allow_exclusion) {

	assert(br>0); //todo: adding appropriate assertions

	//initialize M
	init_M_E_F(E_NO_NO, al, ar, bl, br, def_scoring_view);

	if (trace_debugging_output)	cout << "init_M finished" << endl;
	//    if (al==0 && bl==0)	stopwatch.start("compute_m entries top level");
	//iterate through valid entries
	for (matidx_t i_index = 1; i_index < mapperA.number_of_valid_mat_pos(al); i_index++) {
	    /*
	    //tomark: constraints
	    // limit entries due to trace controller
	    pos_type min_col = std::max(bl+1,params->trace_controller.min_col(i));
	    pos_type max_col = std::min(br-1,params->trace_controller.max_col(i));
	    */
	    for (matidx_t j_index = 1; j_index < mapperB.number_of_valid_mat_pos(bl); j_index++) {
		// E and F matrix entries will be computed by compute_M_entry
		Ms[E_NO_NO](i_index,j_index) = compute_M_entry(E_NO_NO,al,bl,i_index,j_index,def_scoring_view); //toask: where should we care about non_default scoring views
		//	    if (trace_debugging_output) std::cout << "M["<< i_index << "," << j_index << "]=" << Ms[E_NO_NO](i_index,j_index) << std::endl;

	    }
	}
	//    if (al==0 && bl==0) 	stopwatch.stop("compute_m entries top level");


	assert ( ! allow_exclusion );

	//    if (trace_debugging_output)	std::cout << "align_M aligned M is :" << std::endl << Ms[E_NO_NO] << std::endl;
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

	//iterate through arcs beginning at al,bl
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(al,bl).begin();  //tocheck:toask:todo: IMPORTANT! can we use arc_matches to get the common endlist?? arcA,arcB may not be matched!
	    arc_matches.common_left_end_list(al,bl).end() != it; ++it ) {

	    const ArcMatch &am = arc_matches.arcmatch(*it);

	    const Arc &arcA=am.arcA();
	    const Arc &arcB=am.arcB();


	    //define variables for sequence positions & sparsed indices
	    seq_pos_t ar_seq_pos = arcA.right();
	    seq_pos_t br_seq_pos = arcB.right();

	    if (trace_debugging_output)
		cout << "arcA:" << arcA << " arcB:" << arcB  << endl;


	    matidx_t ar_prev_mat_idx_pos = mapperA.first_valid_mat_pos_before(al, ar_seq_pos);
	    matidx_t br_prev_mat_idx_pos = mapperB.first_valid_mat_pos_before(bl, br_seq_pos);


	    seq_pos_t ar_prev_seq_pos = mapperA.get_pos_in_seq_new(al, ar_prev_mat_idx_pos);
	    infty_score_t jumpGapCostA = getGapCostBetween(ar_prev_seq_pos, ar_seq_pos, true);

	    seq_pos_t br_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, br_prev_mat_idx_pos);
	    infty_score_t jumpGapCostB = getGapCostBetween(br_prev_seq_pos, br_seq_pos, false);

	    if (trace_debugging_output)
		{
		    cout << " ar_prev_mat_idx_pos:" << ar_prev_mat_idx_pos << " br_prev_mat_idx_pos:" << br_prev_mat_idx_pos << endl;
		    cout << " ar_prev_seq_pos:" << ar_prev_seq_pos << " br_prev_seq_pos:" << br_prev_seq_pos << endl;
		}

	    //M,IA,IB scores
	    //	infty_score_t m= Ms[0](ar_prev_mat_idx_pos, br_prev_mat_idx_pos) + jumpGapCostA + jumpGapCostB;

	    //-----three cases for gap extension/initiation ---

	    infty_score_t opening_cost_A;
	    if (ar_prev_seq_pos < (ar_seq_pos - 1)) //implicit base deletion because of sparsification
		opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
	    else
		opening_cost_A = (infty_score_t)0;

	    infty_score_t opening_cost_B;
	    if (br_prev_seq_pos < (br_seq_pos - 1)) //implicit base insertion because of sparsification
		opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
	    else
		opening_cost_B  = (infty_score_t)0;

	    infty_score_t gap_score = jumpGapCostA + jumpGapCostB;
	    infty_score_t mdel = (infty_score_t)(gap_score + opening_cost_B + Emat(ar_prev_mat_idx_pos, br_prev_mat_idx_pos)) ;
	    infty_score_t mins = (infty_score_t)(gap_score + opening_cost_A + Fmat(ar_prev_mat_idx_pos, br_prev_mat_idx_pos)) ;
	    infty_score_t mm = (infty_score_t)(gap_score + opening_cost_A + opening_cost_B +  Ms[0](ar_prev_mat_idx_pos, br_prev_mat_idx_pos) ) ;

	    if (trace_debugging_output)	cout << "mdel=" << mdel << " mins=" << mins << " mm=" << mm << endl;


	    infty_score_t m = std::max( mm, std::max(mdel, mins));

	    //------------------------------



	    infty_score_t ia= IAmat(ar_prev_mat_idx_pos,arcB.idx()) + jumpGapCostA;
	    infty_score_t ib= IBmat(arcA.idx(),br_prev_mat_idx_pos) + jumpGapCostB;

	    assert(IADmat(arcA.idx(),arcB.idx()) == infty_score_t::neg_infty || IADmat(arcA.idx(),arcB.idx()) == ia);
	    assert(IBDmat(arcA.idx(),arcB.idx()) == infty_score_t::neg_infty || IBDmat(arcA.idx(), arcB.idx()) == ib);

	    IADmat(arcA.idx(),arcB.idx()) = ia; //TODO: avoid recomputation
	    IBDmat(arcA.idx(),arcB.idx()) = ib; //TODO: avoid recomputation
	    if (trace_debugging_output)	cout << "m=" << m << " ia=" << ia << " ib=" << ib << endl;


	    //	assert(ia == iad);

	    //	cout << "IBDmat" << endl  << IBDmat << endl;
	    //	assert(ib == ibd);



	    assert (! params->struct_local_);

	    D(am) = std::max(m, ia);
	    D(am) = std::max(D(am), ib );

	    //			std::cout <<"D["<< am.arcA() << "," <<am.arcB() <<"]:" << D(am) << std::endl;

	    assert(! scoring->stacking());
	}
    }


    // compute all entries D
    void
    AlignerN::align_D() {
	computeGapCosts(true, def_scoring_view);//gap costs A //tocheck:always def_score view!
	computeGapCosts(false, def_scoring_view);//gap costs B //tocheck:always def_score view!

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
		assert(params->no_lonely_pairs_==false);
		arc_matches.get_max_right_ends(al,bl,&max_ar,&max_br,params->no_lonely_pairs_);
	    
		// check whether there is an arc match at all
		if (al==max_ar || bl == max_br) continue;


		//compute matrix M
		//	    stopwatch.start("compM");
		fill_M_entries(al,max_ar,bl,max_br,params->struct_local_);
		//	    stopwatch.stop("compM");


		//compute IA
		//	    stopwatch.start("compIA");
		for (BasePairs::LeftAdjList::const_iterator arcB = adjlB.begin();
		     arcB != adjlB.end(); arcB++)
		    {
			fill_IA_entries(al, *arcB, max_ar );
		    }
		//	    stopwatch.stop("compIA");

		//comput IB
		//	    stopwatch.start("compIB");
		for (BasePairs::LeftAdjList::const_iterator arcA = adjlA.begin();
		     arcA != adjlA.end(); arcA++)
		    {
			fill_IB_entries(*arcA, bl, max_br );
		    }
		//	    stopwatch.stop("compIB");


		// ------------------------------------------------------------
		// now fill matrix D entries
		//
		assert(! params->no_lonely_pairs_);
		fill_D_entries(al,bl);
	    }
	}
	if (trace_debugging_output) std::cout << "M matrix:" << std::endl << Ms[0] << std::endl;
	if (trace_debugging_output) std::cout << "D matrix:" << std::endl << Dmat << std::endl;

	D_created=true; // now the matrix D is built up
    }


    // compute the alignment score
    infty_score_t
    AlignerN::align() {
	// ------------------------------------------------------------
	// computes D matrix (if not already done) and then does the alignment on the top level
	// ------------------------------------------------------------
	if (!D_created)
	    {
		stopwatch.start("alignD");
		align_D();
		stopwatch.stop("alignD");
	    }

	if (params->sequ_local_) {
	    throw failure("sequ_local is not supported by locarna_n");
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
		cout << "Align top level with ps_al:" << ps_al << ", last_index_A:" << last_index_A << "/last_seq_posA:" << last_valid_seq_pos_A << ", ps_bl:" << ps_bl << ", last_index_B:" << last_index_B << "/last_seq_posB:" << last_valid_seq_pos_B <<endl;

	    stopwatch.start("align top level");
	    fill_M_entries(ps_al, last_index_A, ps_bl, last_index_B, false); //tocheck: always use get_startA-1 (not zero) in sparsification_mapper and other parts
	    stopwatch.stop("align top level");
	    if (trace_debugging_output) std::cout << "M matrix:" << std::endl << Ms[0] << std::endl;
	    if (trace_debugging_output) std::cout << "M(" << last_index_A << "," << last_index_B << ")=" << Ms[E_NO_NO]( last_index_A, last_index_B) << " getGapCostBetween are:"<< getGapCostBetween( last_valid_seq_pos_A, ps_ar, true) << std::endl;//"  "  << getGapCostBetween( last_valid_seq_pos_B, ps_br, false) << std::endl;
	    return Ms[E_NO_NO]( last_index_A, last_index_B)
		+ getGapCostBetween( last_valid_seq_pos_A, ps_ar, true)  //toask: where should we care about non_default scoring views
		+ getGapCostBetween( last_valid_seq_pos_B, ps_br, false) ; //no free end gaps
	}
    }

    // ------------------------------------------------------------

    template <class ScoringView>
    void AlignerN::trace_IX (pos_type xl, matidx_t i_index, const Arc &arcY, bool isA, ScoringView sv)
    {
	const BasePairs &bpsX = isA? bpsA : bpsB;
	const SparsificationMapper &mapperX = isA ? mapperA : mapperB;
	bool constraints_aligned_pos = false;

	seq_pos_t i_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index);
	seq_pos_t i_prev_seq_pos = mapperX.get_pos_in_seq_new(xl, i_index-1);

	if (trace_debugging_output) std::cout << "****trace_IX****" << (isA?"A ":"B ") << " (" << xl << ","<< i_seq_pos << "] , " << arcY << std::endl;


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

	// base del and ins

	if( !constraints_aligned_pos )
	    {
		infty_score_t gap_score =  getGapCostBetween(i_prev_seq_pos, i_seq_pos, isA)  + sv.scoring()->gapX(i_seq_pos, isA) ;
		if( gap_score.is_finite() )
		    {    	// convert the base gap score to the loop gap score
			gap_score  = (infty_score_t)(sv.scoring()->loop_indel_score( gap_score.finite_value())); // todo: unclean interface and casting
			if (IX(i_index, arcY, isA) == IX(i_index-1, arcY, isA) + gap_score )
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
		    }
	    }



	const ArcIdxVec &arcIdxVecX = mapperX.valid_arcs_right_adj(xl, i_index);

	for (ArcIdxVec::const_iterator arcIdx = arcIdxVecX.begin(); arcIdx != arcIdxVecX.end(); ++arcIdx)
	    {
		const Arc& arcX = bpsX.arc(*arcIdx);
		if (trace_debugging_output) std::cout << "arcX=" << arcX  << std::endl;

	
		infty_score_t gap_score =  getGapCostBetween(xl, arcX.left(), isA);
		if (gap_score.is_finite())
		    {    // convert the base gap score to the loop gap score
			gap_score = (infty_score_t)(sv.scoring()->loop_indel_score( gap_score.finite_value()));

			infty_score_t arc_indel_score_extend = IXD(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA) + gap_score ;

			if ( IX(i_index, arcY, isA) == arc_indel_score_extend) {

			    if (trace_debugging_output) std::cout << "Arc Deletion extension for X " << (isA?"A ":"B ") << "arcX=" << arcX << " arcY=" << arcY << std::endl;
			    if (isA)
				{
				    alignment.add_basepairA(arcX.left(), arcX.right());
				    for (size_type k = xl+1; k <= arcX.left(); k++) {
					alignment.append(k, -2);
				    }

				    trace_IXD(arcX, arcY, isA, sv);

				    alignment.append(arcX.right(), -2);
				}
			    else
				{
				    alignment.add_basepairB(arcX.left(), arcX.right());
				    for (size_type k = xl+1; k <= arcX.left(); k++) {
					alignment.append(-2, k);
				    }

				    trace_IXD(arcY, arcX, isA, sv);

				    alignment.append(-2, arcX.right());

				}
			    return;
			}

			infty_score_t arc_indel_score_open = sv.D(arcX, arcY, isA) + sv.scoring()->arcDel(arcX, isA) + gap_score + sv.scoring()->indel_opening_loop();

			if ( IX(i_index, arcY, isA) == arc_indel_score_open) {

			    if (trace_debugging_output) std::cout << "Arc Deletion opening for X " << (isA?"A ":"B ") << std::endl;
			    if (isA)
				{
				    alignment.add_deleted_basepairA(arcX.left(), arcX.right());
				    for (size_type k = xl+1; k <= arcX.left(); k++) {
					alignment.append(k, -2);
				    }

				    trace_D(arcX, arcY, sv);

				    alignment.append(arcX.right(), -2);
				}
			    else
				{
				    alignment.add_deleted_basepairB(arcX.left(), arcX.right());
				    for (size_type k = xl+1; k <= arcX.left(); k++) {
					alignment.append(-2, k);
				    }

				    trace_D(arcY, arcX, sv);

				    alignment.append(-2, arcX.right());

				}
			    return;
			}

		    }
	    }
	if (trace_debugging_output) std::cout << "WARNING: trace_IX No trace found!" << std::endl;

    }
    // AlignerN: traceback
    template<class ScoringView>
    void AlignerN::trace_IXD(const Arc &arcA, const Arc &arcB, bool isA, ScoringView sv) {

	if (trace_debugging_output) std::cout << "****trace_IXD****" << (isA?"A ":"B ") << arcA << " " << arcB <<std::endl;
	assert(! params->struct_local_);

	seq_pos_t al = arcA.left();
	seq_pos_t ar_seq_pos = arcA.right();
	seq_pos_t bl = arcB.left();
	seq_pos_t br_seq_pos = arcB.right();
	seq_pos_t ar_prev_mat_idx_pos = mapperA.first_valid_mat_pos_before(al, ar_seq_pos);
	seq_pos_t ar_prev_seq_pos = mapperA.get_pos_in_seq_new(al, ar_prev_mat_idx_pos);
	infty_score_t jumpGapCostA = getGapCostBetween(ar_prev_seq_pos, ar_seq_pos, true);

	matidx_t br_prev_mat_idx_pos = mapperB.first_valid_mat_pos_before(bl, br_seq_pos); //tocheck: ar or ar-1?
	seq_pos_t br_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, br_prev_mat_idx_pos);
	infty_score_t jumpGapCostB = getGapCostBetween(br_prev_seq_pos, br_seq_pos, false);

	// --------------------
	// case of stacking: not supported
	assert(! scoring->stacking());

	// --------------------
	// handle the case that arc match is not stacked

	if (isA)//trace IAD
	    {
		//first compute IA
		fill_IA_entries(al, arcB, ar_seq_pos);
		if ( IADmat(arcA.idx(), arcB.idx()) == IA( ar_prev_mat_idx_pos, arcB ) + jumpGapCostA )
		    {
			trace_IX(al, ar_prev_mat_idx_pos, arcB, true, sv);
			for ( size_type k = ar_prev_seq_pos + 1; k < ar_seq_pos; k++)
			    {
				alignment.append(k, -1);
			    }
			return;
		    }
	    }
	else //trace IBD
	    {
		fill_IB_entries(arcA, bl, br_seq_pos);
		if (trace_debugging_output)
		    cout << "IXD(" << arcA.idx() << "," << arcB.idx() << ")="  << IBDmat(arcA.idx(), arcB.idx()) << " ?== " << IB(arcA, br_prev_mat_idx_pos ) + jumpGapCostB << endl;

		if (IBDmat(arcA.idx(), arcB.idx()) ==  IB(arcA, br_prev_mat_idx_pos ) + jumpGapCostB )
		    {
			trace_IX(bl, br_prev_mat_idx_pos, arcA, false, sv);
			for ( size_type k = br_prev_seq_pos + 1; k < br_seq_pos; k++)
			    {
				alignment.append(-1, k);
			    }
			return;
		    }
	    }
	if (trace_debugging_output) std::cout << "WARNING: trace_IXD No trace found!" << std::endl;

	return;
    }


    // AlignerN: traceback
    template<class ScoringView>
    void AlignerN::trace_D(const Arc &arcA, const Arc &arcB, ScoringView sv) {

	if (trace_debugging_output) std::cout << "****trace_D****" << arcA << " " << arcB <<std::endl;
	assert(! params->struct_local_);

	seq_pos_t al = arcA.left();
	seq_pos_t ar_seq_pos = arcA.right();
	seq_pos_t bl = arcB.left();
	seq_pos_t br_seq_pos = arcB.right();
	seq_pos_t ar_prev_mat_idx_pos = mapperA.first_valid_mat_pos_before(al, ar_seq_pos);
	seq_pos_t ar_prev_seq_pos = mapperA.get_pos_in_seq_new(al, ar_prev_mat_idx_pos);
	infty_score_t jumpGapCostA = getGapCostBetween(ar_prev_seq_pos, ar_seq_pos, true);

	matidx_t br_prev_mat_idx_pos = mapperB.first_valid_mat_pos_before(bl, br_seq_pos); //tocheck: ar or ar-1?
	seq_pos_t br_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, br_prev_mat_idx_pos);
	infty_score_t jumpGapCostB = getGapCostBetween(br_prev_seq_pos, br_seq_pos, false);

	// --------------------
	// case of stacking: not supported
	assert(! scoring->stacking());

	// --------------------
	// now handle the case that arc match is not stacked

	//first compute IA
	fill_IA_entries(al, arcB, ar_seq_pos);
	if ( sv.D(arcA, arcB) == IA( ar_prev_mat_idx_pos, arcB ) + jumpGapCostA )
	    {
		assert(IADmat(arcA.idx(),arcB.idx()) == infty_score_t::neg_infty || IADmat(arcA.idx(),arcB.idx()) == sv.D(arcA, arcB));
		IADmat(arcA.idx(),arcB.idx()) = sv.D(arcA, arcB); //tocheck: prevent recomputation

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
		assert(IBDmat(arcA.idx(),arcB.idx()) == infty_score_t::neg_infty || IBDmat(arcA.idx(), arcB.idx()) == sv.D(arcA, arcB));
		IBDmat(arcA.idx(),arcB.idx()) = sv.D(arcA, arcB); //tocheck: prevent recomputation

		trace_IX(bl, br_prev_mat_idx_pos, arcA, false, sv);
		for ( size_type k = br_prev_seq_pos + 1; k < br_seq_pos; k++)
		    {
			alignment.append(-1, k);
		    }
		return;
	    }

	// first recompute M
	fill_M_entries(al,ar_seq_pos, bl, br_seq_pos,params->struct_local_);


	//-----three cases for gap extension/initiation ---

	infty_score_t opening_cost_A;
	if (ar_prev_seq_pos < (ar_seq_pos - 1)) //implicit base deletion because of sparsification
	    opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
	else
	    opening_cost_A = (infty_score_t)0;

	infty_score_t opening_cost_B;
	if (br_prev_seq_pos < (br_seq_pos - 1)) //implicit base insertion because of sparsification
	    opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
	else
	    opening_cost_B  = (infty_score_t)0;

	infty_score_t gap_score = jumpGapCostA + jumpGapCostB;

	if (sv.D(arcA, arcB) == (infty_score_t)(gap_score + opening_cost_B + Emat(ar_prev_mat_idx_pos, br_prev_mat_idx_pos)))
	    {
		trace_E(0, al, ar_prev_mat_idx_pos, bl, br_prev_mat_idx_pos, false, def_scoring_view);
	    }
	else if (sv.D(arcA, arcB) == (infty_score_t)(gap_score + opening_cost_A + Fmat(ar_prev_mat_idx_pos, br_prev_mat_idx_pos)))
	    {
		trace_F(0, al, ar_prev_mat_idx_pos, bl, br_prev_mat_idx_pos, false, def_scoring_view);
	    }
	else if (sv.D(arcA, arcB) == (infty_score_t)(gap_score + opening_cost_A + opening_cost_B +  Ms[0](ar_prev_mat_idx_pos, br_prev_mat_idx_pos) ))
	    {
		trace_M(0, al, ar_prev_mat_idx_pos, bl, br_prev_mat_idx_pos, false, def_scoring_view);
	    }
	else //todo: throw exception?
	    std::cerr << "No Trace was found! ****trace_D****" << arcA << " " << arcB  << std::endl;
	//------------------------------


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

    // do the trace for base deletion within one arc match.
    // only the cases without exclusions are considered
    template <class ScoringView>
    void AlignerN::trace_E(int state,pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool top_level, ScoringView sv)
    {
	assert (state == E_NO_NO);
	seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
	if (trace_debugging_output) std::cout << "******trace_E***** " << " al:" << al << " bl:"<< bl << " i:" << i_seq_pos << " :: " <<  Emat(i_index,j_index) << std::endl;

	assert (i_seq_pos > al );

	seq_pos_t i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1); //TODO: Check border i_index==1,0


	// base del
	infty_score_t gap_cost = getGapCostBetween(i_prev_seq_pos, i_seq_pos, true) + sv.scoring()->gapA(i_seq_pos) ;
	if (Emat(i_index, j_index) == gap_cost + Emat(i_index-1,j_index) )
	    {
		if (trace_debugging_output) std::cout << "base deletion E" << i_index-1 << " , " << j_index << std::endl;
		trace_E(state, al, i_index-1, bl, j_index, top_level, sv);
		alignment.append(i_seq_pos, -1);
		return;

	    }
	else  if (Emat(i_index, j_index) == Ms[state](i_index-1, j_index) + gap_cost + sv.scoring()->indel_opening())
	    {
		if (trace_debugging_output) std::cout << "base deletion M" << i_index-1 << " , " << j_index << std::endl;
		trace_M(state, al, i_index-1, bl, j_index, top_level, sv);
		alignment.append(i_seq_pos, -1);
		return;
	    }
	if (trace_debugging_output) std::cout << "WARNING: trace_E No trace found!" << std::endl;
    }

    template <class ScoringView>
    void AlignerN::trace_F(int state,pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool top_level, ScoringView sv)
    {

	assert (state == E_NO_NO);
	seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);

	if (trace_debugging_output) std::cout << "******trace_F***** " << " al:" << al << " bl:"<< bl << " j:" << j_seq_pos << " :: " <<  Fmat(i_index,j_index) << std::endl;

	assert (j_seq_pos > bl );

	seq_pos_t j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1); //TODO: Check border j_index==1,0

	// base ins
	infty_score_t gap_cost = getGapCostBetween(j_prev_seq_pos, j_seq_pos, false) + sv.scoring()->gapB(j_seq_pos);
	if (Fmat(i_index, j_index) == Fmat(i_index,j_index-1) + gap_cost)
	    {
		if (trace_debugging_output) std::cout << "base insertion F" << i_index << " , " << j_index-1 << std::endl;
		trace_F(state, al, i_index, bl, j_index-1, top_level, sv);
		alignment.append(-1, j_seq_pos);
		return;
	    }
	else if (Fmat(i_index, j_index) == Ms[state](i_index, j_index-1) + gap_cost  + sv.scoring()->indel_opening())
	    {
		if (trace_debugging_output) std::cout << "base insertion M" << i_index << " , " << j_index-1 << std::endl;
		trace_M(state, al, i_index, bl, j_index-1, top_level, sv);
		alignment.append(-1, j_seq_pos);
		return;
	    }
	if (trace_debugging_output) std::cout << "WARNING: trace_F No trace found!" << std::endl;
    }


    // trace and handle all cases that do not involve exclusions
    template<class ScoringView>
    void AlignerN::trace_M_noex(int state, pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool top_level, ScoringView sv)
    {

	assert (state == E_NO_NO);
	M_matrix_t &M=Ms[state];


	seq_pos_t i_seq_pos = mapperA.get_pos_in_seq_new(al, i_index);
	seq_pos_t j_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index);

	assert (i_seq_pos >= al );
	assert (j_seq_pos >= bl );

	if ( i_seq_pos == al && j_seq_pos == bl )
	    return;


	seq_pos_t i_prev_seq_pos = al;//tocheck: Important
	if ( i_seq_pos > al )
	    i_prev_seq_pos = mapperA.get_pos_in_seq_new(al, i_index-1); //TODO: Check border i_index==1,0
	seq_pos_t j_prev_seq_pos = bl;
	if (j_seq_pos > bl )
	    j_prev_seq_pos = mapperB.get_pos_in_seq_new(bl, j_index-1); //TODO: Check border j_index==1,0
	bool constraints_alowed_edge = true; // constraints are not considered,  params->constraints_->allowed_edge(i_seq_pos, j_seq_pos)
	bool constraints_aligned_pos_A = false; // TOcheck: Probably unnecessary, constraints are not considered
	bool constraints_aligned_pos_B = false; // TOcheck: Probably unnecessary, constraints are not considered
	// determine where we get M(i,j) from



	// std::cout << i << " " << j << " " << sv.scoring()->basematch(i,j)<<std::endl;

	if (  i_seq_pos > al && j_seq_pos > bl &&
	      constraints_alowed_edge )
	    {
		//------------------------------------
		// calculate possible opening gap costs
		infty_score_t opening_cost_A;
		if (i_prev_seq_pos < (i_seq_pos - 1)) //implicit base deletion because of sparsification
		    opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
		else
		    opening_cost_A = (infty_score_t)0;

		infty_score_t opening_cost_B;
		if (j_prev_seq_pos < (j_seq_pos - 1)) //implicit base insertion because of sparsification
		    opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
		else
		    opening_cost_B  = (infty_score_t)0;
		//------------------------------------
		// base match

		infty_score_t gap_match_score = getGapCostBetween(i_prev_seq_pos, i_seq_pos, true) + getGapCostBetween(j_prev_seq_pos, j_seq_pos, false) + (sv.scoring()->basematch(i_seq_pos, j_seq_pos));
		//base match and continue with deletion
		if (M(i_index,j_index) == (infty_score_t)(gap_match_score + opening_cost_B + Emat(i_index-1, j_index-1)) )
		    {
			if (trace_debugging_output) std::cout << "base match E" << i_index << " , " << j_index << std::endl;
			trace_E(state, al, i_index-1, bl, j_index-1, top_level, sv );
			alignment.append(i_seq_pos,j_seq_pos);
			return;

		    }
		else   	//base match and continue with insertion
		    if (M(i_index,j_index) == (infty_score_t)(gap_match_score + opening_cost_A + Fmat(i_index-1, j_index-1)) )
			{
			    if (trace_debugging_output) std::cout << "base match F" << i_index << " , " << j_index << std::endl;
			    trace_F(state, al, i_index-1, bl, j_index-1, top_level, sv );
			    alignment.append(i_seq_pos,j_seq_pos);
			    return;
			}
		    else	//base match, then continue with M case again, so both gap opening costs(if possible) should be included
			if (M(i_index,j_index) == (infty_score_t)(gap_match_score + opening_cost_A + opening_cost_B + M(i_index-1, j_index-1)) )
			    {
				if (trace_debugging_output) std::cout << "base match M" << i_index << " , " << j_index << std::endl;
				trace_M(state, al, i_index-1, bl, j_index-1, top_level, sv);
				alignment.append(i_seq_pos,j_seq_pos);
				return;
			    }

		/*	for ( size_type k = i_prev_seq_pos + 1; k < i_seq_pos; k++)
			{
			alignment.append(k, -1);
			}
			for ( size_type k = j_prev_seq_pos + 1; k < j_seq_pos; k++)
			{
			alignment.append(-1, k);
			}
		*/
	    }

	// base deletion
	if (  i_seq_pos > al &&
	      !constraints_aligned_pos_A
	      && M(i_index,j_index) == Emat(i_index, j_index) )
	    {
		if (trace_debugging_output) std::cout << "base deletion E" << i_index << " , " << j_index << std::endl;

		trace_E(state, al, i_index, bl, j_index, top_level, sv);
		/* for ( size_type k = i_prev_seq_pos + 1; k <= i_seq_pos; k++)
		   {
		   alignment.append(k, -1);
		   }*/
		return;
	    }

	// base insertion
	if (  j_seq_pos > bl &&
	      !constraints_aligned_pos_B
	      && M(i_index,j_index) == Fmat(i_index, j_index) )
	    {
		if (trace_debugging_output) std::cout << "base insertion F" << i_index << " , " << j_index << std::endl;

		trace_F(state, al, i_index, bl, j_index, top_level, sv);
		/*for ( size_type k = j_prev_seq_pos + 1; k <= j_seq_pos; k++)
		  {
		  alignment.append(-1, k);
		  }*/
		return;
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
		infty_score_t opening_cost_A;
		if (arcA_left_seq_pos_before < (arcA.left() - 1)) //implicit base deletion because of sparsification
		    opening_cost_A = (infty_score_t)(sv.scoring()->indel_opening());
		else
		    opening_cost_A = (infty_score_t)0;

		for (ArcIdxVec::const_iterator arcBIdx = arcsB.begin(); arcBIdx != arcsB.end(); ++arcBIdx)
		    {
			const Arc& arcB = bpsB.arc(*arcBIdx);
			matidx_t arcB_left_index_before = mapperB.first_valid_mat_pos_before(bl, arcB.left());
			seq_pos_t arcB_left_seq_pos_before = mapperB.get_pos_in_seq_new(bl, arcB_left_index_before);

			infty_score_t opening_cost_B;
			if (arcB_left_seq_pos_before < (arcB.left() - 1)) //implicit base insertion because of sparsification
			    opening_cost_B = (infty_score_t)(sv.scoring()->indel_opening());
			else
			    opening_cost_B = (infty_score_t)0;

			infty_score_t gap_match_score = getGapCostBetween( arcA_left_seq_pos_before, arcA.left(), true) + getGapCostBetween( arcB_left_seq_pos_before, arcB.left(), false)
			    + sv.D( arcA, arcB ) + sv.scoring()->arcmatch(arcA, arcB);


			//arc match, then continue with deletion
			if ( M(i_index, j_index) ==	(infty_score_t)(gap_match_score + opening_cost_B + Emat (arcA_left_index_before, arcB_left_index_before)) )
			    {

				if (trace_debugging_output) std::cout << "arcmatch E"<< arcA <<";"<< arcB << " :: "   << std::endl;

				trace_E(state, al, arcA_left_index_before, bl, arcB_left_index_before, top_level, sv);


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

				assert(! params->no_lonely_pairs_);
				trace_D(arcA, arcB, sv);
				alignment.append(arcA.right(),arcB.right());
				return;

			    }
			//arc match, then continue with insertion case
			else if ( M(i_index, j_index) ==
				  (infty_score_t)(gap_match_score + opening_cost_A + Fmat (arcA_left_index_before, arcB_left_index_before)) )
			    {

				if (trace_debugging_output) std::cout << "arcmatch F"<< arcA <<";"<< arcB << " :: "   << std::endl;

				trace_F(state, al, arcA_left_index_before, bl, arcB_left_index_before, top_level, sv);


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

				assert(! params->no_lonely_pairs_);
				trace_D(arcA, arcB, sv);
				alignment.append(arcA.right(),arcB.right());
				return;

			    }
			//arc match, then continue with general M case
			else if ( M(i_index, j_index) == gap_match_score  + opening_cost_A + opening_cost_B + M(arcA_left_index_before, arcB_left_index_before) )
			    {

				if (trace_debugging_output) std::cout << "arcmatch M"<< arcA <<";"<< arcB << " :: "   << std::endl;

				trace_M(state, al, arcA_left_index_before, bl, arcB_left_index_before, top_level, sv);

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

				assert(! params->no_lonely_pairs_);
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
    AlignerN::trace_M(int state,pos_type al, matidx_t i_index, pos_type bl, matidx_t j_index, bool top_level, ScoringView sv) {
	//pre: M matrices for arc computed
	M_matrix_t &M=Ms[state];
	assert (state == E_NO_NO);
	assert(! params->sequ_local_); //Local seq alignment not implemented yet.

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
	    trace_M_noex(state, al, i_index, bl, j_index, top_level, sv);
	    break;
	}
    }

    template<class ScoringView>
    void
    AlignerN::trace(ScoringView sv) {
	// pre: last call align_in_arcmatch(r.get_startA()-1,r.get_endA()+1,r.get_startB()-1,r.get_endB()+1);
	//      or align_top_level_locally for sequ_local alignent

	// reset the alignment strings (to empty strings)
	// such that they can be written again during the trace
	alignment.clear();

	// free end gap version: trace_M(E_NO_NO,r.get_startA()-1,max_i,r.get_startB()-1,max_j,true,sv);
	seq_pos_t ps_al = r.get_startA() - 1;
	matidx_t last_mat_idx_pos_A = mapperA.number_of_valid_mat_pos(ps_al) -1;//tocheck: check the correctness
	//seq_pos_t last_seq_pos_A = mapperA.get_pos_in_seq_new(ps_al, last_mat_idx_pos_A);

	seq_pos_t ps_bl = r.get_startB() - 1;
	matidx_t last_mat_idx_pos_B = mapperB.number_of_valid_mat_pos(ps_bl) -1;//tocheck: check the correctness
	//seq_pos_t last_seq_pos_B = mapperB.get_pos_in_seq_new(ps_bl, last_mat_idx_pos_B);

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
