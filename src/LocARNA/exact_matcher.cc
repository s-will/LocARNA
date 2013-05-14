#include "sequence.hh"
#include "arc_matches.hh"
#include "exact_matcher.hh"
#include <iostream>
#include <fstream>

using namespace std;

namespace LocARNA {

//todo: remove
bool debug = false;
bool debug_trace_G = false;
bool debug_VG = false;
bool debug_trace_LGLR_subopt = false;
bool debug_trace_F_heuristic = false;

	// Constructor
    ExactMatcher::ExactMatcher(const Sequence &seqA_,
			       const Sequence &seqB_,
			       const ArcMatches &arc_matches_,
			       const SparseTraceController &sparse_trace_controller_,
			       PatternPairMap &foundEPMs_,
			       const int &alpha_1_,
			       const int &alpha_2_,
			       const int &alpha_3_,
			       const int &difference_to_opt_score_,
			       const int &min_subopt_score_,
			       const int &easier_scoring_par_,
			       const double &subopt_range_,
			       const int &am_threshold_,
			       const double &cutoff_coverage_
			       )
	: seqA(seqA_),
	  seqB(seqB_),
	  arc_matches(arc_matches_),
	  bpsA(arc_matches_.get_base_pairsA()),
	  bpsB(arc_matches_.get_base_pairsB()),
	  sparse_trace_controller(sparse_trace_controller_),
	  sparse_mapperA(sparse_trace_controller.get_sparse_mapperA()),
	  sparse_mapperB(sparse_trace_controller.get_sparse_mapperB()),
	  foundEPMs(foundEPMs_),
	  alpha_1(alpha_1_),
	  alpha_2(alpha_2_),
	  alpha_3(alpha_3_),
	  difference_to_opt_score(difference_to_opt_score_*100),
	  min_subopt_score(min_subopt_score_*100),
	  easier_scoring_par(easier_scoring_par_),
	  subopt_range(subopt_range_),
	  am_threshold(am_threshold_*100),
	  cutoff_coverage(cutoff_coverage_),
	  pseudo_arcA(bpsA.num_bps(),0,seqA.length()),
	  pseudo_arcB(bpsB.num_bps(),0,seqB.length())
    {

    	// set size of matrices
    	if(debug) cout << "max dimensions " << sparse_mapperA.get_max_info_vec_size()
    			       <<  "x" << sparse_mapperB.get_max_info_vec_size() << endl;

    	L.resize(sparse_mapperA.get_max_info_vec_size(),sparse_mapperB.get_max_info_vec_size());
    	L.fill(infty_score_t::neg_infty);
    	L.set(0,0,infty_score_t(0));

    	G_A.resize(sparse_mapperA.get_max_info_vec_size(),sparse_mapperB.get_max_info_vec_size());
    	G_AB.resize(sparse_mapperA.get_max_info_vec_size(),sparse_mapperB.get_max_info_vec_size());

    	LR.resize(sparse_mapperA.get_max_info_vec_size(),sparse_mapperB.get_max_info_vec_size());
    	LR.fill(infty_score_t::neg_infty);
    	LR.set(0,0,infty_score_t(0));

    	F.resize(seqA.length()+1,seqB.length()+1);
    	F.fill(infty_score_t(0));

    	Dmat.resize(bpsA.num_bps(),bpsB.num_bps());
    	Dmat.fill(infty_score_t::neg_infty); //initialize all arcmatches with -inf

    }

    // Destructor
    ExactMatcher::~ExactMatcher(){}

    // initialization of the gap matrices for the suboptimal traceback
    void ExactMatcher::initialize_gap_matrices(){

    	// initialize first row of G_A with -inf
    	for(pos_type j=0;j<G_A.sizes().second;++j){
    		G_A.set(0,j,infty_score_t::neg_infty);
    	}

    	// initilize first column of G_AB with -inf
    	for(pos_type i=0;i<G_AB.sizes().first;++i){
    		G_AB.set(i,0,infty_score_t::neg_infty);
    	}

    }


    // ---------------------------------------------------------------------------------------------------------
    // fill matrices

    //initializes the F matrix for using the trace controller
    void ExactMatcher::init_Fmat(){
    	if(debug) cout << "init F mat " << endl;


    	// initialize for whole F matrix
    	pos_type al=0;
    	pos_type ar=seqA.length();
    	pos_type bl=0;
    	pos_type br=seqB.length();

    	// al,bl can only be reached in states, where this is legal with cost 0 for empty alignment
    	F(al,bl) = (infty_score_t)0;

    	pos_type i;
    	for (i=al+1; i<ar; ++i) {
    		if (sparse_trace_controller.min_col(i)>bl) break; // fill only as long as column bl is accessible
    		if(debug) cout << "F( " << i << "," << bl << ")=" << 0 << endl;
    		F(i,bl) = (infty_score_t)0;
    	}

    	// fill entries left of valid entries
    	for ( ; i<ar; ++i) {
    		assert(sparse_trace_controller.min_col(i)>bl);
    		if(debug) cout << "F(" << i << "," << sparse_trace_controller.min_col(i)-1 << ")=-inf" << endl;
    		F(i,sparse_trace_controller.min_col(i)-1) = infty_score_t::neg_infty;
    	}

    	// init first row al

    	pos_type j;
    	for (j=bl+1 ; j < std::min(br, sparse_trace_controller.max_col(al)+1) ; j++) {
    		if(debug) cout << "F( " << al << "," << j << ")=" << 0 << endl;
    		F(al,j) = (infty_score_t)0;
    	}

    	// fill entries above valid entries
    	// here j points to one position right of the last initialized entry in row al
    	for (i=al+1; i<ar; ++i) {
    		for (;
    				j<std::min(br,sparse_trace_controller.max_col(i)+1); ++j) {
    			if(debug) cout << "F(" << i-1 << "," << j << ")=-inf" << endl;
    			F(i-1,j)=infty_score_t::neg_infty;
    		}
    	}
    	if(debug) cout << "F" << F << endl;
    }

    // initialize all other compressed matrices (L, G_A, G_AB and LR)
    // for using the trace controller
    void ExactMatcher::init_mat(ScoreMatrix &mat, const Arc &a, const Arc &b,
    		infty_score_t first_entry, infty_score_t first_col, infty_score_t first_row) {

    	size_type num_posA = sparse_mapperA.number_of_valid_mat_pos(a.idx());
    	size_type num_posB = sparse_mapperB.number_of_valid_mat_pos(b.idx());

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	if(debug){
    		for (matidx_t idx_i=0; idx_i<num_posA; ++idx_i) {
    			//size_t min_idx_j = sparse_trace_controller.min_col_idx(a,b,idx_i);
    			matidx_t min_idx_j = sparse_trace_controller.min_col_idx(idxA,idxB,idx_i,b.left());
    			//size_t idx_after_max_idx_j = sparse_trace_controller.idx_after_max_col_idx(a,b,idx_i);
    			matidx_t idx_after_max_idx_j = sparse_trace_controller.idx_after_max_col_idx(idxA,idxB,idx_i,b.left());
    			cout << idx_i<< ":" << min_idx_j << "," <<idx_after_max_idx_j << endl;
    		}
    	}

    	// initialize matrix entry (0,0) with the corresponding score
    	mat(0,0) = (infty_score_t)first_entry;

    	//init first col

    	matidx_t idx_i;
    	matidx_t min_idx_j;
    	for (idx_i=1; idx_i<num_posA; ++idx_i) {

    		min_idx_j = sparse_trace_controller.min_col_idx(idxA,idxB,idx_i,b.left());
    		if (min_idx_j>0) break; // fill only as long as column bl is accessible
    		mat(idx_i,0) = (infty_score_t)first_col;
    	}

    	// fill entries left of valid entries
    	for ( ; idx_i<num_posA; ++idx_i) {
    		matidx_t min_idx_j = sparse_trace_controller.min_col_idx(idxA,idxB,idx_i,b.left());
    		assert(min_idx_j>0);
    		mat(idx_i,min_idx_j-1) = infty_score_t::neg_infty;
    	}

    	// init first row al
    	pos_type idx_j;
    	matidx_t idx_after_max_idx_j = sparse_trace_controller.idx_after_max_col_idx(idxA,idxB,0,b.left());
    	for (idx_j=1 ; idx_j < std::min(num_posB, idx_after_max_idx_j) ; ++idx_j) {
    		mat(0,idx_j) = (infty_score_t)first_row;
    	}

    	// fill entries above valid entries
    	// here j points to one position right of the last initialized entry in row al
    	for (idx_i=1; idx_i<num_posA; ++idx_i) {
    		idx_after_max_idx_j = sparse_trace_controller.idx_after_max_col_idx(idxA,idxB,idx_i,b.left());
    		for (;
    				idx_j<std::min(num_posB,idx_after_max_idx_j); ++idx_j) {
    			mat(idx_i-1,idx_j)=infty_score_t::neg_infty;
    		}
    	}

    }

    // compute arcmatch score by filling matrices L, G and LR (method compute_LGLR)
    // and computing the arcmatch score and store it in matrix D
    // store arcmatch_score with stacking and probs of outermost arcmatch
    void ExactMatcher::compute_arcmatch_score() {

    	matpos_t last_filled_pos; //the last position that was filled in the matrices

    	// for all arc matches from inside to outside
    	for(ArcMatchVec::const_iterator it=arc_matches.begin();it!=arc_matches.end();++it){

    		pos_type al=it->arcA().left();
    		pos_type ar=it->arcA().right();
    		pos_type bl=it->arcB().left();
    		pos_type br=it->arcB().right();

    		// compute the arc match score only for matching arc matches
    		if(seqA[al]==seqB[bl] && seqA[ar]==seqB[br]){

    			last_filled_pos=compute_LGLR(it->arcA(),it->arcB(),false);

    			matidx_t last_i = last_filled_pos.first;
    			matidx_t last_j = last_filled_pos.second;

    			// the arc match score is the maximum of the last matrix entry in
    			// matrices LR, L or G_A (as we used the heuristic computation)
    			D(*it) =max3(LR(last_i,last_j),L(last_i,last_j),G_A(last_i,last_j));
    		}
    	}

    	// compute the best combination of arc matches and unpaired parts in matrix F
    	compute_F();
    }


    // computes for a given pair of arcs the matrices L, G and LR
    // store -inf in last cell of matrix LR if a gap between the
    // last matched positions and the right ends of the arcs exists
    // compute L, G_A (G matrix) and LR matrix
    ExactMatcher::matpos_t ExactMatcher::compute_LGLR(const Arc &a, const Arc &b, bool suboptimal){

    	// initialize matrices for using the sparse trace controller in the heuristic case
    	if(!suboptimal){
    		init_mat(L,a,b,infty_score_t(0),infty_score_t::neg_infty,infty_score_t::neg_infty);
    		init_mat(G_A,a,b,infty_score_t(0),infty_score_t(0),infty_score_t(0));
    		init_mat(LR,a,b,infty_score_t(0),infty_score_t::neg_infty,infty_score_t::neg_infty);
    	}

       	index_t idxA = a.idx();
       	index_t idxB = b.idx();

       	size_type num_posA = sparse_mapperA.number_of_valid_mat_pos(idxA);
       	size_type num_posB = sparse_mapperB.number_of_valid_mat_pos(idxB);

       	matidx_t max_j = sparse_trace_controller.idx_after_max_col_idx(idxA,idxB,0,b.left());
       	assert(max_j>0);
       	pair_seqpos_t last_pos_filled = pair_seqpos_t(0,max_j-1);
       	matidx_t idx_i,idx_j=1;

       	for(idx_i=1;idx_i<num_posA;++idx_i){

       		// determine the correct interval of valid positions in the matrices
       		matidx_t min_idx_j = suboptimal?
       				1 : sparse_trace_controller.min_col_idx(idxA,idxB,idx_i,b.left());
       		matidx_t idx_after_max_idx_j = suboptimal?
       				num_posB : sparse_trace_controller.idx_after_max_col_idx(idxA,idxB,idx_i,b.left());

       		for(idx_j=std::max((matidx_t)1,min_idx_j);idx_j<idx_after_max_idx_j;++idx_j){

       			matpos_t mat_pos = matpos_t(idx_i,idx_j);
       			pair_seqpos_t seq_pos = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,mat_pos);

       			// determine the next position to the left on the diagonal
       			matpos_t idx_pos_diag = suboptimal ?
       					matpos_t(idx_i-1,idx_j-1) :
       					sparse_trace_controller.first_valid_mat_pos_before_with_tc(idxA,idxB,seq_pos,a.left(),b.left());

       			L(idx_i,idx_j)=compute_matrix_entry(a,b,mat_pos,idx_pos_diag,false,suboptimal);

       			if(suboptimal){
       				G_A(idx_i,idx_j)=max(L(idx_pos_diag.first,idx_j),
       									 G_A(idx_pos_diag.first,idx_j));

       				G_AB(idx_i,idx_j)=max3(L(idx_i,idx_pos_diag.second),
       									   G_A(idx_i,idx_pos_diag.second),
       									   G_AB(idx_i,idx_pos_diag.second));
       			}
       			else{
       				G_A(idx_i,idx_j)=max(L(idx_i,idx_j),
       									 G_A(idx_pos_diag.first,idx_pos_diag.second));

       				// check whether border cases access valid positions
       				if(sparse_trace_controller.is_valid_idx_pos(idxA,idxB,matpos_t(idx_pos_diag.first,idx_j))){
       					G_A(idx_i,idx_j) = max(G_A(idx_i,idx_j),
       							               G_A(idx_pos_diag.first,idx_j));
       				}

       				if(sparse_trace_controller.is_valid_idx_pos(idxA,idxB,matpos_t(idx_i,idx_pos_diag.second))){
       					G_A(idx_i,idx_j) = max(G_A(idx_i,idx_j),
       							               G_A(idx_i,idx_pos_diag.second));
       				}
       			}

       			LR(idx_i,idx_j)=compute_matrix_entry(a,b,mat_pos,idx_pos_diag,true,suboptimal);

       			// update last filled position
       			last_pos_filled.first=idx_i;
       			last_pos_filled.second=idx_j;
       		}
       	}

       	if(sparse_trace_controller.get_delta() == (size_type)-1){
       		if(num_posA>1 && num_posB>1) assert(last_pos_filled==matpos_t(num_posA-1,num_posB-1));
       	}

       	// if the right ends of the arcs cannot be reached without creating a gap, we store in the last
       	// filled position of matrix LR -inf
       	if((!sparse_trace_controller.matching_without_gap_possible_with_tc(idxA,idxB,last_pos_filled,
       																	   pair_seqpos_t(a.right(),b.right())))
       		 && last_pos_filled.first>0 && last_pos_filled.second>0){
       			LR(last_pos_filled.first,last_pos_filled.second)=infty_score_t::neg_infty;
       	}
       	return last_pos_filled;
     }

    // computes a matrix entry (i,j) in matrix L or LR for arcs a and b
    // already taking into account the trace controller (not yet implemented for
    // the suboptimal case!)
    infty_score_t ExactMatcher::compute_matrix_entry(const Arc &a, const Arc &b,
    		matpos_t mat_pos, matpos_t mat_pos_diag, bool matrixLR, bool suboptimal){

    	if(debug){cout << "G_A ";if(suboptimal) cout << " and G_B ";cout << " are computed" << endl;}

    	infty_score_t score_seq = infty_score_t::neg_infty; //score from sequential case
    	infty_score_t score_str =infty_score_t::neg_infty; //score from structural case

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	pair_seqpos_t seq_pos = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,mat_pos);

    	seqpos_t i = seq_pos.first;
    	seqpos_t j = seq_pos.second;

    	matidx_t idx_i = mat_pos.first;
    	matidx_t idx_j = mat_pos.second;

    	if(seqA[i]==seqB[j]){
    		// sequential matching
    		// if the position is likely to be unpaired we trace the sequential match
    		if(sparse_trace_controller.pos_unpaired(idxA,idxB,matpos_t(idx_i,idx_j))){

    			/*if(matching_without_gap_possible(a,b,pos(i,j))){
    				score_seq = mat(idx_i_diag,idx_j_diag)+score_for_seq_match();
    			}
    			//if the previous entry in matrix B was taken from matrix G
    			//-> no check if sequential matching is permitted
    			if(matrixLR){
    				score_seq=max(L(idx_i_diag,idx_j_diag)+score_for_seq_match(),score_seq);
    				score_seq = max(G_A(idx_i_diag,idx_j_diag)+score_for_seq_match(),score_seq);
    				if(suboptimal) score_seq = max(G_AB(idx_i_diag,idx_j_diag)+score_for_seq_match(),score_seq);
    			}*/

    			score_seq = seq_str_matching(a,b,mat_pos_diag,seq_pos,
    					score_for_seq_match(),matrixLR,suboptimal);
    		}
    		//structural matching
    		for(ArcIdxVec::const_iterator itA=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).begin();
    				itA!=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).end();++itA){
    			for(ArcIdxVec::const_iterator itB=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).begin();
    					itB!=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).end();++itB){

    				const Arc &inner_a = bpsA.arc(*itA);
    				const Arc &inner_b = bpsB.arc(*itB);

    				const infty_score_t &score_for_inner_am = score_for_am(inner_a,inner_b);
    				if(score_for_inner_am.is_neg_infty()) continue;

    				pair_seqpos_t last_seq_pos_to_be_matched(inner_a.left(),inner_b.left());

    				// determine the first matrix position before the left ends of the arcs
    				matpos_t mat_pos_diag_str =
    						sparse_trace_controller.first_valid_mat_pos_before_with_tc(idxA,idxB,
    								last_seq_pos_to_be_matched,a.left(),b.left());

    				assert(score_for_inner_am.is_finite());

    				score_t score_am_stacking=score_for_inner_am.finite_value()+score_for_stacking(a,b,inner_a,inner_b);

    				/*matching in A and B, respectively without gap possible
    				/if(matching_without_gap_possible(a,b,pos(inner_a.left(),inner_b.left()))){
    					score_str=max(score_str,mat(idx_i_before,idx_j_before)
    							+score_am_stacking);
    				}
    				//-> in matrix LR always allowed
    				if(matrixLR){
    					score_str=max(score_str,L(idx_i_before,idx_j_before)
    					    	+score_am_stacking);
    					score_str=max(score_str,G_A(idx_i_before,idx_j_before)
    							+score_am_stacking);
    					if(suboptimal) score_str=max(score_str,G_AB(idx_i_before,idx_j_before)
    											+score_am_stacking);
    				}*/

    				score_str = max(seq_str_matching(a,b,mat_pos_diag_str,
    						last_seq_pos_to_be_matched,score_am_stacking,matrixLR,suboptimal),score_str);
    			}
    		}
    	}
    	return max(score_seq,score_str);
    }

    // computes the score for a sequential or structural matching by checking the three possibilities
    // (either continue traceback in matrix mat, and if we are currently in matrix LR, we can continue
    // the traceback in L or G_A (the gap matrix)
    infty_score_t ExactMatcher::seq_str_matching(const Arc &a, const Arc &b, matpos_t mat_pos_diag,
    		pair_seqpos_t seq_pos_to_be_matched, score_t add_score, bool matrixLR, bool suboptimal){

    	infty_score_t score = infty_score_t::neg_infty;

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	matidx_t idx_i_diag = mat_pos_diag.first;
    	matidx_t idx_j_diag = mat_pos_diag.second;

    	ScoreMatrix &mat= matrixLR ? LR : L;

    	// if matching without a gap is possible we simply add add_score
    	if(sparse_trace_controller.matching_without_gap_possible_with_tc(idxA,idxB,mat_pos_diag,seq_pos_to_be_matched)){
    		score = mat(idx_i_diag,idx_j_diag)+add_score;
    	}

    	// if an entry in LR is computed, we can also come from L, G_A in the heuristic case and
    	// L, G_A and G_AB in the suboptimal case
    	if(matrixLR){
    		score=max(L(idx_i_diag,idx_j_diag)+add_score,score);
    		score = max(G_A(idx_i_diag,idx_j_diag)+add_score,score);
    		if(suboptimal) score = max(G_AB(idx_i_diag,idx_j_diag)+add_score,score);
    	}
    	return score;
    }

    // computes the final matrix F from which the final EPMs can be traced
    void ExactMatcher::compute_F(){

       	init_Fmat();
       	score_t max_in_F=0; // the maximal score of an EPM
       	infty_score_t score_seq,score_str; // sequential and structural score

       	for(seqpos_t i=1;i<F.sizes().first;++i){
       		// determine the correct interval for row i
       		for(seqpos_t j=std::max(seqpos_t(1),sparse_trace_controller.min_col(i));
       				j<=sparse_trace_controller.max_col(i);++j){

       			score_seq=infty_score_t(0);score_str=infty_score_t(0);

       			//sequential matching
       			if(seqA[i]==seqB[j]){

       				score_seq = F(i-1,j-1)+score_for_seq_match();

       			}

       			//structural matching
       			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
       					arc_matches.common_right_end_list(i,j).end() != it; ++it ) {

       				const ArcMatch &am = arc_matches.arcmatch(*it);
       				const Arc &a = am.arcA();
       				const Arc &b = am.arcB();

       				// if the score for the arc match is less than the threshold,
       				// we go to the next arc match
       				if(score_for_am(a,b)<FiniteInt(am_threshold)) continue;

       				assert(sparse_trace_controller.is_valid(a.left()-1,b.left()-1));

       				score_str = max(score_str,score_for_am(a,b)+
       							F(a.left()-1,b.left()-1));
       			}

       			F(i,j)=max(score_seq,score_str);

       			if(F(i,j)>(infty_score_t)max_in_F){
       				assert(F(i,j).is_finite());
       				max_in_F=F(i,j).finite_value();
       				pos_of_max = pair_seqpos_t(i,j);
       			}
       		}
       	}
       	cout << "max in F " << max_in_F << endl;
       	cout << "pos of max " << pos_of_max << endl;
    }

    // determines the EPMs via traceback
    void ExactMatcher::trace_EPMs(bool suboptimal){
    	//trace EPMs
    	if(suboptimal) trace_EPMs_suboptimal();
    		else trace_EPMs_heuristic();
    }


    // ---------------------------------------------------------------------------------------------------------
    // helper functions

    // returns the score for a sequential match
    score_t ExactMatcher::score_for_seq_match(){
    	return alpha_1*100;
    }

    // returns the score for an arcmatch (the basepair match itself plus the part under it)
    // if easier_scoring_par is set, each basepair match is scored equally
    // otherwise the base pair probabilities are taken into account
    infty_score_t ExactMatcher::score_for_am(const Arc &a, const Arc &b){

    	infty_score_t result=infty_score_t(0);

    	if(easier_scoring_par){
    		result= D(a,b)+FiniteInt(((2*alpha_1)+easier_scoring_par*2)*100);
    	}

    	else{
    		double probArcA = bpsA.get_arc_prob(a.left(),a.right());
    		double probArcB = bpsB.get_arc_prob(b.left(),b.right());

    		result= D(a,b) + FiniteInt(((2*alpha_1)+(probArcA+probArcB)*alpha_2)*100);
    	}
    	return result;
    }


    // todo: use get_arc_stack_prob (conditional probability) instead?
    // returns the stacking score
    // checks whether a and inner_a (and b and inner_b) are stacked and adds the
    // joint (conditional better?) probability that a and inner_a (b and inner_b) occur simultaneously
    score_t ExactMatcher::score_for_stacking(const Arc &a, const Arc &b,
    		const Arc &inner_a,const Arc &inner_b){

    	double prob_stacking_arcA = 0;
    	double prob_stacking_arcB = 0;

    	//stacking arcA
    	if(a.left()+1==inner_a.left() &&
    			a.right()==inner_a.right()+1){
    		prob_stacking_arcA = bpsA.get_arc_2_prob(a.left(),a.right());
    	}

    	//stacking arcB
    	if(b.left()+1==inner_b.left() &&
    			b.right()==inner_b.right()+1){
    		prob_stacking_arcB = bpsB.get_arc_2_prob(b.left(),b.right());
    	}

    	return (prob_stacking_arcA+prob_stacking_arcB)*100*alpha_3;
    }

    //adds a found EPM to the patternPairMap (datastructure used for chaining algorithm)
    void ExactMatcher::add_foundEPM(EPM &cur_epm){

    	static int count = 0;

    	static string seq1_id = seqA.seqentry(0).name();
    	static string seq2_id = seqB.seqentry(0).name();

    	++count;
    	//todo: remove again
    	if(count%10000==0) cout << "count " << count << endl;

    	// sort the pattern vector of the current epm according
    	// to increasing positions
    	cur_epm.sort_patVec();

    	// make sure that the current epm is valid
    	assert(validate_epm(cur_epm));

    	stringstream ss;
    	ss << "pat_" << count;
    	string patId= ss.str();

    	// rewrite information for use in the chaining algorithm
    	intVec pat1Vec;
    	intVec pat2Vec;
    	string structure;

    	for(EPM::pat_vec_t::const_iterator it= cur_epm.begin();it!=cur_epm.end();++it){
    		pat1Vec.push_back(it->first);
    		pat2Vec.push_back(it->second);
    		structure.push_back(it->third);
    	}

    	//SinglePattern pattern1 = SinglePattern(patId,seq1_id,cur_epm.getPat1Vec());
    	//SinglePattern pattern2 = SinglePattern(patId,seq2_id,cur_epm.getPat2Vec());
    	//foundEPMs.add(patId, pattern1, pattern2, cur_epm.getStructure(), cur_epm.get_score() );

    	//todo: check score in chaining -> small example heuristic vs. suboptimal
    	SinglePattern pattern1 = SinglePattern(patId,seq1_id,pat1Vec);
    	SinglePattern pattern2 = SinglePattern(patId,seq2_id,pat2Vec);
    	foundEPMs.add(patId, pattern1, pattern2, structure, cur_epm.get_score() );

    }

    // ---------------------------------------------------------------------------------------------------------
    // heuristic (with TraceController)

    // trace through F matrix and L, G_A and LR matrices to find maximally extended EPMs
    // from each position (i,j)
    void ExactMatcher::trace_EPMs_heuristic(){

    	cout << "compute EPMs heuristic with trace controller " << endl;
    	cout << "min score " << min_subopt_score << endl;

    	for(size_type i=1;i<F.sizes().first;++i){
    		if(debug) cout << "for i range " << sparse_trace_controller.min_col(i) << ","
    				       << sparse_trace_controller.max_col(i) << endl;

    		pos_type min_col = std::max((size_type)1,sparse_trace_controller.min_col(i));
    		pos_type max_col = std::min(F.sizes().second,sparse_trace_controller.max_col(i)+1);

    		for(size_type j=min_col;j<max_col;++j){
    			if(F(i,j)>(infty_score_t)min_subopt_score){

    				// check whether the position needs to be traced (if the corresponding
    				// EPM is maximally extended
    				if(i==F.sizes().first-1 || j==F.sizes().second-1 || seqA[i+1]!=seqB[j+1]){

    					EPM cur_epm;
    					// compute traceback from position (i,j)
    					trace_F_heuristic(i,j,cur_epm);
    					// store the traced epm in the corresponding datastructure
    					add_foundEPM(cur_epm);

    				}
    			}
    		}
    	}
    }

    // traces through the F matrix from position (i,j) to find the best EPM that ends in (i,j)
    void ExactMatcher::trace_F_heuristic(pos_type i, pos_type j, EPM &cur_epm){

    	assert(F(i,j).is_finite());
    	cur_epm.set_score(F(i,j).finite_value());

    	while(!(F(i,j)==(infty_score_t)0)){

    		if(debug_trace_F_heuristic) cout << "i,j " << i << "," << j << endl;

    		assert(i>=1 && j>=1);
    		assert(sparse_trace_controller.is_valid(i,j));

    		if(F(i-1,j-1)+score_for_seq_match()==F(i,j)){
    			cur_epm.add(i,j,'.');
    			i--;j--;
    		}
    		else{
    			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
    					arc_matches.common_right_end_list(i,j).end() != it; ++it){

    				const ArcMatch &am = arc_matches.arcmatch(*it);

    				const Arc &a = am.arcA();
    				const Arc &b = am.arcB();

    				if(score_for_am(a,b)<FiniteInt(am_threshold)) continue;

    				if(F(a.left()-1,b.left()-1)+score_for_am(a,b)==F(i,j)){

    					cur_epm.add_am(am.arcA(),am.arcB());
    					trace_LGLR_heuristic(am.arcA(),am.arcB(),cur_epm);

    					i=am.arcA().left()-1;
    					j=am.arcB().left()-1;
    					break;
    				}
    			}
    		}
    	}
    }

    // traces through the L, G_A and LR matrices for the arcmatch of a and b
    // and stores the result in epm_to_store
    void ExactMatcher::trace_LGLR_heuristic(const Arc &a, const Arc &b, EPM &cur_epm){

    	assert(D(a,b).is_finite());

    	matpos_t cur_pos = compute_LGLR(a,b,false);

    	matidx_t idx_i=cur_pos.first;
    	matidx_t idx_j =cur_pos.second;

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	int state=in_LR;

    	//determine in which matrix to start traceback
    	if(G_A(idx_i,idx_j)==D(a,b)){
    		state=in_G_A;
    	}
    	else if(L(idx_i,idx_j)==D(a,b)){
    		state=in_L;
    	}

    	if(state==in_LR) assert(LR(idx_i,idx_j)==D(a,b));

    	//repeat until we end up in the first row or column
    	while(cur_pos.first!=0 && cur_pos.second!=0){

    		idx_i = cur_pos.first;
    		idx_j = cur_pos.second;

    		pair_seqpos_t cur_seq_pos = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,cur_pos);
    		seqpos_t i = cur_seq_pos.first;
    		seqpos_t j = cur_seq_pos.second;

    		assert(sparse_trace_controller.is_valid(i,j));

    		matpos_t mat_pos_diag = sparse_trace_controller.first_valid_mat_pos_before_with_tc(idxA,idxB,
    				                                               cur_seq_pos,a.left(),b.left());

    		matidx_t idx_i_diag = mat_pos_diag.first;
    		matidx_t idx_j_diag = mat_pos_diag.second;

    		if(sparse_trace_controller.get_delta() == (unsigned int)-1){
    			assert(mat_pos_diag == matpos_t(idx_i-1,idx_j-1));
    		}

    		switch(state){
    		case in_LR: case in_L:
    		{
    			bool seq_matching=false;
    			bool str_matching=false;

    			// check for sequential matching
    			if(sparse_trace_controller.pos_unpaired(idxA,idxB,cur_pos)){

    				/*if(matching_without_gap_possible(a,b,pos(i,j)) &&
    						mat(idx_i_diag,idx_j_diag)+score_for_seq_match()==mat(idx_i,idx_j)){
    					seq_matching=true;
    				}
    				else if(matrixLR){
    					if(L(idx_i_diag,idx_j_diag)+score_for_seq_match()==mat(idx_i,idx_j)){
    						state=in_L; seq_matching=true;
    					}
    					else if(G_A(idx_i_diag,idx_j_diag)+score_for_seq_match()==mat(idx_i,idx_j)){
    						state=in_G_A; seq_matching=true;
    					}
    				}*/

    				seq_matching = trace_seq_str_matching_heuristic(a,b,state,cur_pos,mat_pos_diag,
    						cur_seq_pos,score_for_seq_match());
    			}
    			if(seq_matching){ cur_epm.add(i,j,'.');}

    			else{
    				// check for structural matching
    				for(ArcIdxVec::const_iterator itA=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).begin();
    						itA!=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).end();++itA){
    					for(ArcIdxVec::const_iterator itB=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).begin();
    							itB!=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).end();++itB){

    						const Arc &inner_a = bpsA.arc(*itA);
    						const Arc &inner_b = bpsB.arc(*itB);

    						const infty_score_t &score_for_inner_am = score_for_am(inner_a,inner_b);
    						if(score_for_inner_am.is_neg_infty()) continue;

    						matpos_t last_seq_pos_to_be_matched(inner_a.left(),inner_b.left());

    						matpos_t idx_pos_before =
    								sparse_trace_controller.first_valid_mat_pos_before_with_tc(idxA,idxB,
    										last_seq_pos_to_be_matched,a.left(),b.left());

    						assert(score_for_inner_am.is_finite());
    						score_t score_am_stacking = score_for_inner_am.finite_value()
    																+ score_for_stacking(a,b,inner_a,inner_b);

    						/*if(matching_without_gap_possible(a,b,pos(inner_a.left(),inner_b.left())) &&
    								mat(idx_i_before,idx_j_before)+score_am_stacking==mat(idx_i,idx_j)){
    							str_matching=true;
    						}
    						else if(matrixLR){
    							if(L(idx_i_before,idx_j_before)+score_am_stacking==mat(idx_i,idx_j)){
    								state=in_L; str_matching=true;
    							}
    							else if(G_A(idx_i_before,idx_j_before)+score_am_stacking==mat(idx_i,idx_j)){
    								state=in_G_A; str_matching=true;
    							}
    						}*/

    						str_matching = trace_seq_str_matching_heuristic(a,b,state,cur_pos,idx_pos_before,
    								last_seq_pos_to_be_matched,score_am_stacking);

    						if(str_matching){
    							cur_epm.add_and_store_am(inner_a,inner_b);
    							break;
    						}
    					}
    					if(str_matching) break;
    				}
    				assert(seq_matching || str_matching);
    			}
    		}break;

    		case in_G_A:
    		{
    			assert(idx_i>=1 && idx_j>=1);

    			if(G_A(idx_i,idx_j)==L(idx_i,idx_j)){
    				state=in_L;
    			}

    			else if(G_A(idx_i,idx_j)==G_A(idx_i_diag,idx_j_diag)){
    				cur_pos=matpos_t(idx_i_diag,idx_j_diag);
    			}

    			else if(sparse_trace_controller.is_valid_idx_pos(idxA,idxB,matpos_t(idx_i_diag,idx_j))
    					&& G_A(idx_i,idx_j)==G_A(idx_i_diag,idx_j)){
    				cur_pos=matpos_t(idx_i_diag,idx_j);
    			}

    			else if(sparse_trace_controller.is_valid_idx_pos(idxA,idxB,matpos_t(idx_i,idx_j_diag))  &&
    					G_A(idx_i,idx_j)==G_A(idx_i,idx_j_diag)){
    				cur_pos=matpos_t(idx_i,idx_j_diag);
    			}

    			else{
    				cerr << "no valid traceback found " << endl;
    				return;
    			}
    		}break;
    		}
    	}

    	//make sure that traceback was successful
    	assert((state==in_G_A && G_A(cur_pos.first,cur_pos.second)==infty_score_t(0)) ||
    		   ((state==in_L || state==in_LR) && cur_pos == matpos_t(0,0)));

    	//if there are arcMatches left to process, the last arc match is processed next
    	if(cur_epm.number_of_am()>0){

    		//get index pair of the next arcs that need to be traced
    		const PairArcIdx &arc_idx_pair = cur_epm.next_arcmatch();
    		const Arc &inner_arcA = bpsA.arc(arc_idx_pair.first);
    		const Arc &inner_arcB = bpsB.arc(arc_idx_pair.second);

    		//trace recursively
    		trace_LGLR_heuristic(inner_arcA,inner_arcB,cur_epm);
    	}
    }

    // computes the trace back by checking the three possibilities (either continue
    // traceback in matrix mat, and if we are currently in matrix LR, we can continue
    // the traceback in L or G_A (the gap matrix)
    bool ExactMatcher::trace_seq_str_matching_heuristic(const Arc &a, const Arc &b,int &state,
    		matpos_t &cur_mat_pos, matpos_t mat_pos_diag,
    		pair_seqpos_t seq_pos_to_be_matched,score_t add_score){

    	bool matching = false;
    	bool matrixLR = (state==in_LR);
    	const ScoreMatrix &mat = matrixLR ? LR: L;

    	matidx_t idx_i = cur_mat_pos.first;
    	matidx_t idx_j = cur_mat_pos.second;

    	assert(sparse_trace_controller.is_valid(seq_pos_to_be_matched.first,
    											seq_pos_to_be_matched.second));

    	matidx_t idx_i_diag = mat_pos_diag.first;
    	matidx_t idx_j_diag = mat_pos_diag.second;

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	// if there is no gap on the sequence level, we check if the trace
    	// continues in the matrix mat
    	if(sparse_trace_controller.matching_without_gap_possible_with_tc(idxA,idxB,mat_pos_diag,seq_pos_to_be_matched) &&
    			mat(idx_i_diag,idx_j_diag)+add_score==mat(idx_i,idx_j)){
    		matching=true;
    	}

    	else if(matrixLR){

    		// if we are currently in matrix LR, we check whether the trace continues in
    		// matrix L
    		if(L(idx_i_diag,idx_j_diag)+add_score==mat(idx_i,idx_j)){
    			state=in_L; matching=true;
    		}

    		// if we are currently in matrix LR, we check whether the trace continues in
    		// matrix G_A (the only gap matrix)
    		else if(G_A(idx_i_diag,idx_j_diag)+add_score==mat(idx_i,idx_j)){
    			state=in_G_A; matching=true;
    		}
    	}

    	if(matching) cur_mat_pos = mat_pos_diag;

    	return matching;
    }


    // ---------------------------------------------------------------------------------------------------------
    // suboptimal

    // trace through the F matrix and L, G_A, G_AB and Lr matrices to find all maximally extended
    // suboptimal EPMs up to a certain threshold
    void ExactMatcher::trace_EPMs_suboptimal(){
    	subopt_score = this->min_subopt_score;
    	cout << "compute EPMs suboptimal " << endl;

    	// -----------------------------------------------------------------------------------
    	// for debugging

    	//const ArcMatch &am = arc_matches.arcmatch(128566); //example with many different arcmatches to proc
    	//const ArcMatch &am = arc_matches.arcmatch(215010); //example with high score
    	//const ArcMatch &am = arc_matches.arcmatch(61593); //96,for example 3
    	//const ArcMatch &am = arc_matches.arcmatch(39581);
    	//const Arc &a = bpsA.arc(4);
    	//const Arc &b = bpsB.arc(5);
    	//epm_cont_t found_EPMs;
    	//this->print_matrices_new(a,b,20,20,true);
    	//cout << "trace am " << am.idx() << " with score " << score_for_arc_match(am,true)<< endl;
    	//cout << "number of arcmatches " << arc_matches.num_arc_matches() << endl;
    	/*size_t max_size=0;
        	for(size_t i=0;i<(arc_matches.num_arc_matches());++i){
        		epm_cont_t found_EPMs;
        		const ArcMatch &am = arc_matches.arcmatch(i);
        		if(i%10000==0) cout <<  i << endl;
        		if(i>170000 && i%1000==0) cout << i << endl;
        		if(this->score_for_arc_match(am,true).is_finite()){
        			if(i%10000==0) cout << "trace am " << am.idx() << " with score " << score_for_arc_match(am,true)<< endl;
        			trace_LGLR_suboptimal(am,(infty_score_t)500,found_EPMs,true);
        			//int count=0;
        			if(found_EPMs.size()>max_size) max_size = found_EPMs.size();
        		}
        	}*/
    	//cout << "max number of EPMs found for arcmatch " << max_size << endl;

    	//trace_LGLR_suboptimal(a,b,500,found_EPMs,true);
    	//return;
    	//cout << "number of epms found " << found_EPMs.size() << endl;
    	//cout << "found epms " << found_EPMs << endl;

    	// -------------------------------------------------------------------------------------

    	cout << "maximal score " << F(pos_of_max.first,pos_of_max.second) << endl;
    	score_t max_in_F = F(pos_of_max.first,pos_of_max.second).finite_value();

    	subopt_score = max_in_F;

    	if (difference_to_opt_score > 0){
    		subopt_score = max_in_F-difference_to_opt_score;
    		cout << " Set suboptimal score to allowed difference" << endl;

    	} else if ( (subopt_range > 0.0) && (subopt_range < 1.0 ) ){
    		subopt_score = subopt_range * max_in_F;
    		cout << "trace EPMs within range " << subopt_range << " of best EPM score" << endl;
    	}

    	if (subopt_score < min_subopt_score) {
    		subopt_score = min_subopt_score;
    	}

    	cout << "subopt_score " << subopt_score << endl;

    	//compute length of best EPM
    	EPM best_epm;
    	if(debug) cout << "before trace F heuristic tc " << endl;

    	trace_F_heuristic(pos_of_max.first,pos_of_max.second,best_epm);

    	EPM::pat_vec_t::size_type max_length = best_epm.pat_vec_size();

    	cout << "max EPM length " << max_length << endl;

    	//as trace_F_heuristic_tc overwrites G_A!
    	this->initialize_gap_matrices();

    	double min_seq_length = static_cast<double>(std::min(seqA.length(),seqB.length()));

    	cout << "coverage " << max_length/min_seq_length << endl;

    	if(max_length/min_seq_length>cutoff_coverage){

    		//don't do suboptimal traceback if the coverage of the best EPM is above the coverage cutoff
    		best_epm.set_score(max_in_F);
    		add_foundEPM(best_epm);
    		return;
    	}

    	//compute actual suboptimal traceback
    	for(pos_type i=1;i<F.sizes().first;++i){
    		for(pos_type j=1;j<F.sizes().second;++j){

    			if(F(i,j)>=(infty_score_t)subopt_score){

    				if(i==F.sizes().first-1 || j==F.sizes().second-1 || seqA[i+1]!=seqB[j+1]){

    					score_t max_tol_left=F(i,j).finite_value()-subopt_score;
    					assert(max_tol_left>=0);

    					if(debug) cout << "trace position subopt " << i << "," << j << " with max tol left " << max_tol_left << endl;
    					// compute traceback from position (i,j)
    					trace_F_suboptimal(i,j,max_tol_left,true);
    				}
    			}
    		}
    	}
    }

    // traces through the F matrix from position (i,j) to find all suboptimal EPMs up to a certain
    // threshold that ends in (i,j)
    void ExactMatcher::trace_F_suboptimal(pos_type i,pos_type j,score_t max_tol,bool recurse){

    	if(debug) cout << "trace F suboptimal from pos " << i << "," << j << " with max tol " << max_tol << endl;
    	assert(F(i,j).is_finite());

    	epm_cont_t found_epms;
    	found_epms.push_back(EPM());
    	found_epms.back().set_max_tol_left(max_tol);

    	pos_type pos_cur_epm = 0;
    	pair_seqpos_t cur_pos = pair_seqpos_t(i,j);
    	map_am_to_do_t am_to_do_for_F;
    	infty_score_t cur_max_tol;
    	const PairArcIdx no_am(bpsA.num_bps(),bpsB.num_bps());
    	matpos_t dummy_pos(0,0); // the matrix position is not needed for tracing the F-matrix

    	poss_L_LR poss(-1,(infty_score_t)0,dummy_pos,no_am,cur_pos);

    	score_t min_allowed_score = F(i,j).finite_value()-max_tol;
    	assert(min_allowed_score>=0);

    	bool finished=false;

    	while(!finished){
    		while(!(F(cur_pos.first,cur_pos.second)==(infty_score_t)0)){

    			pos_type i= cur_pos.first;
    			pos_type j =cur_pos.second;
    			assert(i>=1);
    			assert(j>=1);

    			cur_max_tol = (infty_score_t)found_epms.at(pos_cur_epm).get_max_tol_left()-
    					F(i,j)+F(i-1,j-1)+score_for_seq_match();

    			// sequential matching
    			if(seqA[i]==seqB[j] && cur_max_tol>=(infty_score_t)0){

    				store_new_poss(pseudo_arcA,pseudo_arcB,false,poss_L_LR(in_F,cur_max_tol,dummy_pos,no_am,cur_pos),
    						poss,pos_cur_epm,found_epms,am_to_do_for_F);

    			}
    			// structural matching
    			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
    					arc_matches.common_right_end_list(i,j).end() != it; ++it ) {

    				const ArcMatch &am = arc_matches.arcmatch(*it);
    				const Arc &a = am.arcA();
    				const Arc &b = am.arcB();
    				const PairArcIdx pair_arcs(a.idx(),b.idx());

    				if(score_for_am(a,b)<FiniteInt(am_threshold)) continue;
    				cur_max_tol =(infty_score_t)found_epms.at(pos_cur_epm).get_max_tol_left()-
    						F(i,j)+F(am.arcA().left()-1,am.arcB().left()-1)+score_for_am(a,b);

    				if(cur_max_tol>=(infty_score_t)0){

    					store_new_poss(pseudo_arcA,pseudo_arcB,false,
    								poss_L_LR(in_F,cur_max_tol,dummy_pos,pair_arcs,cur_pos),
        							poss,pos_cur_epm,found_epms,am_to_do_for_F);
    				}
    			}
    			assert(poss.first!=-1); // we found at least one possibility

    			// store the first possibility
    			store_new_poss(pseudo_arcA,pseudo_arcB,true,poss,poss,pos_cur_epm,found_epms,am_to_do_for_F);

    			// update current position
    			EPM &cur_epm = found_epms.at(pos_cur_epm);
    			pair_seqpos_t last_matched_pos = cur_epm.last_matched_pos();
    			cur_pos = pair_seqpos_t(last_matched_pos.first-1,last_matched_pos.second-1);
    		}

    		finished=true;

    		// search for next epm to process (epm that is not traced completely)
    		for(;pos_cur_epm<found_epms.size();++pos_cur_epm){

    			pair_seqpos_t last_matched_pos = found_epms.at(pos_cur_epm).last_matched_pos();
    			if(!(F(last_matched_pos.first-1,last_matched_pos.second-1)==infty_score_t(0))){
    				finished=false;
    				cur_pos = pair_seqpos_t(last_matched_pos.first-1,last_matched_pos.second-1);
    				break;
    			}
    		}
    	}

    	// all epms are traced completely
    	// -> fill the missing parts (arcmatches which have been jumped over)
    	if(recurse){ //just for debugging!

    		preproc_fill_epm(am_to_do_for_F,pos_cur_epm,found_epms,min_allowed_score);
    	}

    	// check whether the EPM list doesn't contain duplicates and only maximally extended EPMs
    	assert(validate_epm_list(found_epms));

    	//debugging
    	cout << "number epms " << found_epms.size() << endl;
    	//if(found_epms.size()<1000) cout << "found epms " << found_epms << endl << endl;
    }

    // computes the suboptimal traceback through the L, G_A, G_AB and LR matrices
    // these matrices have to be recomputed
    void ExactMatcher::trace_LGLR_suboptimal(const Arc &a, const Arc &b,
    		score_t max_tol, epm_cont_t &found_epms, bool recurse){

    	if(debug_trace_LGLR_subopt) cout << "trace AGB suboptimal of am :" << a << "," << b << endl;
    	if(debug_trace_LGLR_subopt) cout << "with tol left " << max_tol << endl;
    	if(debug_trace_LGLR_subopt) cout << "D(a,b) " << D(a,b) << endl;

    	compute_LGLR(a,b,true); // recompute matrices L, G_A, G_AB and LR

    	if(debug_trace_LGLR_subopt) this->print_matrices(a,b,30,30,true);

    	assert(D(a,b).is_finite());

    	ArcIdx idxA = a.idx();
    	ArcIdx idxB = b.idx();

    	const PairArcIdx no_am(bpsA.num_bps(),bpsB.num_bps());
    	size_type num_posA = sparse_mapperA.number_of_valid_mat_pos(a.idx());
    	size_type num_posB = sparse_mapperB.number_of_valid_mat_pos(b.idx());

    	found_epms.push_back(EPM());
    	pos_type pos_cur_epm = 0;
    	pair_seqpos_t seq_pos_to_be_matched(a.right(),b.right());
    	poss_L_LR poss(-1,(infty_score_t)0,matpos_t(0,0),no_am,seq_pos_to_be_matched);

    	// a map is used to store the traces for all arcmatches that are part of the trace for the current arcmatch
    	// map: key: ArcIdx, mapped value: pair of max_tol_left and the result (traced EPMs for the ArcIdx)
    	// this map is updated while tracing the current arcmatch
    	map_am_to_do_t map_am_to_do;

    	// check whether traceback can be started in matrix G_A
    	poss_L_LR pot_new_poss(in_G_A,(G_A(num_posA-1,num_posB-1)-D(a,b)+max_tol),
    			               matpos_t(num_posA-1,num_posB-1),no_am,seq_pos_to_be_matched);
    	check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);

    	if(debug_trace_LGLR_subopt) cout << "checked in_G_A " << endl;

    	// check whether traceback can be started in matrix G_AB
    	pot_new_poss = poss_L_LR(in_G_AB,(G_AB(num_posA-1,num_posB-1)-D(a,b)+max_tol),
                matpos_t(num_posA-1,num_posB-1),no_am,seq_pos_to_be_matched);
    	check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);

    	if(debug_trace_LGLR_subopt) cout << "checked in G_AB " << endl;

    	// check whether traceback can be started in matrix LR
    	pot_new_poss = poss_L_LR(in_LR,(LR(num_posA-1,num_posB-1)-D(a,b)+max_tol),
				  matpos_t(num_posA-1,num_posB-1),no_am,seq_pos_to_be_matched);
    	bool LR_matching = check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);

    	if(debug_trace_LGLR_subopt) cout << "matching from LR? " << LR_matching << endl;

    	// if traceback cannot be started in matrix LR, we test whether it can be started in matrix L
    	if(!LR_matching){
    		pot_new_poss = poss_L_LR(in_L,(L(num_posA-1,num_posB-1)-D(a,b)+max_tol),
					  matpos_t(num_posA-1,num_posB-1),no_am,seq_pos_to_be_matched);
    		bool testL = check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);

    		if(debug_trace_LGLR_subopt) if(testL) cout << "beginning in L " << endl;
    	}

    	assert(poss.first!=-1); // we found at least one possibility
    	store_new_poss(a,b,true,poss,poss,pos_cur_epm,found_epms,map_am_to_do); // store first possibility

    	//todo: use unordered_map instead
    	bool finished=false;

    	while(!finished){

    		// continue traceback until we end up at pos (0,0) in matrix L or LR
    		while(found_epms.at(pos_cur_epm).get_cur_pos()!=matpos_t(0,0) ||
    				(found_epms.at(pos_cur_epm).get_state()!=in_L && found_epms.at(pos_cur_epm).get_state()!=in_LR)){

    			matpos_t cur_mat_pos = found_epms.at(pos_cur_epm).get_cur_pos();

    			pos_type idx_i = cur_mat_pos.first;
    			pos_type idx_j = cur_mat_pos.second;

    			assert(idx_i>=1);
    			assert(idx_j>=1);

    			bool matrixLR =  found_epms.at(pos_cur_epm).get_state()==in_LR;
    			const ScoreMatrix &mat = matrixLR ? LR : L;

    			assert(mat(idx_i,idx_j).is_finite());

    			seq_pos_to_be_matched = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,cur_mat_pos);
    			matpos_t mat_pos_diag;

    			if(debug_trace_LGLR_subopt) cout << "cur mat pos " << cur_mat_pos << endl;
    			if(debug_trace_LGLR_subopt) cout << "state " << found_epms.at(pos_cur_epm).get_state() << endl;

    			// sequential matching
    			if(sparse_trace_controller.pos_unpaired(idxA,idxB,matpos_t(idx_i,idx_j))){

    				mat_pos_diag = matpos_t(idx_i-1,idx_j-1); // we don't use the Trace Controller so far

    				// score contribution without the score for the next matrix position
    				score_t score_contr = found_epms.at(pos_cur_epm).get_max_tol_left()
       	       							  +score_for_seq_match()-mat(idx_i,idx_j).finite_value();

    				trace_seq_str_matching_subopt(a,b,score_contr,mat_pos_diag,seq_pos_to_be_matched,
    						no_am,poss,pos_cur_epm,found_epms,map_am_to_do);
    			}

    			// structural matching
    			for(ArcIdxVec::const_iterator itA=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).begin();
    					itA!=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i).end();++itA){
    				for(ArcIdxVec::const_iterator itB=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).begin();
    						itB!=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j).end();++itB){

    					const Arc &inner_a = bpsA.arc(*itA);
    					const Arc &inner_b = bpsB.arc(*itB);

    					if(score_for_am(inner_a,inner_b).is_neg_infty()) continue;

    					matidx_t pos_before_arcA = sparse_mapperA.first_valid_mat_pos_before(a.idx(),inner_a.left(),a.left());
    					matidx_t pos_before_arcB = sparse_mapperB.first_valid_mat_pos_before(b.idx(),inner_b.left(),b.left());

    					const PairArcIdx pair_arcs(*itA,*itB);

    					seq_pos_to_be_matched = pair_seqpos_t(inner_a.left(),inner_b.left());
    					mat_pos_diag=matpos_t(pos_before_arcA,pos_before_arcB);

    					assert(score_for_am(inner_a,inner_b).is_finite());

    					// score contribution without the score for the next matrix position
    					score_t score_contr = score_for_am(inner_a,inner_b).finite_value()
                           			          +score_for_stacking(a,b,inner_a,inner_b)
                           			          +found_epms.at(pos_cur_epm).get_max_tol_left()
                           			          -mat(idx_i,idx_j).finite_value();

    					trace_seq_str_matching_subopt(a,b,score_contr,mat_pos_diag,seq_pos_to_be_matched,
    							pair_arcs,poss,pos_cur_epm,found_epms,map_am_to_do);

    				}
    			}

    			assert(poss.first!=-1); // we found at least one possibility

    			// store the first possibility
    			store_new_poss(a,b,true,poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    		}

    		assert(found_epms.at(pos_cur_epm).get_cur_pos()==matpos_t(0,0) &&
    				(found_epms.at(pos_cur_epm).get_state()==in_L
    						|| found_epms.at(pos_cur_epm).get_state()==in_LR));

    		finished=true;

    		//search for next epm to process (epm that is not at pos(0,0))
    		for(;pos_cur_epm<found_epms.size();++pos_cur_epm){
    			if(found_epms.at(pos_cur_epm).get_cur_pos()!=matpos_t(0,0)){
    				finished=false;
    				break;
    			}
    		}
    	}

    	// all epms are traced completely in the current am (at pos(0,0))
    	// -> fill the missing parts (arcmatches which have been jumped over)
    	if(recurse){ //just for debugging

    		preproc_fill_epm(map_am_to_do,pos_cur_epm,found_epms);
    	}

    	//sort the epms according to the tolerance left in ascending order
    	sort(found_epms.begin(),found_epms.end());
    }

    // traces a sequential or structural match for the suboptimal traceback
    void ExactMatcher::trace_seq_str_matching_subopt(const Arc &a, const Arc &b, score_t score_contr,
    		matpos_t mat_pos_diag, pair_seqpos_t seq_pos_to_be_matched, const PairArcIdx am,
    		poss_L_LR &poss, size_type pos_cur_epm, epm_cont_t &found_epms, map_am_to_do_t &map_am_to_do){

    	bool matrixLR =  found_epms.at(pos_cur_epm).get_state()==in_LR;
    	const ScoreMatrix &mat = matrixLR ? LR : L;

    	matidx_t idx_i_diag = mat_pos_diag.first;
    	matidx_t idx_j_diag = mat_pos_diag.second;

    	seqpos_t i = seq_pos_to_be_matched.first;
    	seqpos_t j = seq_pos_to_be_matched.second;

    	bool matching_in_cur_mat = false; // whether a matching in the current matrix is possible

    	poss_L_LR pot_new_poss(found_epms.at(pos_cur_epm).get_state(),
    				mat(idx_i_diag,idx_j_diag)+score_contr,mat_pos_diag,am,seq_pos_to_be_matched);

    	// check if matching can be continued in the current matrix
    	if(sparse_mapperA.matching_without_gap_possible(a,i) &&
    			sparse_mapperB.matching_without_gap_possible(b,j)) {

    		matching_in_cur_mat = check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    	}

    	if(matrixLR){ // if we are in matrix LR

    		if(!matching_in_cur_mat){ // check if matching can be directly continued in matrix L
    			pot_new_poss = poss_L_LR(in_L,L(idx_i_diag,idx_j_diag)+score_contr,mat_pos_diag,am,seq_pos_to_be_matched);
    			check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    		}

    		// check if traceback can be continued in matrix G_A
    		pot_new_poss = poss_L_LR(in_G_A,G_A(idx_i_diag,idx_j_diag)+score_contr,mat_pos_diag,am,seq_pos_to_be_matched);
    		check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);

    		// check if traceback can be continued in matrix G_AB
    		pot_new_poss = poss_L_LR(in_G_AB,G_AB(idx_i_diag,idx_j_diag)+score_contr,mat_pos_diag,am,seq_pos_to_be_matched);
    		check_poss(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    	}
    }

    // checks whether the new possibility pot_new_poss is valid
    bool ExactMatcher::check_poss(const Arc &a, const Arc &b, const poss_L_LR &pot_new_poss,
    		poss_L_LR &poss, size_type pos_cur_epm, epm_cont_t &found_epms,
    		map_am_to_do_t &map_am_to_do){

    	// the possibility is valid if the tolerance left is at least 0
    	if(pot_new_poss.second>=(infty_score_t)0){

    		if(pot_new_poss.first==in_G_A || pot_new_poss.first==in_G_AB){

    			// trace through the gap-matrices if we are in a gap state
    			trace_G_suboptimal(a,b,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    		}

    		// otherwise store the new possibility
    		else store_new_poss(a,b,false,pot_new_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
    		return true;
    	}
    	return false;
    }

    // stores the new possibility new_poss
    void ExactMatcher::store_new_poss(const Arc &a, const Arc &b, bool last_poss,
    		const poss_L_LR &new_poss, poss_L_LR &poss, size_type pos_cur_epm,
    		epm_cont_t &found_epms, map_am_to_do_t &map_am_to_do){

    	assert(new_poss.first==in_L || new_poss.first==in_LR || new_poss.first==in_F);
    	assert(new_poss.second.is_finite());

    	// if it is the first possibility, we store it in poss and add it later to the
    	// epm (last_poss is true)
    	if(poss.first==-1){ poss = new_poss;}

    	else{

    		const score_t &max_tol = new_poss.second.finite_value();

    		//the epm that is currently processed
    		EPM *cur_epm = &found_epms.at(pos_cur_epm);

    		// if it is not the last possibility, copy the current epm and add
    		// the subsequent extension to the copied epm and reset pointer
    		if(!last_poss){
    			found_epms.push_back(*cur_epm);
    			cur_epm = &found_epms.back();
    		}

    		// sequential match
    		if(new_poss.fourth == PairArcIdx(bpsA.num_bps(),bpsB.num_bps())){
    			const pair_seqpos_t &cur_pos_seq = new_poss.fifth;

    			if(new_poss.first==in_F || cur_pos_seq!=pair_seqpos_t(a.right(),b.right())){
    				// store the sequence positions of the match
    				cur_epm->add(cur_pos_seq.first,cur_pos_seq.second,'.');
    			}
    		}

    		// structural match
    		else{

    			// adds the right and left ends of the arc match and
    			// stores the index of the arc match
    			const PairArcIdx &pair_arc_idx = new_poss.fourth;
    			const Arc &inner_a = bpsA.arc(pair_arc_idx.first);
    			const Arc &inner_b = bpsB.arc(pair_arc_idx.second);
    			cur_epm->add_and_store_am(inner_a,inner_b);

    			// construct map that stores the result of each used
    			// arc match during the traceback
    			const el_map_am_to_do_t &new_mapped_el = el_map_am_to_do_t(max_tol,epm_cont_t());

    			// try to insert the new element in the map
    			pair<map_am_to_do_t::iterator,bool> result=map_am_to_do.insert(map_am_to_do_t::value_type(pair_arc_idx,new_mapped_el));

    			const bool &el_inserted=result.second; // whether element was successfully inserted
    			score_t &max_tol_stored = result.first->second.first;

    			// if inserted element already exists and the tolerance stored is smaller as the
    			// current tol, the tolerance that is stored in the map is updated
    			if((!el_inserted) && max_tol_stored<max_tol){
    				max_tol_stored=max_tol;
    			}
    		}

    		// update information of the current epm
    		cur_epm->set_cur_pos(new_poss.third);
    		cur_epm->set_state(new_poss.first);
    		cur_epm->set_max_tol_left(max_tol);

    		// reset poss for the next iteration if last_poss is true
    		if(last_poss){poss.first=-1;}
    	}
    }

	// traces through the G_A and G_B matrices for the suboptimal case
    // todo: change current possibility instead of adding and removing
	void ExactMatcher::trace_G_suboptimal(const Arc &a, const Arc &b,
			const poss_L_LR &pot_new_poss, poss_L_LR &poss, size_type pos_cur_epm,
			epm_cont_t &found_epms, map_am_to_do_t &map_am_to_do){

		if(debug_trace_G) cout << "trace G suboptimal " << endl;

		list<poss_in_G> poss_G;
		poss_G.push_back(poss_in_G(pot_new_poss.first,pot_new_poss.second,pot_new_poss.third)); // initialization

		while(!poss_G.empty()){

			poss_in_G &cur_poss = poss_G.front();
			infty_score_t max_tol = cur_poss.second;
			matidx_t idx_i = cur_poss.third.first;
			matidx_t idx_j = cur_poss.third.second;
			int cur_state = cur_poss.first;

			if(debug_trace_G) cout << "mat pos " << cur_poss.third << endl;

			switch(cur_state){

			case in_G_A:
			{

				assert(idx_i>0);

				if(G_A(idx_i,idx_j)-G_A(idx_i-1,idx_j)<=max_tol){
					// continue traceback in G_A
					if(debug_trace_G) cout << "G_A (1) " << endl;
					poss_G.push_back(poss_in_G(in_G_A,max_tol-(G_A(idx_i,idx_j)-G_A(idx_i-1,idx_j)),matpos_t(idx_i-1,idx_j)));
				}

				if(G_A(idx_i,idx_j)-L(idx_i-1,idx_j)<=max_tol){
					// continue traceback in matrix L
					if(debug_trace_G) cout << "G_A (2) " << endl;
					poss_L_LR cur_poss(in_L,max_tol-(G_A(idx_i,idx_j)-L(idx_i-1,idx_j)),matpos_t(idx_i-1,idx_j),
							pot_new_poss.fourth,pot_new_poss.fifth);

					// check whether the gap is valid
					if(is_valid_gap(a,b,cur_poss)){
						if(debug_trace_G) cout << "store epm in gap case (1)" << endl;
						// store possibility
						store_new_poss(a,b,false,cur_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
					}
				}
			}break;

			case in_G_AB:
			{
				assert(idx_j>0);

				if(G_AB(idx_i,idx_j)-G_AB(idx_i,idx_j-1)<=max_tol){
					// continue traceback in G_AB
					if(debug_trace_G) cout << "G_AB (1) " << endl;
					poss_G.push_back(poss_in_G(in_G_AB,max_tol-(G_AB(idx_i,idx_j)-G_AB(idx_i,idx_j-1)),matpos_t(idx_i,idx_j-1)));
				}

				if(G_AB(idx_i,idx_j)-G_A(idx_i,idx_j-1)<=max_tol){
					// continue traceback in G_A
					if(debug_trace_G) cout << "G_AB (2) " << endl;
					poss_G.push_back(poss_in_G(in_G_A,max_tol-(G_AB(idx_i,idx_j)-G_A(idx_i,idx_j-1)),matpos_t(idx_i,idx_j-1)));
				}

				if(G_AB(idx_i,idx_j)-L(idx_i,idx_j-1)<=max_tol){
					// continue traceback in L
					if(debug_trace_G) cout << "G_AB (3) " << endl;
					poss_L_LR cur_poss(in_L,max_tol-(G_AB(idx_i,idx_j)-L(idx_i,idx_j-1)),matpos_t(idx_i,idx_j-1),
							pot_new_poss.fourth,pot_new_poss.fifth);

					// check whether the gap is valid
					if(is_valid_gap(a,b,cur_poss)){
						if(debug_trace_G) cout << "store epm in gap case (2)" << endl;
						// store possibility
						store_new_poss(a,b,false,cur_poss,poss,pos_cur_epm,found_epms,map_am_to_do);
					}
				}
			}break;

			default:
			{
				cerr << "no gap state - something went wrong! " << endl;
			}break;

			}
			poss_G.pop_front(); // remove the current possibility from the list
		}
	}

	// checks whether an epm is maximally extended, i.e. if the gap is valid
    bool ExactMatcher::is_valid_gap(const Arc &a, const Arc &b,
    		const poss_L_LR &pot_new_poss){

    	if(debug_VG) cout << "is_valid_gap " << endl;

    	const ArcIdx idxA = a.idx();
    	const ArcIdx idxB = b.idx();

    	// the current position in the L matrix
    	matidx_t cur_idx_i = pot_new_poss.third.first;
    	matidx_t cur_idx_j = pot_new_poss.third.second;

    	// if the current matrix position in matrix L is in the last column or row, the EPM cannot
    	// be extended anymore
    	if(cur_idx_i == sparse_mapperA.number_of_valid_mat_pos(idxA)-1 ||
    	   cur_idx_j == sparse_mapperB.number_of_valid_mat_pos(idxB)-1){
    		if(debug_VG) cout << "last column or row in matrix L " << endl;
    		return true;
    	}

    	// the last pair of sequence positions that was matched from the right
    	const pair_seqpos_t &pos_left_LR_seq = pot_new_poss.fifth;

    	matidx_t idx_i_right_G = sparse_mapperA.first_valid_mat_pos_before(idxA,pos_left_LR_seq.first,a.left());
    	matidx_t idx_j_right_G = sparse_mapperB.first_valid_mat_pos_before(idxB,pos_left_LR_seq.second,b.left());

    	matpos_t pos_right_G = matpos_t(idx_i_right_G,idx_j_right_G); //the right-most matrix position in the gap matrix

    	assert(cur_idx_i<=pos_right_G.first);
    	assert(cur_idx_j<=pos_right_G.second);

    	// length of the gap on the matrix level for sequence A
    	pos_type length_gapA = idx_i_right_G-cur_idx_i;
    	// length of the gap on the matrix level for sequence B
    	pos_type length_gapB = idx_j_right_G-cur_idx_j;

    	// if the length of the gap is 0 in either dimension of the matrix
    	// (i.e. either in the rows or the columns of the matrix is no gap when we compare
    	// the current matrix position (in matrix L) and the right-most position in the
    	// G-matrix), the EPM cannot be extended
    	if(length_gapA==0 || length_gapB==0){
    		if(debug_VG) cout << "gap of length 0 " << endl;
    		return true;
    	}

    	// the next diagonal matrix position from the left
    	matpos_t next_pos_from_left = matpos_t(cur_idx_i+1,cur_idx_j+1);

    	assert(next_pos_from_left.first<sparse_mapperA.number_of_valid_mat_pos(idxA) &&
    		   next_pos_from_left.second<sparse_mapperB.number_of_valid_mat_pos(idxB));

    	// the next pair of sequence positions from the left
    	pair_seqpos_t next_pos_from_left_seq = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,next_pos_from_left);

    	// last sequence position that was matched from the left
    	pair_seqpos_t pos_right_L_seq = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,pot_new_poss.third);

    	pair_seqpos_t pos_for_arcs_left = pair_seqpos_t(pos_right_L_seq.first+1,pos_right_L_seq.second+1);


    	if(debug_VG) cout << "pos right LR seq " << pos_left_LR_seq << endl;
    	if(debug_VG) cout << "pos right G " << pos_right_G << endl;
    	if(debug_VG) cout << "next_pos_from_left " << next_pos_from_left << endl;
    	if(debug_VG) cout << "next_pos_from_left_seq " << next_pos_from_left_seq << endl;
    	if(debug_VG) cout << "pos left L seq " << pos_right_L_seq << endl;
    	if(debug_VG) cout << "pos for arcs left " << pos_for_arcs_left << endl;

    	//----------------------------------------------------------------------------------------------
    	// check for sequential extension on the right side
    	//----------------------------------------------------------------------------------------------

    	if(sparse_mapperA.matching_without_gap_possible(a,pos_left_LR_seq.first) &&
    	   sparse_mapperB.matching_without_gap_possible(b,pos_left_LR_seq.second)){

    		pair_seqpos_t pos_right_G_seq = sparse_trace_controller.get_pos_in_seq_new(idxA,idxB,pos_right_G);

    		if(seqA[pos_right_G_seq.first] == seqB[pos_right_G_seq.second] // matching nucleotides
    		   && sparse_trace_controller.pos_unpaired(idxA,idxB,pos_right_G) ){ // and unpaired

    			if(debug_VG) cout << "not maximally extended to the right sequentially " << endl;
    			return false;
    		}

    		//----------------------------------------------------------------------------------------------
    		// check for structural extension on the right side
    		//----------------------------------------------------------------------------------------------

    		for(ArcIdxVec::const_iterator itA=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i_right_G).begin();
    				itA!=sparse_mapperA.valid_arcs_right_adj(idxA,idx_i_right_G).end();++itA){
    			for(ArcIdxVec::const_iterator itB=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j_right_G).begin();
    					itB!=sparse_mapperB.valid_arcs_right_adj(idxB,idx_j_right_G).end();++itB){

    				const Arc &inner_a = bpsA.arc(*itA);
    				const Arc &inner_b = bpsB.arc(*itB);
    				if(!score_for_am(inner_a,inner_b).is_finite()) continue;

    				if(debug_VG) cout << "inner am " << inner_a << "," << inner_b << endl;
    				if(debug_VG) cout << "pos right L seq " << pos_right_L_seq << endl;

    				// if arc match fits into the gap, the epm is not maximally extended
    				if(inner_a.left()>pos_right_L_seq.first && inner_b.left()>pos_right_L_seq.second){

    					if(debug_VG) cout << " not maximally extended to the right structurally " << endl;
    					return false;
    				}
    			}
    		}
    	}

    	if(debug_VG) cout << "maximally extended to the right " << endl;

    	//----------------------------------------------------------------------------------------------
    	// check for sequential extension on the left side
    	//----------------------------------------------------------------------------------------------

    	if(sparse_mapperA.matching_without_gap_possible(a,next_pos_from_left_seq.first)
    			&& sparse_mapperB.matching_without_gap_possible(b,next_pos_from_left_seq.second) // matching without gap possible
    			&& seqA[next_pos_from_left_seq.first]==seqB[next_pos_from_left_seq.second] // matching nucleotides
    			&& sparse_trace_controller.pos_unpaired(idxA,idxB,next_pos_from_left)){ // unpaired

    		if(debug_VG) cout << "not maximally extended to the left sequentially " << endl;
    		return false;
    	}

    	if(debug_VG) cout << "maximally extended to the left sequentially " << endl;

    	//----------------------------------------------------------------------------------------------
    	// check for structural extension on the left side
    	//----------------------------------------------------------------------------------------------

    	for(ArcIdxVec::const_iterator itA=sparse_mapperA.valid_arcs_left_adj(a,pos_for_arcs_left.first).begin();
    			itA!=sparse_mapperA.valid_arcs_left_adj(a,pos_for_arcs_left.first).end();++itA){
    		for(ArcIdxVec::const_iterator itB=sparse_mapperB.valid_arcs_left_adj(b,pos_for_arcs_left.second).begin();
    				itB!=sparse_mapperB.valid_arcs_left_adj(b,pos_for_arcs_left.second).end();++itB){

    			const Arc &inner_a = bpsA.arc(*itA);
    			const Arc &inner_b = bpsB.arc(*itB);
    			if(!score_for_am(inner_a,inner_b).is_finite()) continue;

    			if(debug_VG) cout << "inner am " << inner_a << "," << inner_b << endl;
    			if(debug_VG) cout << "pos lfet LR seq " << pos_left_LR_seq << endl;

    			// if arc match fits into the gap, the epm is not maximally extended
    			if(inner_a.right()<pos_left_LR_seq.first && inner_b.right()<pos_left_LR_seq.second){
    				if(debug_VG) cout << "not maximally extended to the left structurally " << endl;
    				return false;
    			}
    		}
    	}
    	return true;
    }

    // preprocesses the filling of the missing parts of the arc matches of the epm
    void ExactMatcher::preproc_fill_epm(map_am_to_do_t &map_am_to_do,
    		size_type pos_cur_epm, epm_cont_t &found_epms,score_t min_allowed_score){

    	// compute trace for all arc matches, that have been
    	// encountered while tracing the current position
    	// store all possible EPMs in the map map_am_to_do
    	for(map_am_to_do_t::iterator it = map_am_to_do.begin();
    			it!=map_am_to_do.end();++it){

    		const Arc &inner_a = bpsA.arc(it->first.first);
    		const Arc &inner_b = bpsB.arc(it->first.second);
    		if(debug) cout << "\t recurse: trace am " << inner_a << ","
    				<< inner_b << " with tolerance " << it->second.first << endl;

    		const score_t &tol = it->second.first;
    		epm_cont_t &epms = it->second.second; // the list of traced EPMs will be stored in the map
    		trace_LGLR_suboptimal(inner_a,inner_b,tol,epms,true);
    	}

    	// fill the missing parts of the EPMs
    	pos_type end = found_epms.size(); // we only want to search through the original found_epms vector (all newly
    	// inserted EPMs have the missing parts already filled)

    	for(pos_cur_epm=0;pos_cur_epm<end;++pos_cur_epm){
    		if(!found_epms[pos_cur_epm].number_of_am() == 0){

    			// fill the missing parts of the EPM
    			//preproc_fill_epm(map_am_to_do, min_allowed_score,
    			//		pos_cur_epm, found_epms);
    			vector<const EPM*> epms_to_insert;

    			assert(pos_cur_epm<found_epms.size());
    			assert(found_epms.at(pos_cur_epm).number_of_am()>0);

    			size_type number_of_am = found_epms.at(pos_cur_epm).number_of_am();
    			epms_to_insert.resize(number_of_am);

    			//max_tol_left_up_to_pos(vec_idx) gives the maximal tolerance that is left up to arc match vec_idx-1
    			vector<score_t> max_tol_left_up_to_pos;
    			max_tol_left_up_to_pos.resize(number_of_am);
    			// initialize the first entry in max_tol_left_up_to_pos with the tolerance left for the current epm
    			max_tol_left_up_to_pos[0]=found_epms.at(pos_cur_epm).get_max_tol_left();

    			size_type vec_idx=0;
    			fill_epm(map_am_to_do, vec_idx, max_tol_left_up_to_pos, epms_to_insert, min_allowed_score,
    					pos_cur_epm, found_epms);

    			// insert the missing parts of the first possibility
    			for(PairArcIdxVec::const_iterator arc_pairs = found_epms.at(pos_cur_epm).am_begin();
    					arc_pairs!=found_epms.at(pos_cur_epm).am_end();++arc_pairs){

    				const epm_cont_t &cur_epm_list = map_am_to_do.find(*arc_pairs)->second.second;

    				// make sure that the optimal solution is inserted
    				assert(map_am_to_do.find(*arc_pairs)->second.first == cur_epm_list.begin()->get_max_tol_left());

    				found_epms.at(pos_cur_epm).insert(*cur_epm_list.begin());
    			}
    			found_epms.at(pos_cur_epm).clear_am_to_do(); //clear arc matches to do for the epm at position pos_cur_epm
    			// the tolerance stays the same as we inserted only optimal solutions!

    			// if we came from the F-Matirx (min_score !=-1), we set the score and add the current epm
    			// at position pos_cur_epm to the patternPairMap
    			if(min_allowed_score!=-1){
    				found_epms[pos_cur_epm].set_score(min_allowed_score+found_epms[pos_cur_epm].get_max_tol_left());
    				add_foundEPM(found_epms[pos_cur_epm]);
    			}
    		}
    	}
    }

    // fills the missing parts of the arc matches of the epm
    void ExactMatcher::fill_epm(const map_am_to_do_t &map_am_to_do, size_type vec_idx,
    		vector<score_t> &max_tol_left_up_to_pos, vector<const EPM*> &epms_to_insert,
    		score_t min_score, size_type pos_cur_epm, epm_cont_t &found_epms){

    	assert(pos_cur_epm<found_epms.size());
    	assert(found_epms.at(pos_cur_epm).number_of_am()>0);
    	assert(vec_idx<found_epms.at(pos_cur_epm).number_of_am());

    	// arc match that is filled in the current epm
    	const PairArcIdx &cur_arcs_idx =  found_epms.at(pos_cur_epm).get_am(vec_idx);

    	map_am_to_do_t::const_iterator res = map_am_to_do.find(cur_arcs_idx);
    	assert(res!=map_am_to_do.end());
    	const epm_cont_t &cur_epm_list = res->second.second; // corresponding epm list of current arc match

    	const score_t &tol_traced_for_cur_am = res->second.first; // tolerance for which the current arc match was traced

    	for(epm_cont_t::const_iterator epm_it = cur_epm_list.begin();epm_it!=cur_epm_list.end();++epm_it){

    		// tolerance that is required for the current epm of the current arc match
    		score_t tol_required_for_cur_am = tol_traced_for_cur_am - epm_it->get_max_tol_left();

    		// maximal tolerance left after taking into account the current arc match
    		score_t max_tol_left=max_tol_left_up_to_pos[vec_idx] - tol_required_for_cur_am;

    		// not valid possibility as the maximal tolerance is exceeded
    		// all subsequent epms in cur_epm_list have a lower tolerance left (epm list is sorted in descending order)
    		// -> we don't need to look at these possibilities
    		if(max_tol_left < 0){break;}

    		// store a pointer to the current epm in the current arc match such that
    		// one can insert it later into the current epm
    		epms_to_insert.at(vec_idx) = &(*epm_it);

    		// update maximal tolerance left up to pos vec_idx+1
    		if(vec_idx != found_epms[pos_cur_epm].number_of_am()-1){
    			max_tol_left_up_to_pos[vec_idx+1]=max_tol_left;
    		}

    		// if we haven't filled all arc matches, we go to the next arc match indexed by vec_idx+1 in the current epm
    		if(vec_idx+1<found_epms[pos_cur_epm].number_of_am()){

    			size_type next_vec_idx = vec_idx+1;
    			fill_epm(map_am_to_do, next_vec_idx, max_tol_left_up_to_pos,
    					epms_to_insert, min_score, pos_cur_epm, found_epms);
    		}

    		// if all arc matches are filled
    		else{

    			// if this is not the first insertion for the current epm
    			if(!found_epms.at(pos_cur_epm).get_first_insertion()){

    				// copy the current epm
    				found_epms.push_back(found_epms.at(pos_cur_epm));

    				// insert the parts for the missing arc matches
    				for(vector<const EPM*>::const_iterator epm_to_insert = epms_to_insert.begin();
    						epm_to_insert!=epms_to_insert.end();++epm_to_insert){
    					found_epms.back().insert(**epm_to_insert);
    				}

    				found_epms.back().set_max_tol_left(max_tol_left); // update tolerance left
    				found_epms.back().clear_am_to_do(); // delete arc matches to do

    				if(min_score!=-1){ // we came from the F-matrix

    					// set the final score of the epm
    					found_epms.back().set_score(min_score+max_tol_left);
    					// store epm also in the patternPairMap
    					add_foundEPM(found_epms.back());
    				}
    			}
    			else{
    				// first insertion, we skip this possibility and insert it later
    				found_epms.at(pos_cur_epm).set_first_insertion(false);
    			}
    		}
    	}
    }

    // ---------------------------------------------------------------------------------------------------------
    // for debugging

    // print the matrices in the condensed form
    void ExactMatcher::print_matrices(const Arc &a, const Arc &b, size_type offsetA,size_type offsetB,bool suboptimal){
    	size_type num_posA = sparse_mapperA.number_of_valid_mat_pos(a.idx());
    	size_type num_posB = sparse_mapperB.number_of_valid_mat_pos(b.idx());
    	if(offsetA>num_posA){offsetA=num_posA;}
    	if(offsetB>num_posB){offsetB=num_posB;}
    	cout << "am " << a << "," << b << endl;
    	cout << "number of pos A " << num_posA << endl;
    	cout << "number of pos B " << num_posB << endl;
    	cout << "L" << endl;
    	assert(num_posA>=offsetA);
    	assert(num_posB>=offsetB);
    	for(size_type i=0;i<offsetA;++i){
    		for(size_type j=0;j<offsetB;++j){
    			cout << L(i,j) << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    	cout << "G_A" << endl;
    	for(size_type i=0;i<offsetA;++i){
    		for(size_type j=0;j<offsetB;++j){
    			cout << G_A(i,j) << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    	if(suboptimal){
    		cout << "G_AB" << endl;
    		for(size_type i=0;i<offsetA;++i){
    			for(size_type j=0;j<offsetB;++j){
    				cout << G_AB(i,j) << " ";
    			}
    			cout << endl;
    		}
    		cout << endl;
    	}
    	cout << "LR" << endl;
    	for(size_type i=0;i<offsetA;++i){
    		for(size_type j=0;j<offsetB;++j){
    			cout << LR(i,j) << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    }

    // checks whether an epm is valid, i.e. only one gap per arc match etc.
    bool ExactMatcher::validate_epm(const EPM &epm_to_test){

    	EPM::pat_vec_t::size_type pat_vec_size = epm_to_test.pat_vec_size();

    	if(pat_vec_size==0){
    		if(debug) cout << "empty epm " << endl;
    		return true;
    	}

    	for(EPM::pat_vec_t::size_type i =0;i<epm_to_test.pat_vec_size();++i){

    		EPM::el_pat_vec cur_el = epm_to_test.pat_vec_at(i);

    		//two matched positions in the EPMs have the same nucleotide
    		if(seqA[cur_el.first]!=seqB[cur_el.second]){
    			cerr << "two matched positions have a different nucleotides " << endl;
    			return false;
    		}

    		//check whether all matched positions are valid due to the TraceController
    		if(sparse_trace_controller.get_delta()!= (size_type)-1 &&
    		   !sparse_trace_controller.is_valid(cur_el.first,cur_el.second)){
    			cerr << "not all positions are valid (TraceController) " << endl;
    			return false;
    		}

    		//check whether both patVecs are ascending
    		if(i<pat_vec_size-1 && (cur_el.first>=epm_to_test.pat_vec_at(i+1).first
    				|| cur_el.second>=epm_to_test.pat_vec_at(i+1).second)){
    			cerr << "pattern Vecs are not ascending " << endl;
    			return false;
    		}

    	}

    	//validate connectivity of the epm
    	for(int k=0;k<2;++k){

    		std::vector<EPM::pair_size_t_pat_vec > arcmatches_to_validate;
    		bool gap = true;
    		arcmatches_to_validate.push_back(EPM::pair_size_t_pat_vec(0,pat_vec_size-1));

    		while(arcmatches_to_validate.size()!=0){

    			EPM::pair_size_t_pat_vec part_under_am = arcmatches_to_validate.back();
    			arcmatches_to_validate.pop_back();
    			if(part_under_am != EPM::pair_size_t_pat_vec(0,pat_vec_size-1)) gap=false; //in the F matrix no gap is allowed

    			//go over the part under the am (including the right end!)
    			for(EPM::pat_vec_t::size_type i=part_under_am.first;i<=part_under_am.second;++i){

    				EPM::el_pat_vec cur_pat_vec = epm_to_test.pat_vec_at(i);
    				unsigned int cur_el = (k==0) ? epm_to_test.pat_vec_at(i).first : epm_to_test.pat_vec_at(i).second;

    				unsigned int el_before=0;
    				if(i>0)	el_before = (k==0) ? epm_to_test.pat_vec_at(i-1).first : epm_to_test.pat_vec_at(i-1).second;

    				if(cur_pat_vec.third=='.' || cur_pat_vec.third==')'){

    					if(i>0 && (el_before+1!=cur_el)){
    						if(gap){ cerr << "more than one gap in EPM (1) " << endl; return false;}
    						gap=true;
    					}

    				}
    				else if(cur_pat_vec.third=='('){

    					pos_type pos_after_left_end = i+1; //the first position that needs to be validated of the inner am

    					//check validity for '('
    					if(i>0 && (el_before+1!=cur_el)){
    						if(gap){cerr << "more than one gap in EPM (2) " << endl; return false;}
    						gap=true;
    					}

    					int balance = 1;

    					while(balance!=0){ //find corresponding closing bracket

    						++i;
    						assert(i<epm_to_test.pat_vec_size());
    						if(epm_to_test.pat_vec_at(i).third=='(') ++balance;
    						if(epm_to_test.pat_vec_at(i).third==')') --balance;
    					} // i is the position of the corresponding ')' bracket

    					//if there are positions to check in the inner arcmatch, we store the inner arcmatch to check later,
    					//including the position of the ')'
    					if(i>=pos_after_left_end) arcmatches_to_validate.push_back(EPM::pair_size_t_pat_vec(pos_after_left_end,i));

    					//after incrementation of the for-loop, the next position to be checked is the one after the ')'
    				}
    			}
    		}
    	}
    	return true;
    }

    // checks the validity of the epm list, i.e. that no epm is contained in another epm (all
    // epms are maximally extended)
    bool ExactMatcher::validate_epm_list(epm_cont_t &found_epms){

    	cout << "validate epm list " << found_epms.size() << endl;

    	//todo: remove again!
    	if(found_epms.size()>1000){
    		cout << "epm list too long for validation " << endl;
    		return true;
    	}

    	int counter=0;

    	for(epm_cont_t::const_iterator cur_epm = found_epms.begin();cur_epm!=found_epms.end();++cur_epm){

    		if(counter%1000==0) cout << counter << endl; // show progress
    		++counter;

    		//go through all epms from the current epm on
    		epm_cont_t::const_iterator it = cur_epm;

    		++it; //compare to all other epms after cur_epm

    		for(;it!=found_epms.end();++it){

    			bool result = false;

    			if(cur_epm->pat_vec_size()>= it->pat_vec_size()){
    				result = cur_epm->includes(*it);
    			}

    			else{result = it->includes(*cur_epm);}

    			if(result){
    				cerr << "one pattern is included in another one! " << endl;
    				cerr << *cur_epm << endl;
    				cerr << *it << endl;
    				return false;
    			}
    		}
    	}

    	return true;
    }


    //--------------------------------------------------------------------------
    // class PatternPair
    //    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
    //--------------------------------------------------------------------------
    void PatternPair::resetBounds()
    {
	insideBounds.clear();
    }

    void PatternPair::setOutsideBounds(intPPair myPPair)
    {
	outsideBounds = myPPair;
    };

    void PatternPair::addInsideBounds(intPPair myPPair)
    {
	insideBounds.push_back(myPPair);
    };

    void PatternPair::setEPMScore(int myScore)
    {
	score = myScore;
    };

    //const string& PatternPair::get_struct() const
   // {
	//return structure;
    //};

    //--------------------------------------------------------------------------
    // class PatternPairMap
    //    is able to manage a set of PatternPairs(EPMs), each with 2 SinglePatterns
    //--------------------------------------------------------------------------
    PatternPairMap::PatternPairMap()
    {
	idMap.clear();
	patternList.clear();
	patternOrderedMap.clear();
	minPatternSize = 100000;
    }

    PatternPairMap::~PatternPairMap()
    {
	idMap.clear();
	int size= patternList.size();
	for(int i=0; i<size; i++){
	    delete patternList.front();
	    patternList.pop_front();
	}
	patternList.clear();
	patternOrderedMap.clear();
    }

    void PatternPairMap::add(const string& id,
			     const SinglePattern& first,
			     const SinglePattern& second,
			     const string& structure,
			     int score
			     )
    {
	PatternPair* p= new PatternPair(id,first,second,structure,score);
	SelfValuePTR myP = SelfValuePTR(p);
	patternList.push_back(myP);
	idMap.insert(make_pair(id,myP));
	if (p->getSize() < minPatternSize)  { minPatternSize = p->getSize(); }
    }

    void PatternPairMap::add(const SelfValuePTR value)
    {
	SelfValuePTR myP = SelfValuePTR(new PatternPair(*value));
	patternList.push_back(myP);
	idMap.insert(make_pair(value->getId(),myP));
	if (myP->getSize() < minPatternSize)  { minPatternSize = myP->getSize(); }
    }

    void  PatternPairMap::makeOrderedMap()
    {
	patternOrderedMap.clear();
	for(patListITER i = patternList.begin();i!=patternList.end();++i)
	    {
		patternOrderedMap.insert(make_pair((*i)->getSize(),*i));
	    }
    }

    void PatternPairMap::updateFromMap()
    {
	if (!patternOrderedMap.empty())
	    {
		idMap.clear();
		patternList.clear();
		for (orderedMapITER i=patternOrderedMap.begin();i!=patternOrderedMap.end();++i)
		    {
			add(i->second);
		    }
	    }
    }
    const PatternPair& PatternPairMap::getPatternPair(const string& id)const
    {
	return *(idMap.find(id)->second);
    }

    const    PatternPairMap::SelfValuePTR  PatternPairMap::getPatternPairPTR(const string& id)const
    {
	return (idMap.find(id)->second);
    }

    const PatternPairMap::patListTYPE& PatternPairMap::getList()const
    {
	return patternList;
    }
    const PatternPairMap::orderedMapTYPE& PatternPairMap::getOrderedMap() const
    {
	return patternOrderedMap;
    }

    PatternPairMap::orderedMapTYPE& PatternPairMap::getOrderedMap2()
    {
	return patternOrderedMap;
    }

    const int PatternPairMap::size() const
    {
	return idMap.size();
    }

    int  PatternPairMap::getMapBases()
    {
	int bases = 0;
	for(patListITER i = patternList.begin();i!=patternList.end();++i)
	    {
		bases += (*i)->getSize();
	    }
	return bases;
    }

    int  PatternPairMap::getMapEPMScore()
    {
	int EPMscore = 0;
	for(patListITER i = patternList.begin();i!=patternList.end();++i)
	    {
		EPMscore += (*i)->getEPMScore();
	    }
	return EPMscore;
    }

    std::ostream &operator << (std::ostream &out, const PatternPairMap::patListTYPE &pat_pair_list){
    	size_type i=0;
    	out << "epm_id\t score\t structure\t positions" << endl;
    	for(PatternPairMap::patListCITER it = pat_pair_list.begin(); it != pat_pair_list.end();++it,++i){
    		const PatternPair &pat_pair = **it;
    		out << i << "\t" << pat_pair.getScore() << "\t" << pat_pair.get_struct() <<  "\t"; //todo: fix
    		const intVec &pat1 = pat_pair.getFirstPat().getPat();
    		const intVec &pat2 = pat_pair.getSecPat().getPat();
    		intVec::const_iterator it_pat1 = pat1.begin();
    		intVec::const_iterator it_pat2 = pat2.begin();
    		for(; it_pat1 != pat1.end(), it_pat2 != pat2.end(); ++it_pat1, ++it_pat2){
    			out << *it_pat1 <<  ":" << *it_pat2 << " ";
    		}
    		out << std::endl;
    	}
    	return out;
    }

    LCSEPM::~LCSEPM()
    {
	//cout << endl << " execute destructor..." << endl;

	EPM_Table2.clear();
	// todo: delete pointers in EPM_table
	holeOrdering2.clear();
    }

    void    LCSEPM::calculateLCSEPM()
    {
	cout << " LCSEPM preprocessing..."  <<endl;
	cout << "     found #EPMs = " << patterns.size() << endl;
	cout << "    min EPM size = "<< patterns.getMinPatternSize()<< endl;
 	preProcessing();
	cout << " LCSEPM calculate holes..."  <<endl;
	cout << "   holes to calculate = " << holeOrdering2.size() << endl;
	calculateHoles3();
	cout << " LCSEPM calculate outmost D_rec..."  <<endl;
	int i = 1;
	int k = 1;
	vector < vector<int> > last_vec;
	int LCSEPMscore = D_rec2(i,seqA.length(),k,seqB.length(),last_vec,false);
	cout << "    Score LCS-EPM: "<< LCSEPMscore <<endl;
	cout << " LCSEPM calculate traceback..."  <<endl;
	calculateTraceback2(i,seqA.length(),k,seqB.length(),last_vec);
	int LCSEPMsize = matchedEPMs.getMapBases();
	cout << "    #EPMs: "<< matchedEPMs.size() << " / matched Bases: "<< LCSEPMsize <<endl;
    }

    void    LCSEPM::calculatePatternBoundaries(PatternPair*   myPair)
    {
	const vector<unsigned int>& myPatStr1 = myPair->getFirstPat().getPat();
	const vector<unsigned int>& myPatStr2 = myPair->getSecPat().getPat();

	myPair->resetBounds();

	for (unsigned int k=1;k < (myPatStr1.size());++k)
	    {
		if ( (myPatStr1[k]-patterns.getMinPatternSize() > myPatStr1[k-1])  &&
		     (myPatStr2[k]-patterns.getMinPatternSize() > myPatStr2[k-1]) ) {
		    myPair->addInsideBounds(std::make_pair(make_pair(myPatStr1[k-1],myPatStr1[k]),make_pair(myPatStr2[k-1],myPatStr2[k])));
		}
	    }

	// insert global min/max of the pattern
	myPair->setOutsideBounds(make_pair(make_pair(myPatStr1.front(),myPatStr1.back()),make_pair(myPatStr2.front(),myPatStr2.back())));
    }

    void LCSEPM::preProcessing()
    {
	// set EPM_Table size
	EPM_Table2.resize(seqA.length()+1);
        for (unsigned int i = 0; i < EPM_Table2.size();++i)
	    EPM_Table2[i].resize(seqB.length()+1);

	for (PatternPairMap::patListCITER myPair = patterns.getList().begin(); myPair != patterns.getList().end(); ++myPair)
	    {
		calculatePatternBoundaries(*myPair);

		// add EPM to EPM_table
		EPM_Table2[(*myPair)->getOutsideBounds().first.second][(*myPair)->getOutsideBounds().second.second].push_back(*myPair);

		// add all inside Holes from current EPM to holeOrdering multimap, sorted by holes size and exact position
		for(IntPPairCITER h = (*myPair)->getInsideBounds().begin(); h != (*myPair)->getInsideBounds().end(); ++h)
		    {
			// insert hole in multimap
			intPPairPTR myH = &(*h);
			holeOrdering2.insert(make_pair(myH,*myPair));
		    }
	    }
    }


    int LCSEPM::D_rec2(const int& i,const  int& j,const int& k,const int& l,vector < vector<int> >& D_h,const bool debug)
    {

	// initialize D_h matrix with 0
	D_h.clear();
	D_h.resize(j - i + 2);
	for (unsigned int a = 0; a < D_h.size();++a)
	    D_h[a].resize(l - k + 2,0);

	for(unsigned int j_1 = 1; j_1 < (j-i+2); ++j_1)
	    for (unsigned int l_2 = 1; l_2 < (l-k+2); ++l_2)
		{
		    if (debug==true){
			//	cout << "debug " << j_1 << "," << l_2 << endl;
		    }
		    // check if EPMs ending at current position
		    if (EPM_Table2[i + j_1-1][k + l_2-1].size() == 0)
			{
			    D_h[j_1][l_2] = (D_h[j_1-1][l_2]>D_h[j_1][l_2-1])? D_h[j_1-1][l_2]:D_h[j_1][l_2-1] ;
			    // bug in old version? this is new: - No!
			    //D_h[j_1][l_2] = max3(D_h[j_1-1][l_2],D_h[j_1][l_2-1],D_h[j_1-1][l_2-1]) ;
			}
		    else
			{
			    // get list of all EPMS ending at current pos
			    vector<PatternPairMap::SelfValuePTR> EPM_list = EPM_Table2[i + j_1-1][k + l_2-1];
			    int maxScore = 0;

			    // iterate over all EPMS to get best score
			    for (vector<PatternPairMap::SelfValuePTR>::iterator myIter = EPM_list.begin(); myIter < EPM_list.end(); ++myIter){

				//cout << i+j_1-1 << "," << k+l_2-1 << " patid: " <<  (*myIter)->getId() << endl;

				int pos_before_EPM_Str1 = (*myIter)->getOutsideBounds().first.first  - i;
				int pos_before_EPM_Str2 = (*myIter)->getOutsideBounds().second.first - k;

				int	score_EPM = 0;

				// check if EPM fits into cuurent hole
				if ((pos_before_EPM_Str1 >= 0) && (pos_before_EPM_Str2 >= 0)){
				    score_EPM = D_h[pos_before_EPM_Str1][pos_before_EPM_Str2] + (*myIter)->getScore();
				    //cout << (*myIter)->getId() << " FITS - EPM max score "<< score_EPM << " before " << pos_before_EPM_Str1+i <<","<< pos_before_EPM_Str2+k << " " << D_h[pos_before_EPM_Str1][pos_before_EPM_Str2] <<  endl;
				}


				if (score_EPM > maxScore) { maxScore = score_EPM; }
				//cout << (*myIter)->getId() << " EPM max score "<< score_EPM << " before " << pos_before_EPM_Str1+i <<","<< pos_before_EPM_Str2+k <<  endl;
			    }
			    //cout << "score hole max "<< maxScore << endl;
			    D_h[j_1][l_2] = max3(maxScore,D_h[j_1-1][l_2],D_h[j_1][l_2-1]);
			}

		}
	return (D_h[j - i + 1][l - k + 1]);
    }


    void LCSEPM::calculateHoles3()
    {
	intPPairPTR lastHole 			= NULL;
	//PatternPairMap::SelfValuePTR lastEPM 	= NULL;
	int lastHoleScore 			= 0;
	int skippedHoles			= 0;
	for (HoleMapCITER2 t = holeOrdering2.begin();t != holeOrdering2.end();++t)
	    {
		// check if current hole is exactly teh same as last hole
		// then we do not need to calculate again the same hole
		// ordering of "holeOrdering" ensures that similar holes are next to each other
		if ((lastHole == NULL) || (lastHole->first.first   != (*t).first->first.first) ||
		    (lastHole->first.second  != (*t).first->first.second) ||
		    (lastHole->second.first  != (*t).first->second.first) ||
		    (lastHole->second.second != (*t).first->second.second) ) {

		    //cout << endl << (*t).second->getId() << endl <<  " new current hole " << (*t).first->first.first << "," << (*t).first->first.second;
		    //cout << " - " << (*t).first->second.first << "," << (*t).first->second.second << endl;
		    //cout << "score old " << (*t).second->getScore() << " " << (*t).second->get_struct() << endl;

		    // calculate best score of hole
		    bool deb=false;
		    vector < vector<int> > vec;
		    int holeScore = D_rec2((*t).first->first.first+1,(*t).first->first.second-1,(*t).first->second.first+1,(*t).first->second.second-1,vec,deb);
		    (*t).second->setEPMScore(	(*t).second->getScore() + holeScore );

		    //cout << "score new " << (*t).second->getScore() << endl;
		    
		    lastHole = (*t).first;
		    //lastEPM = (*t).second;
		    lastHoleScore = holeScore;
		} else{
		    // add score of last hole to current EPM
		    (*t).second->setEPMScore((*t).second->getScore() + lastHoleScore);
		    skippedHoles++;
		    //cout << endl << (*t).second->getId() << endl <<  " new current hole " << (*t).first->first.first << "," << (*t).first->first.second;
		    //cout << " - " << (*t).first->second.first << "," << (*t).first->second.second <<  " " << (*t).second->get_struct() << endl;
		    //cout << "score:"<< lastHoleScore << "-"<< (*t).second->getEPMScore() << "-" << (*t).second->getScore() << " - current hole is same as last hole. skip!" << endl;
		}
	    }
	cout << "   skipped holes = " << skippedHoles << endl;
    }


    void LCSEPM::calculateTraceback2(const int i,const  int j,const int k,const int l,vector < vector<int> > holeVec)
    {
	int j_1 = holeVec.size()-1;
	int l_2 = holeVec[0].size()-1;

	while ((j_1 >= 1)&&(l_2 >= 1) && (holeVec[j_1][l_2]>0) )
	    {
		//cout << "traceback " << i+j_1-1 <<","<< k+l_2-1 << " score: "<<holeVec[j_1][l_2] << endl;
		if (holeVec[j_1][l_2-1] == holeVec[j_1][l_2])
		    --l_2;
		else
		    if (holeVec[j_1-1][l_2] == holeVec[j_1][l_2])
			--j_1;
		    else
			{
			    // get all EPMs which end at (i + j_1-1,k + l_2-1)
			    vector<PatternPairMap::SelfValuePTR> EPM_list = EPM_Table2[i + j_1-1][k + l_2-1];


			    // over all EPMs which end at (i+j_1-1,k+l_2-1)
			    for (vector<PatternPairMap::SelfValuePTR>::iterator myIter = EPM_list.begin(); myIter < EPM_list.end(); ++myIter){
				//cout << "here " << (*myIter)->getId() << endl;

				// check if current EPM fits inside current hole
				int x1 = (*myIter)->getOutsideBounds().first.first - i;
				int x2 = (*myIter)->getOutsideBounds().second.first - k;
				if  ( ( x1 >= 0 ) && ( x2 >= 0)){
				    // check score

				    //cout << "(j_1,l_2)=(" << j_1<< "," << l_2 <<")  "<< i <<","<< j << "-" << k << "," << l << "  outsidebounds first " <<  (*myIter)->getOutsideBounds().first.first << "," << (*myIter)->getOutsideBounds().second.first << endl;
				    //cout << "score (j1,l2)="<< holeVec[j_1][l_2] << " score=" << (*myIter)->getScore() << "  EPM_score=" << (*myIter)->getEPMScore() << " before="<<holeVec[(*myIter)->getOutsideBounds().first.first-1][(*myIter)->getOutsideBounds().second.first-1];
				    //cout << " " << (*myIter)->getOutsideBounds().first.first-1 << "," << (*myIter)->getOutsideBounds().second.first-1<< endl;
				    int check = (*myIter)->getScore() + holeVec[x1][x2];
				    if (holeVec[j_1][l_2] == check){

					// add current EPM to traceback
					//cout << "added traceback EPM "<< (*myIter)->getId() << endl;
					matchedEPMs.add( *myIter );

					// recurse with traceback into all holes of best EPM
					for(IntPPairCITER h = (*myIter)->getInsideBounds().begin(); h != (*myIter)->getInsideBounds().end(); ++h)
					    {
						vector < vector<int> > tmpHoleVec;
						tmpHoleVec.clear();
						//cout << (*myIter)->getId() << " D_rec2 hole " << (*h).first.first+1 << "," << (*h).first.second-1 << "-" << (*h).second.first+1 << "," << (*h).second.second-1 << endl;
						int sc = D_rec2((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec,true);
						// call traceback only if there is an EPM within hole
						//cout << (*myIter)->getId() << "score "<< sc << " " << (*h).first.first+1 << "," << (*h).first.second-1 << "-" << (*h).second.first+1 << "," << (*h).second.second-1 << " hole traceback..." << endl;
						if (sc > 0) {
						    calculateTraceback2((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec);
						}
					    }
					// jump with traceback to position before EPM
					j_1 = ( (*myIter)->getOutsideBounds().first.first ) - i;
					l_2 = ( (*myIter)->getOutsideBounds().second.first) - k;
					break;
				    }
				} // if EPM fits hole
			    } // for

			}
	    }
    }

    char* LCSEPM::getStructure(PatternPairMap& myMap, bool firstSeq, int length)
    {
	char* s= (char*) space(sizeof(char) * (length+1));
	for(int i= 0; i<length; i++)
	    s[i]='.';
	intVec patternVec;
	string structure;
	char x;
	for (PatternPairMap::patListCITER i=myMap.getList().begin();i != myMap.getList().end();++i)
	    {
		if(firstSeq)
		    patternVec= (*i)->getFirstPat().getPat();
		else patternVec= (*i)->getSecPat().getPat();
		structure= (*i)->get_struct();
		for(int j= 0; j<(int)patternVec.size(); j++)
		    {
			if(structure[j]=='(')
			    x= '(';
			else if(structure[j]==')')
			    x= ')';
			else 
			    x= '.';
			s[patternVec[j]-1]= x;
		    }
    
	    }
	return s;
    }
    void LCSEPM::MapToPS(const string& sequenceA, const string& sequenceB, PatternPairMap& myMap, const string& file1, const string& file2)
    {
	string func_str="\
   /drawpattern {\n\
      /Panz pattern length def\n\
      0 1 pattern length 1 sub {\n\
         /i exch def\n\
         pattern i get\n\
         newpath\n\
         {\n\
            1 Panz div i mul 0 add 1 1 sethsbcolor\n\
            coor exch 1 sub get aload pop fsize 2.1 div 0 360 arc\n\
            fill\n\
         } forall\n\
      } for\n\
   } bind def\n\
   \n\
   /pattern [\n";
	string clus1_str,clus2_str;



	stringstream label1Str, label2Str;

	for (unsigned int i=1;i<=sequenceA.length();++i)
	    {
		if (i % 50 == 0)
		    label1Str << i << " 0.5 0.5 (" << i << ") Label\n";
	    }

	for (unsigned int i=1;i<=sequenceB.length();++i)
	    {
		if (i % 50 == 0)
		    label2Str << i << " 0.5 0.5 (" << i << ") Label\n";
	    }

	for (PatternPairMap::patListCITER i=myMap.getList().begin();i != myMap.getList().end();++i)
	    {
		intVec tmpvec1=(*i)->getFirstPat().getPat();
      
		clus1_str+="["+intvec2str(tmpvec1," ")+"]\n";

		intVec tmpvec2=(*i)->getSecPat().getPat();
     
		clus2_str+="["+intvec2str(tmpvec2," ")+"]\n";
	    }
	clus1_str+="] def\n\n";
	clus2_str+="] def\n\n";
	clus1_str=func_str+clus1_str;
	clus2_str=func_str+clus2_str;

   
	string psfilename = file1;
	string pos1= "drawpattern\ndrawbases\n";
	pos1+=label1Str.str();
   
	fold_constrained= 1;
	char* structure= getStructure(myMap,true,sequenceA.length());
	fold(upperCase(sequenceA).c_str(), structure);
   
	PS_rna_plot_a(const_cast<char*>(sequenceA.c_str()),
		      const_cast<char*>(structure),
		      const_cast<char*>(psfilename.c_str()),
		      const_cast<char*>(clus1_str.c_str()),
		      const_cast<char*>(pos1.c_str()));

	pos1= "drawpattern\ndrawbases\n";
	pos1+= label2Str.str();
   
	psfilename = file2;
	//free(structure);
	//structure= NULL;
	//if(base_pair){ free(base_pair); base_pair= NULL;}
	//free_arrays();
	structure= getStructure(myMap,false,sequenceB.length());
	fold(upperCase(sequenceB).c_str(), structure);
   
	PS_rna_plot_a(const_cast<char*>(sequenceB.c_str()),
		      const_cast<char*>(structure),
		      const_cast<char*>(psfilename.c_str()),
		      const_cast<char*>(clus2_str.c_str()),
		      const_cast<char*>(pos1.c_str()));
	//free(structure);
	//structure= NULL;
	//if(base_pair){ free(base_pair); base_pair= NULL;}
	//free_arrays();
    }

    void LCSEPM::output_locarna(const string& sequenceA, const string& sequenceB, const string& outfile){

	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	for (PatternPairMap::patListCITER i=matchedEPMs.getList().begin();i != matchedEPMs.getList().end();++i)
	    {
		positionsSeq1LCSEPM.insert(positionsSeq1LCSEPM.end(),(*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.insert(positionsSeq2LCSEPM.end(),(*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());
		//SinglePattern my1 = (*i)->getFirstPat();
		//my1.print();
		//my1 = (*i)->getSecPat();
		//my1.print();
	    }

	sort(positionsSeq1LCSEPM.begin(),positionsSeq1LCSEPM.end());
	sort(positionsSeq2LCSEPM.begin(),positionsSeq2LCSEPM.end());;

	for (unsigned int i=0;i<positionsSeq1LCSEPM.size();++i)
	    {
		matchingsLCSEPM.push_back(make_pair(positionsSeq1LCSEPM[i],positionsSeq2LCSEPM[i]));
	    }
	//string outname = "locarna_constraints_input.txt"; //"/home/radwan/Exparna_P/LocARNA/src/locarna_constraints_input.txt";
	ofstream outLocARNAfile (outfile.c_str());

	int last_edge_seq1,last_edge_seq2;
	last_edge_seq1=0;
	last_edge_seq2=0;

	string seq1_1,seq1_2,seq1_3,seq2_1,seq2_2,seq2_3;
	int edge = 100;

	for (vector<intPair>::iterator i_edge = matchingsLCSEPM.begin(); (i_edge != matchingsLCSEPM.end() && seq1_1.size()<seqA.length() && seq2_2.size()<seqB.length());++i_edge)
	    {
		//cout << "first: " << (*i_edge).first << " second: " << (*i_edge).second << endl;

		for (int i=last_edge_seq1+1;(i<(int)((*i_edge).first)&&seq1_1.size()<seqA.length());++i)
		    {
			seq1_1.push_back('.');
			seq1_2.push_back('.');
			seq1_3.push_back('.');
		    }

		for (int j=last_edge_seq2+1;(j<(int)((*i_edge).second)&&seq2_2.size()<seqB.length());++j)
		    {
			seq2_1.push_back('.');
			seq2_2.push_back('.');
			seq2_3.push_back('.');
		    }

		ostringstream edge_st_;
		edge_st_ << edge;
		string edge_st;
		edge_st = edge_st_.str();
		const char *c_str_edge = edge_st.c_str();

		seq1_1.push_back(c_str_edge[0]);
		seq1_2.push_back(c_str_edge[1]);
		seq1_3.push_back(c_str_edge[2]);

		seq2_1.push_back(c_str_edge[0]);
		seq2_2.push_back(c_str_edge[1]);
		seq2_3.push_back(c_str_edge[2]);

		++edge;

		last_edge_seq1= (*i_edge).first;
		last_edge_seq2 = (*i_edge).second;
	    }

	// end stuff
	for (size_type i=last_edge_seq1+1;i<=seqA.length() && seq1_1.size()<seqA.length();++i)
	    {
		seq1_1.push_back('.');
		seq1_2.push_back('.');
		seq1_3.push_back('.');
	    }

	for (size_type j=last_edge_seq2+1;j<=seqB.length() && seq2_1.size()<seqB.length();++j)
	    {
		seq2_1.push_back('.');
		seq2_2.push_back('.');
		seq2_3.push_back('.');
	    }

	seq1_1 += "#1";
	seq1_2 += "#2";
	seq1_3 += "#3";

	seq2_1 += "#1";
	seq2_2 += "#2";
	seq2_3 += "#3";

	outLocARNAfile << ">"<< seqA.seqentry(0).name() << endl << upperCase(sequenceA) << endl;
	outLocARNAfile << seq1_1 << endl << seq1_2 << endl << seq1_3 << endl;
	outLocARNAfile << ">"<<seqB.seqentry(0).name() << endl << upperCase(sequenceB) << endl;
	outLocARNAfile << seq2_1 << endl << seq2_2 << endl << seq2_3 << endl << endl;

	outLocARNAfile.close();


    }

    void	LCSEPM::output_clustal(const string& outfile_name)
    {
	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	for (PatternPairMap::patListCITER i=matchedEPMs.getList().begin();i != matchedEPMs.getList().end();++i)
	    {
		positionsSeq1LCSEPM.insert(positionsSeq1LCSEPM.end(),(*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.insert(positionsSeq2LCSEPM.end(),(*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());
	    }

	sort(positionsSeq1LCSEPM.begin(),positionsSeq1LCSEPM.end());
	sort(positionsSeq2LCSEPM.begin(),positionsSeq2LCSEPM.end());;

	for (unsigned int i=0;i<positionsSeq1LCSEPM.size();++i)
	    {
		matchingsLCSEPM.push_back(make_pair(positionsSeq1LCSEPM[i],positionsSeq2LCSEPM[i]));
	    }

	//string outname = ensOptions.out_dir + "/" + ensOptions.align_file;
	ofstream outfile (outfile_name.c_str());

	string seq1_aln,seq2_aln; //,seq1_aln_str,seq2_aln_str;

	int last_edge_seq1,last_edge_seq2;
	last_edge_seq1=0;
	last_edge_seq2=0;

	for (vector<intPair>::iterator i_edge = matchingsLCSEPM.begin(); i_edge != matchingsLCSEPM.end();++i_edge)
	    {
		for (size_type i=last_edge_seq1+1;i<(*i_edge).first;++i)
		    {
			seq1_aln.push_back(seqA[i][0]);
			seq2_aln.push_back('-');
			//seq1_aln_str.push_back(myMol1.getStructure(i));
			//seq2_aln_str.push_back('-');
		    }
		for (size_type j=last_edge_seq2+1;j<(*i_edge).second;++j)
		    {
			seq1_aln.push_back('-');
			seq2_aln.push_back(seqB[j][0]);
			//seq1_aln_str.push_back('-');
			//seq2_aln_str.push_back(myMol2.getStructure(j));
		    }

		seq1_aln.push_back(seqA[(*i_edge).first][0]);
		seq2_aln.push_back(seqB[(*i_edge).second][0]);

		//seq1_aln_str.push_back(myMol1.getStructure((*i_edge).first));
		//seq2_aln_str.push_back(myMol2.getStructure((*i_edge).second));

		last_edge_seq1= (*i_edge).first;
		last_edge_seq2 = (*i_edge).second;
	    }

	// for the part after the last edge
	for (size_type i=last_edge_seq1+1;i<=seqA.length();++i)
	    {
		seq1_aln.push_back(seqA[i][0]);
		seq2_aln.push_back('-');
		//seq1_aln_str.push_back(myMol1.getStructure(i));
		//seq2_aln_str.push_back('-');
	    }
	for (size_type j=last_edge_seq2+1;j<=seqB.length();++j)
	    {
		seq1_aln.push_back('-');
		seq2_aln.push_back(seqB[j][0]);
		//seq1_aln_str.push_back('-');
		//seq2_aln_str.push_back(myMol2.getStructure(j));
	    }

	outfile << "CLUSTAL W (1.83) multiple sequence alignment --- expaRNA 0.7.2 - exact pattern Alignment of RNA --- Score: "<< matchingsLCSEPM.size() << endl <<endl;

	string tmp1 = seqA.seqentry(0).name() +"      ";
	string tmp2 = seqB.seqentry(0).name() +"      ";
	if (tmp1.length() < tmp2.length())
	    tmp1.resize(tmp2.length(),' ');
	else
	    if (tmp2.length() < tmp1.length())
		tmp2.resize(tmp1.length(),' ');
	string tmp3;
	tmp3.resize(tmp1.length(),' ');
	//outfile << tmp3 << seq1_aln_str <<endl;
	outfile << endl;
	outfile << tmp1 << seq1_aln << endl;
	outfile << tmp2 << seq2_aln << endl;
	//outfile << tmp3 << seq2_aln_str << endl << endl;
	outfile.close();
    }
} //end namespace
