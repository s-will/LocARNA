#include "alphabet.hh"
#include "sequence.hh"
#include "scoring.hh"
#include "rna_data.hh"
#include "arc_matches.hh"
#include "match_probs.hh"
#include "ribosum.hh"


#include <math.h>
#include <fstream>

namespace LocARNA {

    //! vector for unpaired probabilities temporarily required by
    //! sigma_
    //! @note initialized in precompute sigma, if needed
    std::vector<double> punA_tab;

    //! vector for unpaired probabilities
    //! @see punB_tab
    std::vector<double> punB_tab;

    Scoring::Scoring(const Sequence &seqA_,
		     const Sequence &seqB_,
		     const RnaData &rna_dataA_,
		     const RnaData &rna_dataB_,
		     const ArcMatches &arc_matches_,
		     const MatchProbs *match_probs_,
		     const ScoringParams &params_,
		     bool exp_scores
		     ):
	params(&params_),
	arc_matches(&arc_matches_),
	match_probs(match_probs_),
	rna_dataA(rna_dataA_),
	rna_dataB(rna_dataB_),
	seqA(seqA_),
	seqB(seqB_),
	lambda_(0)
    {

#ifndef NDEBUG
	if (params->ribosum) {
	    // check sequences
	    if (!seqA.checkAlphabet(Alphabet<char>((char*)"ACGUN-",6))
		||
		!seqB.checkAlphabet(Alphabet<char>((char*)"ACGUN-",6))
		) {
		std::cerr << "WARNING: unsupported sequence characters found." <<std::endl;
	    }
	}
#endif


	precompute_sigma();
	precompute_gapcost();
	precompute_weights();

	apply_unpaired_penalty();
	if (exp_scores) {
	    exp_indel_opening_score =
		boltzmann_weight(params->indel_opening);
	    precompute_exp_sigma();
	    precompute_exp_gapcost();
	}
    }


    /*
      not need: copy constructor is implicit
      Scoring(const Scoring &s):
      params(s.params),
      arc_matches(s.arc_matches),
      match_probs(s.match_probs),
      rna_dataA(s.rna_dataA),
      rna_dataB(s.rna_dataB),
      seqA(s.seqA),
      seqB(s.seqB),
      lambda_(s.lambda_),
      // precomputed tables:
      sigma_tab(s.sigma_tab),
      gapcost_tabA(s.gapcost_tabA),
      gapcost_tabB(s.gapcost_tabB),
      weightsA(s.weightsA),
      weightsB(s.weightsB),
      stack_weightsA(s.stack_weightsA),
      stack_weightsB(s.stack_weightsB),
      exp_sigma_tab(s.exp_sigma_tab),
      exp_indel_opening_score(s.exp_indel_opening_score),
      exp_gapcost_tabA(s.exp_gapcost_tabA),
      exp_gapcost_tabB(s.exp_gapcost_tabB)
      {
      }
    */


    void
    Scoring::subtract(std::vector<score_t> &v,score_t x) const {
	std::transform(v.begin(),v.end(),v.begin(),std::bind2nd(std::minus<score_t>(),x));
    }
    void
    Scoring::subtract(Matrix<score_t> &m,score_t x) const {
	m.transform(std::bind2nd(std::minus<score_t>(),x));
    }



    void
    Scoring::apply_unpaired_penalty() {


	// subtract unpaired_penalty from precomputed tables
	// * simga_tab
	// * gapcost_tabA
	// * gapcost_tabB

	subtract(sigma_tab, 2*params->unpaired_penalty);
	subtract(gapcost_tabA, params->unpaired_penalty);
	subtract(gapcost_tabB, params->unpaired_penalty);

    }


    void
    Scoring::modify_by_parameter(score_t lambda) {

	score_t delta_lambda = lambda-lambda_;

	lambda_ = lambda;

	// subtract delta_lambda from precomputed tables
	// * simga_tab
	// * gapcost_tabA
	// * gapcost_tabB

	subtract(sigma_tab,2*delta_lambda);
	subtract(gapcost_tabA,delta_lambda);
	subtract(gapcost_tabB,delta_lambda);

    }

    void Scoring::precompute_sigma() {
	size_type lenA = seqA.length();
	size_type lenB = seqB.length();

	sigma_tab.resize(lenA+1,lenB+1);

	// precompute the unpaired probabilities and store in vectors
	if (params->mea_scoring) {
	    punA_tab.resize(lenA+1);
	    for (size_type i=1; i<=lenA; ++i) {
		punA_tab[i]=rna_dataA.prob_unpaired(i);
	    }
	    punB_tab.resize(lenB+1);
	    for (size_type i=1; i<=lenB; ++i) {
		punB_tab[i]=rna_dataB.prob_unpaired(i);
	    }
	}

	for (size_type i=1; i<=lenA; ++i) {
	    for (size_type j=1; j<=lenB; ++j) {
		sigma_tab(i,j) = sigma_(i,j);
	    }
	}
    }


    void Scoring::precompute_exp_sigma() {
	size_type lenA = seqA.length();
	size_type lenB = seqB.length();

	exp_sigma_tab.resize(lenA+1,lenB+1);

	for (size_type i=1; i<=lenA; ++i) {
	    for (size_type j=1; j<=lenB; ++j) {
		exp_sigma_tab(i,j) = boltzmann_weight(sigma_tab(i,j));
	    }
	}
    }


    /**
       returns similarity of two alignment columns
    */
    score_t
    Scoring::sigma_(int ia, int ib) const {

	if (params->mea_scoring) {

	    // for our mea alignment, we score basematchs
	    // by the sum of
	    //   the edge probability
	    //   and a term for structural accuracy

	    return
		round2score
		(
		 params->probability_scale *
		 ( match_probs->prob(ia,ib)
		   +
		   (params->alpha_factor/100.0)
		   * (
		      punA_tab[ia]
		      +
		      punB_tab[ib]
		      )
		   )
		 );
	} else {
	    // compute average score for aligning the two alignment columns

	    const Sequence::AliColumn &colA=seqA[ia];
	    const Sequence::AliColumn &colB=seqB[ib];

	    score_t score=0;

	    // the sum of pairs is quite ad hoc
	    // e.g. matching - and - counts as match ...
	    // N is a wildcard, matchs with N don't count
	    //
	    for (size_t i=0; i<colA.size(); i++) {
		for (size_t j=0; j<colB.size(); j++) {
		    // if we have a ribosum matrix, use it!
		    // !!! if ribosum and characters are not in the ribosum matrix
		    // we fall back to match/mismatch scoring !!!

		    //
		    if (params->ribosum
			&& params->ribosum->alphabet().in(colA[i])
			&& params->ribosum->alphabet().in(colB[j])) {

			score +=
			    round2score(100.0 * params->ribosum->basematch_score_corrected(colA[i],colB[j]));
		    } else {
			if (colA[i]!='N' && colB[j]!='N') {
			    score +=
				(colA[i]==colB[j])
				?params->basematch
				:params->basemismatch;
			}
		    }
		}
	    }

	    return  round2score(score / (int)(colA.size()*colB.size())) ;
	}
    }

    void
    Scoring::precompute_weights(const RnaData &rna_data,
				const BasePairs &bps,
				double exp_prob,
				std::vector<score_t> &weights,
				std::vector<score_t> &stack_weights) {
	size_type s   = bps.num_bps();
	weights.resize(s);
	if (params->stacking || params->new_stacking) stack_weights.resize(s);
	for (size_type i=0; i<s; i++) {
	    const Arc &a=bps.arc(i);

	    double p = rna_data.arc_prob(a.left(),a.right());
	    weights[i] = probToWeight(p,exp_prob);
	    
	    if (params->stacking) {
		if (rna_data.arc_prob(a.left()+1,a.right()-1)>0) {
		    double stack_p = rna_data.stacked_arc_prob(a.left(),a.right());
		    stack_weights[i] = probToWeight(stack_p,exp_prob);
		    // std::cout << i << " " << stack_p << " " << stack_weights[i] << std::endl;
		}
	    }
	    if (params->new_stacking) {
		if (!params->stacking) {
		    // if the old stacking is not used together with new one
		    stack_weights[i] = weights[i];
		}
		if (rna_data.arc_prob(a.left()+1,a.right()-1)>0) {
		    double joint_arc_p = rna_data.joint_arc_prob(a.left(),a.right());
		    stack_weights[i] += probToWeight(joint_arc_p,exp_prob);
		}
	    }
	}
    }

    void
    Scoring::precompute_weights() {
	//score_t weight = 
	//score_t cond_weight = probToWeight(cond_prob);
	
	precompute_weights(rna_dataA,arc_matches->get_base_pairsA(),params->exp_probA,weightsA,stack_weightsA);
	precompute_weights(rna_dataB,arc_matches->get_base_pairsB(),params->exp_probB,weightsB,stack_weightsB);
    }

    score_t
    Scoring::probToWeight(double p, double pe) const
    {
	if (params->mea_scoring) {
	    return round2score(params->probability_scale*p);
	} else {
	    return
		round2score(
			    round(params->struct_weight*
				  (1-log(p)/log(pe)))
			    );
	    /*
	      Score soll 0 sein für p=pe und struct_weight für p=1

	      ==> Normalisierung log(p/pe) / log(1/pe) 
	      //                  [=1-log(p)/log(pe)]

	      */
	}
    }

    // here we implement a special treatment of gap cost in multiple alignment
    // gap cost is varied position specific in order to favor gapping columns
    // that are already gapped.
    // This can be refined further!
    void
    Scoring::precompute_gapcost() {

	size_type lenA = seqA.length();
	size_type lenB = seqB.length();

	// resize and create tables
	gapcost_tabA.resize(lenA+1);
	gapcost_tabB.resize(lenB+1);
	std::vector<float> gapfreqA(lenA+1,0);
	std::vector<float> gapfreqB(lenB+1,0);

	// determine gap frequencies for each column in A and B

	for (size_type i=1; i<lenA+1; i++) {
	    const Sequence::AliColumn &colA=seqA[i];
	    for (size_type k=0; k<colA.size(); k++) {
		gapfreqA[i] += (colA[k]=='-')?1:0;
	    }
	    gapfreqA[i] /= colA.size();
	}

	for (size_type i=1; i<lenB+1; i++) {
	    const Sequence::AliColumn &colB=seqB[i];
	    for (size_type k=0; k<colB.size(); k++) {
		gapfreqB[i] += (colB[k]=='-')?1:0;
	    }
	    gapfreqB[i] /= colB.size();
	}
	// gap freqs determined

	// ----------------------------------------
	// compute position specific gap cost

	for (size_type i=1; i<lenA+1; i++) {
	    gapcost_tabA[i] = round2score((1-gapfreqA[i]) * params->indel);
	}

	for (size_type i=1; i<lenB+1; i++) {
	    gapcost_tabB[i] = round2score((1-gapfreqB[i]) * params->indel);
	}

    }

    void
    Scoring::precompute_exp_gapcost() {

	size_type lenA = seqA.length();
	size_type lenB = seqB.length();

	// resize and create tables
	exp_gapcost_tabA.resize(lenA+1);
	exp_gapcost_tabB.resize(lenB+1);

	for (size_type i=1; i<lenA+1; i++) {
	    exp_gapcost_tabA[i] = boltzmann_weight(gapcost_tabA[i]);
	}

	for (size_type i=1; i<lenB+1; i++) {
	    exp_gapcost_tabB[i] = boltzmann_weight(gapcost_tabB[i]);
	}
    }



    /*
      ATTENTION: handling of unknown nucleotide symbols (e.g. IUPAC) too simplistic
    */
    double
    Scoring::ribosum_arcmatch_prob(const Arc &arcA, const Arc &arcB) const {
	assert(params->ribosum != 0);
	// compute average ribosum score

	RibosumFreq *ribosum = static_cast<RibosumFreq *>(params->ribosum);


	const size_type rowsA = seqA.num_of_rows();
	const size_type rowsB = seqB.num_of_rows();

	double score=0;
	int gapless_combinations=0;

	const Alphabet<char> &alphabet = ribosum->alphabet();

	// compute geometric mean
	for(size_type i=0; i<rowsA; i++) { // run through all combinations of rows in A and B
	    for(size_type j=0; j<rowsB; j++) {
		// how to handle gaps?
		// current solution: ignore gap entries

		if (seqA[arcA.left()][i]!='-'
		    && seqA[arcA.right()][i]!='-'
		    && seqB[arcB.left()][j]!='-'
		    && seqB[arcB.right()][j]!='-') {

		    gapless_combinations++;

		    if (alphabet.in(seqA[arcA.left()][i])
			&& alphabet.in(seqA[arcA.right()][i])
			&& alphabet.in(seqB[arcB.left()][j])
			&& alphabet.in(seqB[arcB.right()][j])) {

			score+=
			    log(
				ribosum->arcmatch_prob( seqA[arcA.left()][i],seqA[arcA.right()][i],
							seqB[arcB.left()][j],seqB[arcB.right()][j] )
				/
				( ribosum->basematch_prob(seqA[arcA.left()][i],seqB[arcB.left()][j])
				  *
				  ribosum->basematch_prob(seqA[arcA.right()][i],seqB[arcB.right()][j]))
				);
		    }
		    else {
			// score += 0.0; // undetermined nucleotides
		    }
		}
	    }
	}

	return exp(score / gapless_combinations);
    }

    score_t
    Scoring::ribosum_arcmatch_score(const Arc &arcA, const Arc &arcB) const {
	assert(params->ribosum != 0);

	// compute average ribosum score

	RibosumFreq *ribosum = static_cast<RibosumFreq *>(params->ribosum);

	const size_type rowsA = seqA.num_of_rows();
	const size_type rowsB = seqB.num_of_rows();

	double score=0;
	int considered_combinations=0;

	const Alphabet<char> &alphabet = ribosum->alphabet();

	for(size_type i=0; i<rowsA; i++) { // run through all combinations of rows in A and B
	    for(size_type j=0; j<rowsB; j++) {
		// how to handle gaps?
		// current solution: ignore gap entries

		if (seqA[arcA.left()][i]!='-'
		    && seqA[arcA.right()][i]!='-'
		    && seqB[arcB.left()][j]!='-'
		    && seqB[arcB.right()][j]!='-') {

		    if (alphabet.in(seqA[arcA.left()][i])
			&& alphabet.in(seqA[arcA.right()][i])
			&& alphabet.in(seqB[arcB.left()][j])
			&& alphabet.in(seqB[arcB.right()][j])) {
			
			considered_combinations++;
			
			score+=
			    log(
				ribosum->arcmatch_prob( seqA[arcA.left()][i],seqA[arcA.right()][i],
							seqB[arcB.left()][j],seqB[arcB.right()][j] )
				/
				( ribosum->basepair_prob(seqA[arcA.left()][i],seqA[arcA.right()][i])
				  *
				  ribosum->basepair_prob(seqB[arcB.left()][j],seqB[arcB.right()][j]))
				)
			    / log(2);
		    } else {
			// score += 0.0; // undetermined nucleotides
		    }
		}
	    }
	}
	
	if (considered_combinations==0) return 0;
	
	return round2score(100.0 * score / considered_combinations);
    }


    score_t
    Scoring::arcmatch(const Arc &arcA, const Arc &arcB, bool stacked) const {
	// this method is disallowed with explicit arcmatch scores
	assert(!arc_matches->explicit_scores());

	// assert: if stacking score requested, inner arcs must have
	// probability > 0; moreover, there must be a non-zero joint
	// probability of the arc and the inner arc in each RNA
	assert(!stacked || 
	       (rna_dataA.arc_prob(arcA.left()+1,arcA.right()-1)>0
		&& rna_dataB.arc_prob(arcB.left()+1,arcB.right()-1)>0));
	
	assert(!stacked ||
	       (rna_dataA.joint_arc_prob(arcA.left(),arcA.right())>0 
		&& rna_dataB.joint_arc_prob(arcB.left(),arcB.right())>0));

	score_t sequence_contribution=0;


	// if tau is non-zero
	// we add a sequence contribution
	// for mea or if we don't have a ribosum matrix, we use sigma
	// otherwise we use a ribosum-like score 
	if ( params->tau_factor != 0 ) {
	    if (!params->mea_scoring && params->ribosum) {
		sequence_contribution = ribosum_arcmatch_score(arcA,arcB);
	    } else {
		sequence_contribution =
		    sigma_tab(arcA.left(),arcB.left())
		    + sigma_tab(arcA.right(),arcB.right())
		    + 2*lambda_; // when a scoring is modified by parameter
		// lambda, we undo the change of the
		// sequence contribution due to modified
		// sigma tabs.

	    }
	}

	if (! params->mea_scoring) {
	    return
		// base match contribution
		// (also for arc-match add terms for the base match on both ends, weighted by tau_factor)
		(params->tau_factor * sequence_contribution / 100)
		+
		// base pair weights
		(
		 stacked
		 ?
		 (stack_weightsA[arcA.idx()] + stack_weightsB[arcB.idx()])
		 :
		 (weightsA[arcA.idx()] + weightsB[arcB.idx()])
		 )
		 -4*lambda_
		;
	} else { //mea scoring
	    return
		static_cast<score_t>
		(
		 (
		  // base match contribution
		  // (also for arc-match add terms for the base match on both ends, weighted by tau_factor)
		  params->probability_scale *
		  (params->tau_factor/100.0) * sequence_contribution
		  +

#ifdef MEA_SCORING_OLD
		  params->probability_scale *
		  // structure contribution, weighted by beta
		  (params->beta_factor/100.0)
		  *
		  (
		   stacked
		   ?
		   (stack_weightsA[arcA.idx()] + stack_weightsB[arcB.idx()])
		   :
		   (weightsA[arcA.idx()] + weightsB[arcB.idx()])
		   )
		  +
		  // consensus contribution
		  params->probability_scale *
		  (params->gamma_factor/100.0) *
		  (
		   match_probs==NULL
		   ?
		   1.0
		   :
		   (
		    match_probs->prob(arcA.left(),arcB.left())
		    * match_probs->prob(arcA.right(),arcB.right())
		    )
		   )
		  *
		  ribosum_arcmatch_prob(arcA,arcB)

#else // NEW MEA SCORING
		  params->probability_scale *
		  // structure and consensus contribution, weighted by beta
		  (params->beta_factor/100.0)
		  *
		  (
		   stacked
		   ?
		   (rna_dataA.stacked_arc_prob(arcA.left(),arcA.right())
		    + rna_dataB.stacked_arc_prob(arcB.left(),arcB.right()))
		   :
		   (rna_dataA.arc_prob(arcA.left(),arcA.right())
		    + rna_dataB.arc_prob(arcB.left(),arcB.right()))
		   )
		  *
		  //
		  (
		   match_probs==NULL
		   ?
		   1.0
		   :
		   (
		    match_probs->prob(arcA.left(),arcB.left())
		    * match_probs->prob(arcA.right(),arcB.right())
		    )
		   )
		  *
		  ribosum_arcmatch_prob(arcA,arcB)
#endif
		  )
		 )
		 -4*lambda_;
	}
    }



    score_t
    Scoring::arcmatch(const ArcMatch &am, bool stacked) const {
	score_t score;
	if (arc_matches->explicit_scores()) { // does not take stacking into account!!!
	    score = arc_matches->get_score(am)  - 4*lambda_;
	} else {	
	    const  Arc &arcA = am.arcA();
	    const  Arc &arcB = am.arcB();

	    score = arcmatch(arcA,arcB,stacked);
	}

	return score; // modify for normalized alignment
    }

    // Very basic interface
    score_t
    Scoring::arcDel(const Arc &arcX, bool isA, bool stacked) const { //TODO Important Scoring scheme for aligning an arc to a gap is not defined and implemented!

	if (arc_matches->explicit_scores()) { // will not take stacking into account!!!
	    std::cerr << "ERROR sparse explicit scores is not supported!" << std::endl; //TODO: Supporting explicit scores for arcgap
	    assert( ! arc_matches->explicit_scores());
	}

	if (! params->mea_scoring) {
	    return
		// base pair weights
		loop_indel_score(gapX(arcX.left(), isA) + gapX(arcX.right(), isA)) + //score of aligining base-pair ends wo gap, such that it is multiply by gamma_loop ratio
	        (
		 stacked
		 ?
		 (isA? stack_weightsA[arcX.idx()] : stack_weightsB[arcX.idx()])
		 :
		 (isA? weightsA[arcX.idx()] : weightsB[arcX.idx()])
		 );
	}
	else
	    {
		std::cerr << "ERROR sparse mea_scoring is not supported!" << std::endl; //TODO: Supporting mea_scoring for arcgap
		assert( ! params->mea_scoring);
		return 0;
	    }
    }

    
    bool
    Scoring::is_stackable_arcA(const Arc &a) const {
	return
	    rna_dataA.joint_arc_prob(a.left(),a.right()) > 0;
    } 
    
    bool
    Scoring::is_stackable_arcB(const Arc &a) const {
	return
	    rna_dataB.joint_arc_prob(a.left(),a.right()) > 0;
    }

    bool
    Scoring::is_stackable_am(const ArcMatch &am) const {
	return
	    is_stackable_arcA(am.arcA())
	    &&
	    is_stackable_arcB(am.arcB());
    }
}
