#include "scoring.hh"
#include "rna_data.hh"
#include "sequence.hh"
#include "alphabet.hh"
#include "arc_matches.hh"


#include <math.h>
#include <fstream>

//! vectors for unpaired probabilities
//! temporarily required by sigma_
//! initialized in precompute sigma, if needed
//!
std::vector<double> punA_tab;
std::vector<double> punB_tab;

Scoring::Scoring(const Sequence &seqA_,
		 const Sequence &seqB_,
		 const ArcMatches *arc_matches_,
		 const MatchProbs *match_probs_,
		 const ScoringParams *params_,
		 bool exp_scores
		 ):
    params(params_),
    arc_matches(arc_matches_),
    match_probs(match_probs_),
    bpsA(&arc_matches_->get_base_pairsA()),
    bpsB(&arc_matches_->get_base_pairsB()),
    seqA(seqA_),
    seqB(seqB_),
    lambda_(0)
{

#ifndef NDEBUG
    if (params->ribosum) {
	// check sequences
	seqA.checkAlphabet(Alphabet<char>((char*)"ACGUN-",6),true);
	seqB.checkAlphabet(Alphabet<char>((char*)"ACGUN-",6),true);
    }
#endif
    

    precompute_sigma();
    precompute_gapcost();
    precompute_weights();

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
    bpsA(s.bpsA),
    bpsB(s.bpsB),
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
	    punA_tab[i]=bpsA->prob_unpaired(i);
	}
	punB_tab.resize(lenB+1);
	for (size_type i=1; i<=lenB; ++i) {
	    punB_tab[i]=bpsB->prob_unpaired(i);
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

	const AliColumn &colA=seqA[ia];
	const AliColumn &colB=seqB[ib];
	
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
Scoring::precompute_weights(const BasePairs &bps,
			    size_type len,
			    std::vector<score_t> &weights,
			    std::vector<score_t> &stack_weights) {
    size_type s   = bps.num_bps();
    weights.resize(s);
    if (params->stacking) stack_weights.resize(s);
    for (size_type i=0; i<s; i++) {
	const Arc &a=bps.arc(i);
	
	double p = bps.get_arc_prob(a.left(),a.right());
	weights[i] = probToWeight(p,len);
	
	if (params->stacking) {
	    if (bps.get_arc_prob(a.left()+1,a.right()-1)>0) {
		double stack_p = bps.get_arc_stack_prob(a.left(),a.right());
		stack_weights[i] = probToWeight(stack_p,len);
		// std::cout << i << " " << stack_p << " " << stack_weights[i] << std::endl;
	    }
	}
    }
}

void
Scoring::precompute_weights() {
    //score_t weight = 
    //score_t cond_weight = probToWeight(cond_prob);
    
    precompute_weights(*bpsA,seqA.length(),weightsA,stack_weightsA);
    precompute_weights(*bpsB,seqB.length(),weightsB,stack_weightsB);
}

//! convert probability to weight for scoring. In standard case,
//! return normlized log of probability. In case of mea,
//! return score_res*probability.
score_t
Scoring::probToWeight(double p, size_t len) const
{
    if (params->mea_scoring) {
	return round2score(params->probability_scale*p);
    } else {
	double pe=prob_exp(len);
	
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

double Scoring::prob_exp(size_t len) const {
    if (params->exp_prob >= 0) {
	return  params->exp_prob;
    } else {
	return prob_exp_f(len);
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
	const AliColumn &colA=seqA[i];
	for (size_type k=0; k<colA.size(); k++) {
	    gapfreqA[i] += (colA[k]=='-')?1:0;
	}
	gapfreqA[i] /= colA.size();
    }

    for (size_type i=1; i<lenB+1; i++) {
	const AliColumn &colB=seqB[i];
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
    
	
    const size_type rowsA = seqA.get_rows();
    const size_type rowsB = seqB.get_rows();
    
    double score=0;
    int gapless_combinations=0;
    
    const Alphabet<char> &alphabet = ribosum->alphabet();
    
    // compute geometric mean
    for(size_type i=0; i<rowsA; i++) { // run through all combinations of rows in A and B
	for(size_type j=0; j<rowsB; j++) { 
	    // wie sollen gaps behandelt werden?
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
    
    const size_type rowsA = seqA.get_rows();
    const size_type rowsB = seqB.get_rows();
	
    double score=0;
    int gapless_combinations=0;

    const Alphabet<char> &alphabet = ribosum->alphabet();
    
    for(size_type i=0; i<rowsA; i++) { // run through all combinations of rows in A and B
	for(size_type j=0; j<rowsB; j++) { 
	    // wie sollen gaps behandelt werden?
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
			    ( ribosum->basepair_prob(seqA[arcA.left()][i],seqA[arcA.right()][i])
			      *
			      ribosum->basepair_prob(seqB[arcB.left()][j],seqB[arcB.right()][j]))
			    );
		} else {
		    // score += 0.0; // undetermined nucleotides   
		}
	    }
	}
    }
	
    return round2score(100.0 * score / gapless_combinations);
}


score_t 
Scoring::arcmatch(const Arc &arcA, const Arc &arcB, bool stacked) const {
    // assert: if stacking score requested, inner arcs must have probability > 0,
    // and there must be a non-zero joint probability of the arc and the inner arc in each RNA 
    
    assert(!stacked || (bpsA->get_arc_prob(arcA.left()+1,arcA.right()-1)>0 && bpsB->get_arc_prob(arcB.left()+1,arcB.right()-1)>0));
    assert(!stacked || (bpsA->get_arc_2_prob(arcA.left(),arcA.right())>0 && bpsB->get_arc_2_prob(arcB.left(),arcB.right())>0));
    
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
	     );
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
	       (bpsA->get_arc_stack_prob(arcA.left(),arcA.right())
		+ bpsB->get_arc_stack_prob(arcB.left(),arcB.right()))			    
	       :
	       (bpsA->get_arc_prob(arcA.left(),arcA.right())
		+ bpsB->get_arc_prob(arcB.left(),arcB.right()))
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
	     );
    }
}



score_t
Scoring::arcmatch(const ArcMatch &am, bool stacked) const {
    score_t score;
    if (arc_matches->explicit_scores()) { // will not take stacking into account!!!
	score = arc_matches->get_score(am);
    } else {	
	const  Arc &arcA = am.arcA(); 
	const  Arc &arcB = am.arcB(); 
	
	score = arcmatch(arcA,arcB,stacked);
    }

    return score -4*lambda_; // modify for normalized alignment
}
