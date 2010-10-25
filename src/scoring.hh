#ifndef SCORING_HH
#define SCORING_HH

#include <math.h>

#include "matrices.hh"
#include "basepairs.hh"
#include "rna_data.hh"
#include "ribosum.hh"
#include "match_probs.hh"

#include "infty_arith_int.hh"

typedef long int score_t;
//! an extended score_t that can store and calculate with
//! infinite values (i.p. we use -infty for invalid matrix entries)
typedef InftyArithInt infty_score_t;

typedef Matrix<infty_score_t> ScoreMatrix;
typedef std::vector<infty_score_t> ScoreVector;

//! type of partition functions

#ifdef LARGE_PF
  typedef long double pf_score_t;
#else
  typedef double pf_score_t;
#endif

typedef std::vector<pf_score_t> PFScoreVector;
typedef Matrix<pf_score_t> PFScoreMatrix;
typedef Matrix<double> ProbMatrix;



//#define MEA_SCORING_OLD

class RibosumFreq;
class Scoring;
class ArcMatches;
class ArcMatch;



//! Contains all parameters
//! for doing the scoring of alignments
//! in class Scoring.
//! The class encapsulates the configuration
//! of the score.
//!
class ScoringParams {
    friend class Scoring;

    //! constant bonus for a base match.
    //! together with basemismatch yields the simplest form
    //! of base match/mismatch scoring.
    //! Can be replaced by RIBOSUM scores (planned).
    score_t basematch;
    
    //! constant cost of a base mismatch
    score_t basemismatch;
   
    //! cost per indel (for linear or affine gap cost).
    score_t indel;

    //! cost per gap (for affine gap-cost). Use affine gap cost if non-zero.
    score_t indel_opening;
        
    //! the ribosum matrix, if non-null it is used
    //! for base match/mismatch instead of constant values
    //! and as contribution for arc-matchs (tau_factor)
    RibosumFreq *ribosum;
    
    //! Factor for structure contribution in the classic score.
    //! Maximal contribution of 1/2 arc match
    score_t struct_weight;
    
    //! Factor for the contribution of sequence score or
    //! ribosum score to arc matchs
    score_t tau_factor;
    
    //! cost of one exclusion.
    score_t exclusion;
    
    double exp_prob;
  
    double temperature;
    
    
    //! turn on/off stacking terms
    bool stacking;

    //! turn on/off mea scoring
    bool mea_scoring;

    //! weight for mea contribution "unstructured"
    score_t alpha_factor;

    //! weight for mea contribution "structure" 
    score_t beta_factor;

    //! weight for mea contribution "consensus"
    score_t gamma_factor;

    
    //! resolution of the mea score.
    //! Since the mea score is composed from
    //! probabilities and scores are integral
    //! there must be some factor to scale the probabilities
    //! in order to get scores.
    score_t probability_scale;
    
    /*
      The mea score is composed in the following way:
      
      params->probability_scale *
      (
       sum_basematchs (i,j)
           P(i~j)
           + alpha_factor/100 * (P(i unstructured) + P(j unstructured))
       + 
       sum_arcmatchs (a,b)
           beta_factor/100 * (P(a)+P(b)) * P(al~bl) * P(ar~br) * ribosum_arcmatch(a,b)
      )
    */
    
public:
    ScoringParams(score_t basematch_,
		  score_t basemismatch_,
		  score_t indel_,
		  score_t indel_opening_,
		  RibosumFreq *ribosum_,
		  score_t struct_weight_,
		  score_t tau_factor_,
		  score_t exclusion_,
		  double exp_prob_,
		  double temp_,
		  bool stacking_,
		  bool mea_scoring_,
		  score_t alpha_factor_,
		  score_t beta_factor_,
		  score_t gamma_factor_,
		  score_t probability_scale_
		  )
	: basematch(basematch_),
	  basemismatch(basemismatch_),
	  indel(indel_),
	  indel_opening(indel_opening_),
	  ribosum(ribosum_),
	  struct_weight(struct_weight_),
	  tau_factor(tau_factor_),
	  exclusion(exclusion_),
	  exp_prob(exp_prob_),
	  temperature(temp_),
	  stacking(stacking_),
	  mea_scoring(mea_scoring_),
	  alpha_factor(alpha_factor_),
	  beta_factor(beta_factor_),
	  gamma_factor(gamma_factor_),
	  probability_scale(probability_scale_)
	{
	}
};

//! Maintains and provides the scores for the alignment of two specific
//! Rnas.
//! Configurable via class ScoringParams, supports
//! the classic LocARNA-score as well as the LocARNA MEA-Score.
//! for "Classic Score" it supports
//!  * match/mismatch cost or RIBOSUM (planned) for bases
//!  * linear or affine gap cost
//!  * arc match scores build of arc probabilitiesmatch_probs
//!    and optionally base match contribution or RIBOSUM (planned)
//! for the "MEA Score" it supports
//!  * base match/mismatch scores derived from base match probabilities
//!  * RIBOSUM for base and arc match contribution
//!  * arc match score contributions from arc probabilities
//!  * weighting of single score contributions
//! 
//! The scoring class supports only characters ACGU, gap -, and N for ribosum scoring.
//! ALL other characters are treated as unknown characters and are basically ignored.
//! In consequence, IUPAC codes are handled badly. Even more important,
//! for ribosum scoring, sequences have to be upper case and Ts should be
//! converted to Us!!!
//! 
class Scoring {
public:
    
private:
    const ScoringParams *params;
    
    const ArcMatches *arc_matches;
    
    const MatchProbs *match_probs;

    const BasePairs *bpsA;
    const BasePairs *bpsB;
    const Sequence &seqA;
    const Sequence &seqB;
     
    score_t lambda_; // parameter for modified scoring in normalized
		    // local alignment
    
public:
    typedef Sequence::size_type size_type;

    //! @brief construct scoring object
    //!
    //! @param seqA first sequence
    //! @param seqB second sequence
    //! @param arc_matches the (significant) arc matches between the sequences
    //! @param ScoringParams a collection of parameters
    //! @param exp_scores only if true, the results of the exp_* scoring functions are defined, otherwise precomputations can be ommitted.
    //!
    Scoring(const Sequence &seqA,
	    const Sequence &seqB,
	    const ArcMatches *arc_matches,
	    const MatchProbs *match_probs,
	    const ScoringParams *params,
	    bool exp_scores=false
	    );

    /*
      implicit copy constructor
      Scoring(const Scoring &s);
    */
    
    //! modify scoring by a parameter lambda. Used in the
    //! computational of normalized local alignment by fractional
    //! programming in the form of the Dinkelbach algorithm.  The
    //! modification is destructive.  modify_by_parameter(0) results
    //! in unmodified scoring. Repeated modifications are valid.  The
    //! Outcome is *undefined* when the Scoring object is used to
    //! obtain exponentiated scores (methods exp_...) as in partition
    //! function alignment.  @param lambda the parameter
    void
    modify_by_parameter(score_t lambda);

    score_t lambda() const {return lambda_;}
    
private:
    // ------------------------------
    // tables for precomputed score contributions
    //
    Matrix<score_t> sigma_tab; //!< precomputed table of base match similarities 
    
    std::vector<score_t> gapcost_tabA; //!< table for gapcost in A
    std::vector<score_t> gapcost_tabB; //!< table for gapcost in B
    
    std::vector<score_t> weightsA; //<! weights of base pairs in A
    std::vector<score_t> weightsB; //<! weights of base pairs in B

    std::vector<score_t> stack_weightsA; //<! weights of stacked base
					 //<! pairs in A
    std::vector<score_t> stack_weightsB; //<! weights of stacked base
                                         //<! pairs in B
    
    // ------------------------------
    // tables for precomputed exp score contributions for partition function 
    //
    Matrix<pf_score_t> exp_sigma_tab; //!< precomputed table of exp base match similarities 
    pf_score_t exp_indel_opening_score; //!< precomputed value for exp of indel opening cost
    std::vector<pf_score_t> exp_gapcost_tabA; //!< table for exp gapcost in A
    std::vector<pf_score_t> exp_gapcost_tabB; //!< table for exp gapcost in B


    //! round a double to score_t.
    //! more robust than truncate
    score_t round2score(double d) const {
      return (score_t)((d<0) ? (d-0.5) : (d+0.5));
    }

    //! compute base similarity
    score_t
    sigma_(int i, int j) const;

    //! precompute all base similarities
    //! and put them in the sigma table
    void
    precompute_sigma();

    //! precompute all exp base similarities
    //! and put them in the exp_sigma table
    void
    precompute_exp_sigma();

    //! precompute the tables for gapcost
    void
    precompute_gapcost();

    //! precompute the tables for exp gapcost
    void
    precompute_exp_gapcost();
    
    //! compute weights/stacked weights for all arcs in A and B
    void
    precompute_weights();

    //! helper for precompute_weights (does job for one rna)
    void
    precompute_weights(const BasePairs &bps,
		       size_type len,
		       std::vector<score_t> &weights,
		       std::vector<score_t> &stack_weights);
    
    //! returns weight that corresponds to probability p
    //! in non-mea score
    score_t
    probToWeight(double p, size_type len) const;
    
    double
    prob_exp_f(int seqlen) const {return 1.0/(2.0*seqlen);}//!<expected probability of a base pair (null-model)
    

    //! returns probability of matching the concrete bases in an arcmatch,
    //! probability is based on ribosum data
    double
    ribosum_arcmatch_prob(const Arc &arcA, const Arc &arcB) const;    
    
    //! returns score of matching the concrete bases in an arcmatch
    //! based on ribosum data (only ribosum contribution)
    score_t
    ribosum_arcmatch_score(const Arc &arcA, const Arc &arcB) const;
    
    
    pf_score_t
    boltzmann_weight(score_t s) const { return exp(s/(pf_score_t)params->temperature); }


    //! subtract from each element of a score_t vector v a value x
    void 
    subtract(std::vector<score_t> &v,score_t x) const;
    
    //! subtract from each element of a score_t Matrix m a value x
    void
    subtract(Matrix<score_t> &m,score_t x) const;


public:
    // ------------------------------------------------------------
    // SCORE CONTRIBUTIONS
    
    //! score an match of bases (without structure)
    score_t basematch(size_type i, size_type j) const {
	return sigma_tab(i,j);
    }
    //! score an match of bases with exponential function (without structure)
    pf_score_t exp_basematch(size_type i, size_type j) const {
	return exp_sigma_tab(i,j);
    }
    
    //! score arc match, does *NOT* support explicit arc match scores
    score_t arcmatch(const Arc &arcA, const Arc &arcB, bool stacked=false) const;

    //! score arc match, support explicit arc match scores
    score_t arcmatch(const ArcMatch &am, bool stacked=false) const;

    //! partition function contribution of arc match
    pf_score_t exp_arcmatch(const ArcMatch &am) const {
	return boltzmann_weight(arcmatch(am));
    }
  
    //! stacked arc match
    score_t arcmatch_stacked(const ArcMatch &am) const {
	return arcmatch(am, true);
    }

    //! score of deletion of posA after posB
    //! returns gapA score (normal)
    score_t gapA(size_type posA, size_type posB) const {
	assert(1<=posA && posA <= seqA.length());

	return gapcost_tabA[posA];
    }
    //! returns gapA scores with exponent function
    pf_score_t exp_gapA(size_type posA, size_type posB) const {
	assert(1<=posA && posA <= seqA.length());
	return exp_gapcost_tabA[posA];
    }

    //! score of insertion of posB after posA
    //! returns gapB scores (normal)
    score_t gapB(size_type posA, size_type posB) const {
	assert(1<=posB && posB <= seqB.length());

	return gapcost_tabB[posB];
    }
    
    //! returns gapA scores with exponent function
    pf_score_t exp_gapB(size_type posA, size_type posB) const {
	assert(1<=posB && posB <= seqB.length());

	return exp_gapcost_tabB[posB];
    }
    
    //! cost of an exclusion
    score_t exclusion() const {
	return params->exclusion;
    }   
    
    //! cost to begin a new indel
    score_t indel_opening() const {
	return params->indel_opening;
    }
    
    //! exp of cost to begin a new indel
    pf_score_t exp_indel_opening() const {
	return exp_indel_opening_score;
    }
    
    //
    // ------------------------------------------------------------
    
    //! expected basepair probability,
    //! influences computation of weights for basepairs
    double prob_exp(size_type len) const; //!< return expected probability

    //! returns flag, whether stacking probabilities are used
    bool stacking() const {return params->stacking;}

};

#endif // SCORING_HH
