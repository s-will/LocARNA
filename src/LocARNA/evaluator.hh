#ifndef LOCARNA_EVALUATOR_HH
#define LOCARNA_EVALUATOR_HH

#include <string>
#include "scoring.hh"
#include "multiple_alignment.hh"


namespace LocARNA {

/**
   \brief Evaluate a multiple alignment, given a scoring scheme
      
   Evaluator needs the following  information:
   - scoring parameter
   - pair probabilities or arc match probs/scores
   - whether and how to get base match probs 
   - noLP option
   - whether to evaluate for a fixed structure or optimize
   
   
   It computes the alignment score
   for the given structure or for the best
   consensus structure.

   @todo Implement

*/

class Evaluator {
    MultipleAlignment aln;
    
    const ScoringParams *scoring_params;
    
public:
    //! \brief construct from file
    //! @param alignment_file The name of the alignment file
    //! @param dp_dir            Directory with dot plot information
    //! @param scoring_params The parameters for scoring
    Evaluator(const std::string alignment_file, const std::string dp_dir, const ScoringParams *scoring_params);
    
    //! compute the score of the alignment
    //! @return alignment score
    score_t eval();
};

}

#endif // LOCARNA_EVALUATOR_HH

