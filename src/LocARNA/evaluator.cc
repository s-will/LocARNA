#include "evaluator.hh"

namespace LocARNA {

Evaluator::Evaluator(const std::string alignment_file, 
		     const std::string dir,
		     const ScoringParams *scoring_params_)
    : aln(alignment_file),
      scoring_params(scoring_params_)
{
    
}

score_t Evaluator::eval() {
    return 0;
}

}
