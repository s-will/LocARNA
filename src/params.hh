#ifndef PARAMS_HH
#define PARAMS_HH
/*
  Parameter for Aligner
*/

class AnchorConstraints;
class TraceController;


/**
   Description of free end gaps.
   Decodes the description given by a string of 4 characters '+'/'-'
   and provides methods with reasonable names.
 */
class FreeEndgapsDescription {
    std::vector<bool> desc;
public:
    FreeEndgapsDescription(const std::string d) {
       	desc.resize(4);
	if (d.length()>=4) {
	    for (size_t i=0; i<4; i++) desc[i] = (d[i]=='+');
	} else {
	    for (size_t i=0; i<4; i++) desc[i] = false;
	}
    }
    
    bool
    allow_left_1() const  {
	return desc[0];
    }
    bool
    allow_right_1() const  {
	return desc[1];
    }
    bool
    allow_left_2() const  {
	return desc[2];
    }
    bool
    allow_right_2() const  {
	return desc[3];
    }
};


/**
   Collects the parameters for the aligner object.  These parameters
   controll the kind of alignment (local/global),
   restrictions/constraints on the alignment and certain heuristics.
   Parameters for the score are collected in a different class.
*/
class AlignerParams {  
public:
    typedef size_t size_type;

    const bool no_lonely_pairs; // no lonely pairs option
  
    const bool STRUCT_LOCAL; // allow exclusions for maximizing alignment of connected substructures 
    const bool SEQU_LOCAL; // maximize alignment of subsequences

    const FreeEndgapsDescription free_endgaps;
    
    const bool DO_TRACE;

    const TraceController &trace_controller;

    const int max_diff_am;
    const double min_am_prob;
    const double min_bm_prob; 
    
    const int stacking;

    const AnchorConstraints &constraints;
    
    
    AlignerParams(bool _no_lonely_pairs, 
		  bool _STRUCT_LOCAL, 
		  bool _SEQU_LOCAL, 
		  std::string _free_endgaps,
		  const TraceController &_trace_controller,
		  int _max_diff_am,
		  const double _min_am_prob,
		  const double _min_bm_prob,
		  bool _stacking,
		  const AnchorConstraints &_constraints		  
	):
	no_lonely_pairs(_no_lonely_pairs),
	STRUCT_LOCAL(_STRUCT_LOCAL),
	SEQU_LOCAL(_SEQU_LOCAL),
	free_endgaps(_free_endgaps),
	DO_TRACE(true),
	trace_controller(_trace_controller),
	max_diff_am(_max_diff_am),
	min_am_prob(_min_am_prob), 
	min_bm_prob(_min_bm_prob),	   
	stacking(_stacking),
	constraints(_constraints) 
	{}
  
};

#endif // PARAMS_HH
