#ifndef EXACT_MATCHER_HH
#define EXACT_MATCHER_HH


#include "scoring.hh"
#include <iostream>
#include <sstream>
#include <list>
#include <algorithm>
#include <iterator>
#include "matrices.hh"
#include <tr1/unordered_map>

extern "C"
{
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/fold.h>
    int    PS_rna_plot(char *string, char *structure, char *file);
    int    PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
    float  fold(const char *sequence, char *structure);
}

namespace LocARNA {

    typedef size_t 					size_type;
    typedef std::vector<unsigned int> 			intVec;
    typedef std::pair<unsigned int,unsigned int> 	intPair;
    typedef std::pair<intPair, intPair> 			intPPair;
    typedef const intPPair* 			intPPairPTR;
    typedef std::vector<intPPair>::const_iterator	IntPPairCITER;


    class SinglePattern
    {
    public:
	SinglePattern(){};
	SinglePattern(const std::string& myId_,const std::string& seqId_,const intVec& mySinglePattern_)
	    :myId(myId_),seqId(seqId_),pattern(mySinglePattern_)
	{};

	virtual ~SinglePattern() { pattern.clear(); };

	const std::string&        getmyId()  const { return myId; };
	const std::string&	getseqId() const {return seqId; };
	const intVec&        getPat() const { return pattern; };

    private:

	std::string         myId;
	std::string	     seqId;
	intVec         pattern;
    };

    //--------------------------------------------------------------------------
    // class PatternPair
    //    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
    //--------------------------------------------------------------------------
    class PatternPair
    {
    public:
	PatternPair(){};
	PatternPair(const std::string& myId,const SinglePattern& myFirstPat,const SinglePattern& mySecPat, const std::string& structure_, int& score_)
	    : id(myId),first(myFirstPat),second(mySecPat), structure(structure_), EPMscore(score_)
	{
	    if (first.getPat().size() != second.getPat().size()){
		std::cerr << "Error! PatternPair cannot be constructed due to different sizes of SinglePatterns!" << std::endl;
	    }
	    score = EPMscore;
	    size = first.getPat().size();
	};

	virtual ~PatternPair()
	{
	    insideBounds.clear();
	};

	const std::string& 		getId() const { return id; };
	const int& 		getSize() const { return size; };
	const SinglePattern& 	getFirstPat() const { return first; };
	const SinglePattern& 	getSecPat() const { return second;};
	void		resetBounds();
	void		setOutsideBounds(intPPair myPPair);
	const   intPPair 		getOutsideBounds() const { return outsideBounds; };
	void		addInsideBounds(intPPair myPPair);
	const   std::vector<intPPair>& getInsideBounds() const { return insideBounds; };

	void			setEPMScore(int myScore);
	const		int 			getScore() const { return score;  };
	const		int 			getEPMScore() const { return EPMscore; };
	std::string& get_struct();

    private:
	std::string         	id;
	int            	size;
	SinglePattern  	first;
	SinglePattern  	second;
	  
	std::string 		structure;
	int				score;
	int				EPMscore;
	std::vector<intPPair>   insideBounds;
	intPPair           outsideBounds;
    };
   
    //--------------------------------------------------------------------------
    // class PatternPairMap
    //    manage a set of EPMs (PatternPair)
    //--------------------------------------------------------------------------
    class PatternPairMap
    {
    public:
	typedef  PatternPair                                  selfValueTYPE;
	typedef  PatternPair*				               		SelfValuePTR;

	typedef  std::multimap<int,SelfValuePTR,std::greater<int> >     orderedMapTYPE;
	typedef  orderedMapTYPE::const_iterator               orderedMapCITER;
	typedef  orderedMapTYPE::iterator                     orderedMapITER;
	typedef  std::list<SelfValuePTR>                           patListTYPE;
	typedef  patListTYPE::iterator                        patListITER;
	typedef  patListTYPE::const_iterator                  patListCITER;
	typedef  std::tr1::unordered_map<std::string,SelfValuePTR> PatternIdMapTYPE;


	PatternPairMap();
	PatternPairMap(const PatternPairMap& myPairMap)
	    :patternList(myPairMap.patternList),
	     patternOrderedMap(myPairMap.patternOrderedMap),
	     idMap(myPairMap.idMap)  { minPatternSize = 100000;};

	virtual ~PatternPairMap();

	void              add( const std::string& id,
			       const SinglePattern& first,
			       const SinglePattern& second,
			       const std::string& structure,
			       int score
			       );
	void              add(const SelfValuePTR value);
	void              makeOrderedMap();
	void              updateFromMap();
	const    PatternPair&      getPatternPair(const std::string& id)const;
	const    SelfValuePTR      getPatternPairPTR(const std::string& id)const;
	const    patListTYPE&      getList() const;
	const    orderedMapTYPE&   getOrderedMap() const;
	orderedMapTYPE&   getOrderedMap2();
	const    int               size()   const;
	int		 getMapBases();
	int  	         getMapEPMScore();
	const    int		 getMinPatternSize() const { return minPatternSize; };

    private:

	patListTYPE        patternList;
	orderedMapTYPE     patternOrderedMap;
	PatternIdMapTYPE   idMap;
	int minPatternSize;
    };


    class LCSEPM
    {
    public:

	LCSEPM(const	Sequence& 		seqA_,
	       const	Sequence& 		seqB_,
	       const 	PatternPairMap& myPatterns,
	       PatternPairMap& myLCSEPM,
	       const	int&		EPM_min_size_ )

	    :seqA(seqA_),
	     seqB(seqB_),
	     matchedEPMs(myLCSEPM),
	     patterns(myPatterns),
	     EPM_min_size(EPM_min_size_){};
		
	virtual		~LCSEPM();

	void 		MapToPS(const std::string& sequenceA, const std::string& sequenceB, PatternPairMap& myMap, const std::string& file1, const std::string& file2);
        void		calculateLCSEPM();

        //! outputs anchor constraints to be used as input for locarna
        void		output_locarna(const std::string& sequenceA, const std::string& sequenceB, const std::string& outfile);
        void		output_clustal(const std::string& outfile_name);

    private:

        struct HoleCompare2 {
	    bool operator()(const intPPairPTR & h1, const intPPairPTR & h2) const {
		// first compare size of holes
		if (h1->first.second - h1->first.first-1 < h2->first.second - h2->first.first-1){
		    return true; }
		// compare if holes are identical in both structures
		if (h1->first.second - h1->first.first-1 == h2->first.second - h2->first.first-1){
		    if ((h1->first.first == h2->first.first) && (h1->first.second == h2->first.second) &&
			(h1->second.first==h2->second.first) && (h1->second.second==h2->second.second))
			{ return true; }
		}

		return false;
	    }
        };

        typedef     std::multimap<intPPairPTR,PatternPairMap::SelfValuePTR,HoleCompare2>	HoleOrderingMapTYPE2;
        typedef     HoleOrderingMapTYPE2::const_iterator HoleMapCITER2;


        void    preProcessing			();
        void    calculateHoles3			();
	void    calculatePatternBoundaries	(PatternPair* myPair);
        void 	calculateTraceback2		(const int i,const int j,const int k,const int l,std::vector < std::vector<int> > holeVec);
        int 	D_rec2				(const int& i,const  int& j,const int& k,const int& l,std::vector < std::vector<int> >& D_h,const bool debug);

        int 	max3				(int a, int b, int c)
	{
	    int tmp = a>b? a:b;
	    return (tmp>c? tmp:c);
	};
	
	//!@brief returns the structure of the given sequence
	char* getStructure(PatternPairMap& myMap, bool firstSeq, int length);
	
	std::string intvec2str(const std::vector<unsigned int>& V, const std::string delim){
	    std::stringstream oss;
	    copy(V.begin(), V.end(), std::ostream_iterator<unsigned int>(oss, delim.c_str()));
	    std::string tmpstr;
	    tmpstr = oss.str();
	    if (tmpstr.length()>0) tmpstr.erase(tmpstr.end()-1);
	    return tmpstr;
	}

	std::string upperCase(const std::string &seq){
	    std::string s= "";
	    for(unsigned int i= 0; i<seq.length(); i++)
		s+= toupper(seq[i]);
	    return s;
	}

	std::vector< std::vector <std::vector<PatternPairMap::SelfValuePTR> > >	EPM_Table2;
	HoleOrderingMapTYPE2    		holeOrdering2;
	const 	Sequence&				seqA;
	const 	Sequence& 				seqB;
	PatternPairMap&				matchedEPMs;
        const 	PatternPairMap&         		patterns;
	const int& 					EPM_min_size;
    };


    class Mapping{
	typedef std::vector<int> pos_vec;
	typedef std::vector<pos_vec> bp_mapping;
	

    public:
	//! constructor
	Mapping(const BasePairs &bps_,const RnaData &rnadata_,
		const double &prob_unpaired_in_loop_threshold_,
		const double &prob_unpaired_in_F_threshold_):
	    prob_unpaired_in_loop_threshold(prob_unpaired_in_loop_threshold_),
	    prob_unpaired_in_F_threshold(prob_unpaired_in_F_threshold_),
	    bps(bps_),
	    rnadata(rnadata_)
	{
	    compute_mapping();
	}


    private:
	
	const double &prob_unpaired_in_loop_threshold;
	const double &prob_unpaired_in_F_threshold;
	const BasePairs &bps;
	const RnaData &rnadata;
	bp_mapping pos_vecs; //! mapping from the new positions to the sequence positions (i.e. which positions relative to the beginning of the arc are valid)
	bp_mapping new_pos_vecs; //!mapping from the sequence positions to the new positions;
		                 //!sequence positions are relative to the beginning of the arc
                                 //!vec contains -1 if sequence position isn't valid
		
	void compute_mapping();

    public:
  
	//!is sequential matching from position new_pos-1 to position new_pos possible?
	bool seq_matching(size_type arcIdx,size_type new_pos)const {
	    return pos_vecs.at(arcIdx).at(new_pos-1)+1==pos_vecs.at(arcIdx).at(new_pos);
	}
	
	//!is position k valid (i.e. does probability that the base k is unpaired under the loop exceed some threshold) for the basepair with index arcIdx? 
	bool is_valid_pos(Arc arc,size_type k) const{
	    return rnadata.prob_unpaired_in_loop(k,arc.left(),arc.right())>=prob_unpaired_in_loop_threshold;
	}

	//!returns the sequence position corresponding to the position new_pos in the matrix
	int get_pos_in_seq(const Arc &arc, size_type new_pos) const{
	    return pos_vecs.at(arc.idx()).at(new_pos)+arc.left();
	}

	//!returns the new position in the matrix corresponding to the position pos in the sequence
	//!returns -1 if position pos isn't valid
	int get_pos_in_new_seq(const Arc &arc, size_type pos) const{
	    return new_pos_vecs.at(arc.idx()).at(pos-arc.left());
	}

	//!returns the number of valid positions for a basepair with index arcIdx
	int number_of_valid_pos(size_type arcIdx) const{
	    return pos_vecs.at(arcIdx).size();
	}

	bool unpaired_in_F(size_type k) const{
	    return rnadata.prob_unpaired(k)>=prob_unpaired_in_F_threshold;
	}

	//for debugging
	void print_vec() const;

	//!class distructor
	virtual ~Mapping(){
	    new_pos_vecs.clear();
	    pos_vecs.clear();
	}
	   
    };

    //!a class for the representation of exact pattern matches (EPM)
    class EPM{

	intVec pat1Vec;
	intVec pat2Vec;
	std::string structure;
	int score; //needed for suboptimal trace
	int state; //needed for suboptimal trace in AGB
        std::pair<int,int> curPos; //needed for suboptimal trace
	int score_to_substract_in_F; //needed for suboptimal trace

	//!contains the indices of the arcMatches which need to be considered
	std::vector<ArcMatch::idx_type> arcmatches_to_do;

	void swap(int i, int j){
	    swap(i,j,pat1Vec);
	    swap(i,j,pat2Vec);
	    //swap position i and j in structure
	    char tmp = structure[i];
	    structure[i]=structure[j];
	    structure[j]=tmp;
	}

	void swap(int i, int j, std::vector<unsigned int> &vec){
	    int tmp = vec[i];
	    vec[i]=vec[j];
	    vec[j]=tmp;
	}

	int split(int left,int right){
	    unsigned int pivot=pat1Vec[right];
	    int i=left;
	    int j=right-1;
	    while(i<j){
		while(pat1Vec[i]<=pivot && i<right){i++;}
		while(pat1Vec[j]>=pivot && j>left){j--;}
		if(i<j){swap(i,j);}
	    }
	    if(i<right){
		if(pat1Vec[i] > pivot){swap(i,right);}
	    }
	    return i;
	}

	void quicksort_EPM(int left,int right){
	    if(left<right) {
		int Pivot_idx=split(left,right);
		quicksort_EPM(left,Pivot_idx-1);
		quicksort_EPM(Pivot_idx+1,right);
	    }
	}


    public:
	
	//!Constructor
	EPM(){
	    reset();
	}

	virtual ~EPM(){
	    pat1Vec.clear();
	    pat2Vec.clear();
	    arcmatches_to_do.clear();
	}

	void set_arcmatches_to_do(std::vector<ArcMatch::idx_type> arcmatches_to_do_){
	    arcmatches_to_do = arcmatches_to_do_;
	}

	void set_score(int score_){
	    score=score_;
	}

	void set_state(int state_){
	    state=state_;
	}

	void set_curPos(std::pair<int,int> curPos_){
	    curPos = curPos_;
	}

	void set_score_to_substract_in_F(int score){
	    score_to_substract_in_F=score;
	}

	int get_score(){
	    return score;
	}

	int get_state(){
	    return state;
	}

        std::pair<int,int> get_curPos(){
	    return curPos;
	}

	int get_score_to_substract_in_F(){
	    return score_to_substract_in_F;
	}

	std::vector<ArcMatch::idx_type> get_arcmatches_to_do(){
	    return arcmatches_to_do;
	}

	intVec getPat1Vec() const{
	    return pat1Vec;
	}

	intVec getPat2Vec() const{
	    return pat2Vec;
	}

	std::string getStructure() const{
	    return structure;
	}

	//!reset epm
	void reset(){
	    pat1Vec.clear();
	    pat2Vec.clear();
	    structure.clear();
	    arcmatches_to_do.clear();
	    score = 0;
	}

	void add(int pos1_, int pos2_,char c){
	    pat1Vec.push_back(pos1_);
	    pat2Vec.push_back(pos2_);
	    structure.push_back(c);
	}

	//!appends an arcMatch to the epm
	void add_arcmatch(const ArcMatch &am){
	    add(am.arcA().left(),am.arcB().left(),'(');
	    add(am.arcA().right(),am.arcB().right(),')');
	}

	void store_arcmatch(ArcMatch::idx_type idx){
	    arcmatches_to_do.push_back(idx);
	}

	//!checks if there are arcMatches left which need to be processed
	bool arcmatch_to_process(){
	    return arcmatches_to_do.begin()!=arcmatches_to_do.end();
	}

	//!returns the index of the last arcMatch in the vector arcmatches_to_do
	ArcMatch::idx_type next_arcmatch(){
	    ArcMatch::idx_type arc_idx = arcmatches_to_do.back();
	    arcmatches_to_do.pop_back();
	    return arc_idx;
	}

	void sort_patVec(){
	    quicksort_EPM(0,pat1Vec.size()-1);
	}

	bool isEmpty(){
	    //return epm.empty();
	    return pat1Vec.empty();
	}

	void print_epm(std::ostream &out, int score){
	    out << "epm with score " << score << std::endl;
	    intVec::iterator it2=pat2Vec.begin();
	    out << " ";
	    for(intVec::iterator it=pat1Vec.begin();it!=pat1Vec.end();it++,it2++){
		out << *it;
		out << ":";
		out << *it2 << " " ;
	    }
	    out << std::endl;
	    out << " ";
	    for(std::string::iterator it=structure.begin();it!=structure.end();it++){
		out << *it;
	    }
	    out << std::endl;
	    out << "pos " << this->curPos.first << "," << this->curPos.second << std::endl;
	    out << "state " << this->state << std::endl;
	}
    };

    class ExactMatcher {
  
  
    private:

	typedef size_t size_type;
	const Sequence &seqA;
	const Sequence &seqB;
	const ArcMatches &arc_matches;
	const BasePairs &bpsA;
	const BasePairs &bpsB;
	const Mapping &mappingA;
	const Mapping &mappingB;
    

	EPM epm;

	ScoreMatrix A;
	ScoreMatrix G;
	ScoreMatrix G2; //for suboptimal traceback -> unambig
	ScoreMatrix B;
	ScoreMatrix F;

	ScoreVector arc_match_score; //!vector for the arcMatch scores: score under the arcMatch with potential stacking probabilities

	struct Trace_entry{
	    infty_score_t score;
	    std::pair<int,int> *next_pos;
	    ArcMatch::idx_type *arc_match_idx;
	};


	Matrix<Trace_entry> Trace; //!for traceback for heuristic case

	std::pair<int,int> pos_of_max;
	int alpha_1;
	int alpha_2;
	int alpha_3;
	int difference_to_opt_score;
	int min_subopt_score;
	int easier_scoring_par;
	double subopt_range;
	score_t am_threshold;
	int subopt_score;
	double cutoff_coverage;

	enum{in_B,in_G,in_G2,in_A};

	PatternPairMap& foundEPMs;

	struct info_for_trace_in_G_subopt{
	    std::pair<int,int> curPos;
	    int state;
	    int score;
	};

	// ----------------------------------------
	// evaluate the recursions / fill matrices
    

	//!computes matrices A,G and B
	void compute_AGBmatrices(const ArcMatch &arc_match);
    
	//!helper function for compute_matrices
	infty_score_t seq_str_matching(ScoreMatrix &mat,const ArcMatch &arc_match, size_type i, size_type j,bool matrixB, bool subopt = false);

	//!computes matrix F with filter criteria
	void compute_F_with_prob_unpaired();
    
	//! computes score for arcMatch: score under the arcMatch plus the probability of the two basepairs of the arcMatch
	infty_score_t score_for_arc_match(const ArcMatch &am, bool with_part_under_am);

	//! computes the stacking score: if stacking occurs with respect to a structure, the stacking probability is taken as a score
	infty_score_t score_for_stacking(const ArcMatch &am, const ArcMatch &inner_am);

	//! add current epm to set of all EPMs (PatternPairMap)
	void add_foundEPM();

	//!computes the length of the EPM with the highest score -> calculate coverage
	int compute_length_of_best_EPM();

	//!starts suboptimal traceback from all positions where the EPM cannot be extended
	void find_start_pos_for_traceback();

	//!calculates suboptimal tracebacks
	void trace_in_F_suboptimal(int i, int j);

	//!calculates suboptimal tracebacks in A,G and B
	void trace_AGB_suboptimal(const ArcMatch &am, EPM &epm_to_store,EPM cur_epm_AM, std::list<EPM> &epms_to_proc_AM);

	//!calculates matrices A,G,G2 and B for nonambiguous suboptimal traceback
	void compute_AGB_subopt(const ArcMatch &arc_match);

	//!helper function for trace_AGB_suboptimal
	bool str_trace_AGB_suboptimal(ScoreMatrix &mat, size_type posA, size_type posB,const ArcMatch &am,int cur_score,int state, std::pair<int,int> curPos,std::list<EPM> &epms_to_proc_AM,EPM cur_epm_AM);

	//!calculates suboptimal traceback in G
	void trace_suboptimal_in_G(int state,int cur_score,std::pair<int,int> curPos,std::list<ExactMatcher::info_for_trace_in_G_subopt> &result);

	//!puts together the EPM from the part before the am and the part under the am
	void copy_epm_trace(std::list<EPM> &epms_to_proc_AM,EPM &trace_F, std::list<EPM> &epms_to_proc);

	//!returns if a position (i,j) is likely to be unpaired in F
	bool valid_unpaired_pos_in_F(size_type i, size_type j);

	//!helper function for trace_suboptimal_in_G
	void store_cur_el(std::pair<int,int> pos, int score, int state,std::list<ExactMatcher::info_for_trace_in_G_subopt> &result);

	//!compute the backward score and the forward pointer
	void compute_trace_matrix();

	//!traverses matrix F
	void trace_in_F();
    
	//!outputs the exact matching starting at position $(i,j)$ while setting the processed elements Trace(i,j) to -inf
	void get_matching(size_type i, size_type j);

	void set_el_to_neg_inf();

	//!checks the structural case for the traceback in the matrices A,G and B
	//!if an inner_am is encountered, it is stored for later processing and the left and right endpoint is added to the epm structure
	bool str_traceAGB(const ScoreMatrix &mat, const ArcMatch &am, size_type posA, size_type posB,std::pair<int,int> &curPos, EPM &epm_to_store);

	//!recomputes matrices A,G and B for arcMatch recursively and stores the traceback
	bool trace_AGB(const ArcMatch &am, EPM &epm_to_store);

	//!adds the arcMatch to the epm if the positions aren't contained in another epm (score==pos_infty)
	bool add_arcmatch(const ArcMatch &am, EPM &epm_to_store);

	//!adds the structure and the position corresponding to the position pos_ in the matrix if
	//!the position isn't contained in another epm (score==pos_infty)
	bool add(const ArcMatch &am,std::pair<int,int> pos_, char c, EPM &epm_to_store);
    
	//Debugging
	void print_EPM_start_pos(std::list<std::pair<std::pair<int,int>,infty_score_t> > &EPM_start_pos);
	void output_trace_matrix();
	void output_arc_match_score();
	void print_matrices(const ArcMatch &am, size_type offset_A, size_type offset_B);
	void print_epms_to_proc_AGB(std::list<EPM> &epms_to_proc_AGB);
	void print_info_for_G_matrix(std::list<info_for_trace_in_G_subopt> &result);
	void print_arcmatches_to_do(std::vector<ArcMatch::idx_type> arcmatches_to_do);
	bool validate_epm();

	//!converts std::string to uppercase
	std::string upperCase(std::string seq){
	    std::string s= "";
	    for(unsigned int i= 0; i<seq.length(); i++)
		s+= toupper(seq[i]);
	    return s;
	}
    
    public:

	//! construct with sequences and possible arc matches
	ExactMatcher(const Sequence &seqA_,
		     const Sequence &seqB_,
		     const ArcMatches &arc_matches_,
		     const Mapping &mappingA_,
		     const Mapping &mappingB_,
		     PatternPairMap &foundEPMs_,
		     const int &alpha_1,
		     const int &alpha_2,
		     const int &alpha_3,
		     const int &difference_to_opt_score,
		     const int &min_subopt_score,
		     const int &easier_scoring_par,
		     const double &subopt_range_,
		     const int &am_threshold_,
		     const double &cutoff_coverage_
		     );
	~ExactMatcher();
    
	//! fills the A,G,B and F matrices
	void
	compute_matrices();
    
	//! store all exact matchings in PatternPairMap and call chaining algorithm
	//! pre: call to compute_marices()
	void
	compute_EPMs_heuristic();
    
	void
	compute_EPMs_suboptimal();
    };



} //end namespace

#endif //  EXACT_MATCHER_HH
