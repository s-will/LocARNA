#ifndef EXACT_MATCHER_HH
#define EXACT_MATCHER_HH
#include "scoring.hh"
#include <iostream>
#include <sstream>
#include <list>
#include <algorithm>
#include "matrices.hh"

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

using namespace std;

namespace LocARNA {

typedef size_t size_type;
typedef vector<unsigned int> intVec;
typedef pair<unsigned int,unsigned int> intPair;
typedef pair<intPair, intPair> intPPair;
typedef vector<intPPair>::const_iterator IntPPairCITER;

class StringHash
{
   public:
      size_t operator()(const string &myStr) const
      {
         unsigned long hash = 5381;

         for (unsigned int i = 0; i < myStr.length(); i++)
         {
            hash = ((hash << 5) + hash) + myStr[i]; // hash * 33 + cStateString[i]
         }
         return hash;
   }
};

class StringEq {
public:
   bool operator()(const string &a,const string &b) const
   {
      return (a == b);
   }
};

class SinglePattern
{
public:
      SinglePattern(){};
      SinglePattern(const string& myId,const string& seqId,const intVec& mySinglePattern);

	virtual ~SinglePattern();

   const string&        getmyId()  const;
   const string&	getseqId() const;
   const intVec&        getPat() const;

private:

      string         myId;
      string	     seqId;
      intVec         pattern;
};

//--------------------------------------------------------------------------
// class PatternPair
//    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
//--------------------------------------------------------------------------
class PatternPair
   {
      public:
      PatternPair(const string& myId,const int& mySize,const SinglePattern& myFirstPat,const SinglePattern& mySecPat, const string& structure_, int& score_)
                  : id(myId),size(mySize),first(myFirstPat),second(mySecPat), structure(structure_), EPMscore(score_)
      {  score = EPMscore;  };

      virtual ~PatternPair()
      {
    	insideBounds.clear();
      };

      const string& 		getId() const;
      const int& 			getSize() const;
      const SinglePattern& 	getFirstPat() const;
      const	SinglePattern& 	getSecPat() const;
			void			resetBounds();
			void			setOutsideBounds(intPPair myPPair);
			intPPair 		getOutsideBounds();
			void			addInsideBounds(intPPair myPPair);
      const vector<intPPair>& getInsideBounds();
			void			setEPMScore(int myScore);
			int 			getScore();
			int 			getEPMScore();
      string& get_struct();

      private:
         string         	id;
         int            	size;
         SinglePattern  	first;
         SinglePattern  	second;
	  
	 string 		structure;
         int				score;
         int				EPMscore;
         vector<intPPair>   insideBounds;
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

      typedef  multimap<int,SelfValuePTR,greater<int> >     orderedMapTYPE;
      typedef  orderedMapTYPE::const_iterator               orderedMapCITER;
      typedef  orderedMapTYPE::iterator                     orderedMapITER;
      typedef  list<SelfValuePTR>                           patListTYPE;
      typedef  patListTYPE::iterator                        patListITER;
      typedef  patListTYPE::const_iterator                  patListCITER;
      typedef  __gnu_cxx::hash_map<string,SelfValuePTR,StringHash,StringEq> PatternIdMapTYPE;


         PatternPairMap();
         PatternPairMap(const PatternPairMap& myPairMap)
                           :patternList(myPairMap.patternList),
                            patternOrderedMap(myPairMap.patternOrderedMap),
                            idMap(myPairMap.idMap)  {};

      virtual ~PatternPairMap();

               void              add( const string& id,
                                      const int& mysize,
                                      const SinglePattern& first,
                                      const SinglePattern& second,
				      const string& structure,
				      int score
				    );
               void              add(const SelfValuePTR value);
               void              makeOrderedMap();
               void              updateFromMap();
      const    PatternPair&      getPatternPair(const string& id)const;
      const    SelfValuePTR      getPatternPairPTR(const string& id)const;
      const    patListTYPE&      getList() const;
      const    orderedMapTYPE&   getOrderedMap() const;
               orderedMapTYPE&   getOrderedMap2();
      const    int               size()   const;
			   int  			 getMapBases();

   private:

     patListTYPE        patternList;
     orderedMapTYPE     patternOrderedMap;
     PatternIdMapTYPE   idMap;
};

class LCSEPM
{
    public:

					LCSEPM(const int& size1_, const int& size2_,
					  const 		PatternPairMap& myPatterns,
										PatternPairMap& myLCSEPM,
										int&		    mySize,
										int&			myScore,
						const int& EPM_min_size_
					      )

    	                                   :size1(size1_),
					    size2(size2_),
					    matchedEPMs(myLCSEPM),
    	                                    patterns(myPatterns),
    	                                    LCSEPMsize(mySize),
    	                                    LCSEPMscore(myScore),
    	                                    EPM_min_size(EPM_min_size_){};
		
					void MapToPS(const string& sequenceA, const string& sequenceB, int& mySize, PatternPairMap& myMap, const string& file1, const string& file2);
        virtual 		~LCSEPM();

        		void    calculateLCSEPM();

	private:

		struct	HoleKeyS;
        struct	HoleKeyS
        {
        	PatternPairMap::SelfValuePTR    pattern;
            intPPair                        bounds;
        };

        typedef     HoleKeyS                            HoleKey;
        typedef     HoleKey*                            HoleKeyPTR;
        typedef     multimap<int,HoleKeyPTR>            HoleOrderingMapTYPE;
        typedef     HoleOrderingMapTYPE::const_iterator HoleMapCITER;

        void    preProcessing				();
        void    calculateHoles				();
	    void    calculatePatternBoundaries	(PatternPair* myPair);
        void 	calculateTraceback			(const int i,const int j,const int k,const int l,vector < vector<int> > holeVec);
        int 	D_rec						(const int& i,const  int& j,const int& k,const int& l,vector < vector<int> >& D_h,const bool debug);
        int 	max3						(int a, int b, int c);
	
	//!@brief returns the structure of the given sequence
	char* getStructure(PatternPairMap& myMap, bool firstSeq, int length);
	
	string intvec2str(const std::vector<unsigned int>& V, const std::string delim){
	      stringstream oss;
	      copy(V.begin(), V.end(), ostream_iterator<unsigned int>(oss, delim.c_str()));
	      string tmpstr;
	      tmpstr = oss.str();
	      if (tmpstr.length()>0) tmpstr.erase(tmpstr.end()-1);
	      return tmpstr;
	  }
	string upperCase(string seq){
	  string s= "";
	  for(unsigned int i= 0; i<seq.length(); i++)
	    s+= toupper(seq[i]);
	  return s;
	}
				vector< vector<PatternPairMap::SelfValuePTR> >	EPM_Table;
				const int& size1;
				const int& size2;
				PatternPairMap&									matchedEPMs;
				HoleOrderingMapTYPE    					 		holeOrdering;
        const 	PatternPairMap&         						patterns;
				int& 											LCSEPMsize;
				int& 											LCSEPMscore;
	const int& EPM_min_size;
				
};


class Mapping{
	typedef std::vector<int> pos_vec;
	typedef std::vector<pos_vec> bp_mapping;
	

public:
	  //! constructor
	Mapping(const BasePairs &bps_,const RnaData &rnadata_,
			const double &prob_unpaired_in_loop_threshold_,
			const double &prob_unpaired_external_threshold_,
			const double &prob_basepair_external_threshold_):
	  prob_unpaired_in_loop_threshold(prob_unpaired_in_loop_threshold_),
	  prob_unpaired_external_threshold(prob_unpaired_external_threshold_),
	  prob_basepair_external_threshold(prob_basepair_external_threshold_),
	  bps(bps_),
	  rnadata(rnadata_)
	  {
		compute_mapping();
	}


private:
	
	const double &prob_unpaired_in_loop_threshold;
	const double &prob_unpaired_external_threshold;
	const double &prob_basepair_external_threshold;
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
		//return true; //TODO remove!
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
	
	bool unpaired_external(size_type k) const {
		//cout << "prob unpaired external " << rnadata.prob_unpaired_external(k) << endl;
		//cout << "prob unapired external threshold " << prob_unpaired_external_threshold;
		return rnadata.prob_unpaired_external(k)>=prob_unpaired_external_threshold;
	}

	bool basepair_external(size_type i, size_type j) const{
		cout << "i,j " << i << "," << j << endl;
		cout << "prob basepair external " <<  rnadata.prob_basepair_external(i,j) << endl;
		return rnadata.prob_basepair_external(i,j)>=prob_basepair_external_threshold;
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
	 //const Sequence &seqA_;
	 //const Sequence &seqB_;

	/*struct EPM_entry{
		int pos1;
		int pos2;
		char str;
	};*/

	intVec pat1Vec;
	intVec pat2Vec;
	string structure;

	//! the EPM consists of the positions and the corresponding structure
	//list<EPM_entry> epm;

	//typedef list<EPM_entry>::iterator iter;
	//iter cur_it; //! current position

        //!contains the indices of the arcMatches which need to be considered and a corresponding iterator which
	//!gives the position where the traceback of the arcMatch should be inserted
	//std::vector<pair<ArcMatch::idx_type,iter> > arcmatches_to_do;
	std::vector<ArcMatch::idx_type> arcmatches_to_do;

	void swap(int i, int j){
		swap(i,j,pat1Vec);
		swap(i,j,pat2Vec);
		//swap position i and j in structure
		char tmp = structure[i];
		structure[i]=structure[j];
		structure[j]=tmp;
	}

	void swap(int i, int j, vector<unsigned int> &vec){
		int tmp = vec[i];
		vec[i]=vec[j];
		vec[j]=tmp;
	}

	int split(int left,int right){
		   int pivot=pat1Vec[right];
		   int i=left;
		   int j=right-1;
		   while(i<j){
		   while(pat1Vec[i]<=pivot && i<right){
			   i++;
		   }
		   while(pat1Vec[j]>=pivot && j>left){
			   j--;
		   }
		   if(i<j){
			   swap(i,j);
		   }
		   }
		   if(i<right){
			   if(pat1Vec[i] > pivot){
				   swap(i,right);
			   }
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
	//EPM(const Sequence &seqA, const Sequence &seqB)
	//: seqA_(seqA),seqB_(seqB){//, mcsPatterns() {
		EPM(){
		reset();
	}

	virtual ~EPM(){
		pat1Vec.clear();
		pat2Vec.clear();
		arcmatches_to_do.clear();
	}
	//!reset epm and reset current position (cur_it) to the beginning of epm
	void reset(){
		pat1Vec.clear();
		pat2Vec.clear();
		//epm.clear();
		structure.clear();
		//cur_it=epm.begin();
	}

	//!reset current position to the end position of epm
	/*void reset_cur_it(){
		cur_it=epm.end();
		cur_it--;
	}*/

	void add(int pos1_, int pos2_,char c){
		pat1Vec.push_back(pos1_);
		pat2Vec.push_back(pos2_);
		structure.push_back(c);
		//cout << "add 1 " << pos1_ << endl;
		//cout << "add 2 " << pos2_ << endl;
		/*cout << "pat1Vec " << endl;
		for(int i=0;i<pat1Vec.size(); i++){
			cout << pat1Vec.at(i) << " ";
		}
		cout << endl;
		cout << "pat2Vec " << endl;
		for(int i=0;i<pat2Vec.size(); i++){
			cout << pat2Vec.at(i) << " ";
		}
	    cout << endl;*/
	}

	//!appends an arcMatch to the epm
	void add_arcmatch(const ArcMatch &am){
		//cout << "add_test( " << am.arcA().left() << "," << am.arcB().left() << "'('" << endl;
		//cout << "add_test( " << am.arcA().right() << "," << am.arcB().right() << "'('" << endl;
		//cout << "add arcmatch " << am.arcA().left() << "," << am.arcB().left()<< endl;
		//cout << "add arcmatch " << am.arcA().right() << "," << am.arcB().right()<< endl;
		add(am.arcA().left(),am.arcB().left(),'(');
		add(am.arcA().right(),am.arcB().right(),')');
	}

	//!appends a position and the corresponding structure to the epm and increments the iterator
	/*void append(int pos1_, int pos2_, char c){
		struct EPM_entry entry = {pos1_,pos2_,c};
		epm.push_back(entry);
		//pat1Vec.push_back(pos1_);
		//pat2Vec.push_back(pos2_);
		//structure.push_back(c);
		cur_it++;
	}*/

	//!appends an arcMatch to the epm
	/*void append_arcmatch(const ArcMatch &am){
		append(am.arcA().left(),am.arcB().left(),'(');
		append(am.arcA().right(),am.arcB().right(),')');
	}*/

	//!adds a positions and the corresponding structure to the epm at cur_it and decrements the iterator
	/*void add(int pos1_, int pos2_,char c){
		struct EPM_entry entry = {pos1_,pos2_,c};
		epm.insert(cur_it,entry);
		cur_it--;
		//pat1Vec.push_back(pos1_);
		//pat2Vec.push_back(pos2_);
		//structure.push_back(c);
	}*/


	//!adds an arcmatch at the current position (cur_it) and stores the corresponding arcMatch
	/*void add_arcmatch(const ArcMatch &am){
		add(am.arcA().right(),am.arcB().right(),')');
		store_arcmatch(am.idx());
		add(am.arcA().left(),am.arcB().left(),'(');
	}*/

	//!stores the index of the arcMatch and the current position (cur_it) in the vector arc_matches_to_do
	/*void store_arcmatch(ArcMatch::idx_type idx){
		//list_it++;
		//cout << "store archmatch " << endl;
		//cout << "pos for arcmatch " << list_it << endl;
		arcmatches_to_do.push_back(pair<ArcMatch::idx_type,iter>(idx,cur_it));
		//arcmatches_to_do_test.push_back(pair<ArcMatch::idx_type,int>(idx,list_it));
		//list_it--;
	}*/

	void store_arcmatch(ArcMatch::idx_type idx){
		arcmatches_to_do.push_back(idx);
	}

	//!checks if there are arcMatches left which need to be processed
	bool arcmatch_to_process(){
		return arcmatches_to_do.begin()!=arcmatches_to_do.end();
	}

	//!returns the index of the last arcMatch in the vector arcmatches_to_do
	ArcMatch::idx_type next_arcmatch(){

		//pair<ArcMatch::idx_type,iter> lastEl = arcmatches_to_do.back();

		//cur_it = lastEl.second;

		//ArcMatch::idx_type result=lastEl.first;
		//arcmatches_to_do.pop_back();
		ArcMatch::idx_type arc_idx = arcmatches_to_do.back();
		arcmatches_to_do.pop_back();
		//arcmatches_to_do_test.pop_back();
		return arc_idx;
	}

	void validate_epm(){
		//EPMs have the same size
		if(pat1Vec.size()!=pat2Vec.size()){
			cerr << "wrong pattern " << endl;
			return;
		}

		//two matched positions in the EPMs have the same nucleotide
		/*for(int i=0;i<pat1Vec.size();i++){
			if(seqA_[pat1Vec.at(i)]!=seqB_[pat2Vec.at(i)]){
				cerr << "no EPM " << endl;
			}
	    }*/

		vector<pair<int,int> > arcmatches_to_validate;
		int balance=0;

		for(int i=0;i<pat1Vec.size()-1;i++){
			if(pat1Vec[i]>=pat1Vec[i+1] || pat2Vec[i]>=pat2Vec[i+1]){
				cerr << "ERROR" << endl;
				cout << "hier " << endl;
				print_epm(cout,0);
			}
		}
		bool gap = false;
		arcmatches_to_validate.push_back(pair<int,int>(0,pat1Vec.size()-1));
		//for(int i=0;i<pat1Vec.size()-1;i++){
		while(arcmatches_to_validate.size()!=0){
			pair<int,int> cur_arcmatch = arcmatches_to_validate.back();
			arcmatches_to_validate.pop_back();
			gap=false;
			//cout << "validate cur arcmatch " << cur_arcmatch.first << "," << cur_arcmatch.second << endl;
			for(int i=cur_arcmatch.first+1;i<=cur_arcmatch.second;i++){
			//cout << "i " << i << endl;
			if(structure.at(i)=='.'){
				//cout << "validated " << i-1 << "," << i << endl;
				if(i!=0 && !(pat1Vec.at(i-1)+1==pat1Vec.at(i))){
					//cout << "gap " << endl;
					if(gap) cerr << "ERROR";//cout << "from position " << i-1 << " to " << i << endl; //cerr << "ERROR";
					else gap=true;
				}
			}
			else if(structure.at(i)=='('){
				int start = i;
				//cout << "validated " << i-1 << "," << i << endl;
				if(!(pat1Vec.at(i-1)+1==pat1Vec.at(i))){
					//cout << "gap " << endl;
					if(gap) cerr << "ERROR " << endl;
					gap=true;
				}
				//gap=false;
				//i++;
				//cout << "i " << i << endl;
				//cout << "pat vec " << structure.at(i) << endl;
				int balance = 1; //count first open bracket at position i
				//i--;
				while(balance!=0){
					i++;
				//	cout << "i " << i << endl;
					//cout << "structure " << structure.at(i) << endl;
					if(structure.at(i)=='(') balance++;
					if(structure.at(i)==')') balance--;
				//	cout << "balance " << balance << endl;
				//	cout << "i " << i << endl;
				}
				//cout << "arcmatches_to_validate " << start <<","<< i << endl;
				arcmatches_to_validate.push_back(pair<int,int>(start,i));

			}
			else if(structure.at(i)==')'){
				//cout << "validated " << i-1 << "," << i << endl;
				if(!(pat1Vec.at(i-1)+1==pat1Vec.at(i))){
					//cout << "gap " << endl;
					if(gap) cerr << "ERROR " << endl;
					gap=true;
				}
			}
			}
		}
	}

	void print_epm(ostream &out, int score){
		cout << "epm with score " << score << endl;
		intVec::iterator it2=pat2Vec.begin();
		for(intVec::iterator it=pat1Vec.begin();it!=pat1Vec.end();it++,it2++){
			out << *it;
			out << ":";
			out << *it2 << " " ;
		}
		out << endl;
		for(string::iterator it=structure.begin();it!=structure.end();it++){
					out << *it;
				}
		cout << endl;

	  }


	//for debugging
	/*void print_arcmatches_to_do(ostream &out){
		out << "to do vector " << endl;
		for(std::vector<pair<ArcMatch::idx_type,iter> >::iterator it=arcmatches_to_do.begin();it!=arcmatches_to_do.end();it++){
			out << it->first << endl;
		}
		out << endl;
	}*/
	void sort_patVec(){
		//sort(pat1Vec.begin(), pat1Vec.end());
		quicksort_EPM(0,pat1Vec.size()-1);
		//sort(pat2Vec.begin(), pat2Vec.end());
	}
	intVec getPat1Vec() const{
		return pat1Vec;
	}
	intVec getPat2Vec() const{
		return pat2Vec;
	}
	
	string getStructure() const{
		return structure;
	}

	/*void add_pattern(string patId ,SinglePattern pattern1 ,SinglePattern pattern2, int score){
	  //set_struct();
	  mcsPatterns.add(patId, pat1Vec.size(), pattern1, pattern2, structure, score);
	}*/
	
	/*PatternPairMap& get_patternPairMap() {
	  return mcsPatterns;
	}*/

	bool isEmpty(){
	  //return epm.empty();
		return pat1Vec.empty();
	}

	/*void set_struct(){
	  structure= "";
	  for(iter it=epm.begin();it!=epm.end();it++){
			structure+= it->str;
		}
	}*/
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
    ScoreMatrix B;
    ScoreMatrix F;

    ScoreVector arc_match_score; //!vector for the arcMatch scores: score under the arcMatch with potential stacking probabilities

    struct Trace_entry{
    	infty_score_t score;
    	pair<int,int> *next_pos;
    	ArcMatch::idx_type *arc_match_idx;
    };


    Matrix<Trace_entry> Trace; //!for traceback

    int EPM_threshold;
    int EPM_min_size;
    int alpha_1;
    int alpha_2;
    int alpha_3;
    int easier_scoring_par;
    int subopt_score;
    const string& sequenceA;
    const string& sequenceB;
    enum{in_B,in_G,in_A};
    const string& file1;
    const string& file2;
    
    PatternPairMap myLCSEPM; //!PatternPairMap result of chaining algorithm
    PatternPairMap mcsPatterns;
    
    // ----------------------------------------
    // evaluate the recursions / fill matrices
    
    
    //!computes matrices A,G and B
    void compute_AGBmatrices(const ArcMatch &arc_match);

    //! computes matrix F
    void compute_F();
    
    //!helper function for compute_matrices
    infty_score_t seq_str_matching(ScoreMatrix &mat,const ArcMatch &arc_match, size_type i, size_type j,bool matrixB);
    
    //! computes score for arcMatch: score under the arcMatch plus the probability of the two basepairs of the arcMatch
    infty_score_t score_for_arc_match(const ArcMatch &am);
    
    //! computes the stacking score: if stacking occurs with respect to a structure, the stacking probability is taken as a score
    infty_score_t score_for_stacking(const ArcMatch &am, const ArcMatch &inner_am);

    //!compute the backward score and the forward pointer
    void trace_F();
    
    //!traverses matrix F
    void trace_in_F();
    
    //!adds the structure and the position corresponding to the position pos_ in the matrix if
    //!the position isn't contained in another epm (score==pos_infty)
    bool add(const ArcMatch &am,pair<int,int> pos_, char c, EPM &epm_to_store);
    
    //!adds the arcMatch to the epm if the positions aren't contained in another epm (score==pos_infty)
    bool add_arcmatch(const ArcMatch &am, EPM &epm_to_store);
    
    //!outputs the exact matching starting at position $(i,j)$ while setting the processed elements Trace(i,j) to -inf
    void get_matching(size_type i, size_type j);
    
    //!recomputes matrices A,G and B for arcMatch recursively and stores the traceback
    bool trace_AGB(const ArcMatch &am, EPM &epm_to_store);
    
    //!checks the structural case for the traceback in the matrices A,G and B
    //!if an inner_am is encountered, it is stored for later processing and the left and right endpoint is added to the epm structure 
    bool str_traceAGB(const ScoreMatrix &mat, const ArcMatch &am, size_type posA, size_type posB,pair<int,int> &curPos, EPM &epm_to_store);
    
    void set_el_to_inf();
    void set_el_to_neg_inf();
    
    void compute_F_with_prob_external();
    void find_start_pos_for_traceback(vector<pair<int,int> > &EPM_start_pos);
    void trace_in_F_suboptimal(int i, int j);
    void print_epms_to_proc(list<pair<EPM,int> > &epms_to_proc);
    bool valid_external_arcmatch(const ArcMatch &am);
    bool valid_external_pos(size_type i,size_type j);

    //!converts string to uppercase
    string upperCase(string seq){
	  string s= "";
	  for(unsigned int i= 0; i<seq.length(); i++)
	    s+= toupper(seq[i]);
	  return s;
	}
    
    //Debugging
    void output_trace_matrix();
    void output_arc_match_score();
    void print_EPM_start_pos(list<pair<pair<int,int>,infty_score_t> > &EPM_start_pos);
    void print_matrices(const ArcMatch &am, size_type offset_A, size_type offset_B);
    
public:

    //! construct with sequences and possible arc matches
    ExactMatcher(const Sequence &seqA_,const Sequence &seqB_,const ArcMatches &arc_matches_,const Mapping &mappingA_, const Mapping &mappingB_,
		 const int &threshold_,const int &min_size_,const int &alpha_1,const int &alpha_2, const int &alpha_3, const int &subopt_score,
		 const int &easier_scoring_par, const string& sequenceA_,
		 const string& sequenceB_, const string& file1_, const string& file2_);
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

    //! outputs anchor constraints to be used as input for locarna
    void
    output_locarna();
};



} //end namespace

#endif //  EXACT_MATCHER_HH
