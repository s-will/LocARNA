#include "sequence.hh"
#include "arc_matches.hh"

//#include "basepairs_looptraversal.hh"

#include "exact_matcher.hh"

#include "scoring.hh"

#include <iostream>

namespace LocARNA {

ExactMatcher::ExactMatcher(const Sequence &seqA_,
			   const Sequence &seqB_,
			   const ArcMatches &arc_matches_,
			   const Mapping &mappingA_,
			   const Mapping &mappingB_,
			   const int &EPM_threshold_,
			   const int &EPM_min_size_,
			   const int &alpha_1_,
			   const int &alpha_2_,
			   const int &alpha_3_,
			   const int &subopt_score_,
			   const int &easier_scoring_par_,
			   PatternPairMap &foundEPMs_
			   //const string& sequenceA_,
			   //const string& sequenceB_,
			   //const string& file1_,
			   //const string& file2_
			   )
    : seqA(seqA_),
      seqB(seqB_),
      arc_matches(arc_matches_),
      bpsA(arc_matches_.get_base_pairsA()),
      bpsB(arc_matches_.get_base_pairsB()),

      mappingA(mappingA_),
      mappingB(mappingB_),
    
      //epm(seqA_,seqB_),
      epm(),
      EPM_threshold(EPM_threshold_*100),
      EPM_min_size(EPM_min_size_),
      alpha_1(alpha_1_),
      alpha_2(alpha_2_),
      alpha_3(alpha_3_),
      subopt_score(subopt_score_*100),
      easier_scoring_par(easier_scoring_par_),
      foundEPMs(foundEPMs_)
     // sequenceA(sequenceA_),
     // sequenceB(sequenceB_),
      //file1(file1_),
      //file2(file2_)

{
    // set size of matrices

	A.resize(seqA.length()+1,seqB.length()+1);
	A.fill(infty_score_t::neg_infty);
	A.set(0,0,infty_score_t(0));
	G.resize(seqA.length()+1,seqB.length()+1);
	G.fill(InftyArithInt(0));
	B.resize(seqA.length()+1,seqB.length()+1);
	B.fill(InftyArithInt(0));
	F.resize(seqA.length()+1,seqB.length()+1);
	F.fill(InftyArithInt(0));
	Trace.resize(seqA.length()+2,seqB.length()+2);
	struct Trace_entry el = {InftyArithInt(0),0,0};
	Trace.fill(el);
	//myLCSEPM= PatternPairMap();
}

ExactMatcher::~ExactMatcher(){
	for(size_type i=0;i<seqA.length()+2;i++){
		for(size_type j=0;j<seqB.length()+2;j++){
			if(Trace(i,j).next_pos!=0){
				delete Trace(i,j).next_pos;
				Trace(i,j).next_pos=0;
			}
			if(Trace(i,j).arc_match_idx!=0){
				delete Trace(i,j).arc_match_idx;
				Trace(i,j).arc_match_idx=0;
			}
		}
	}
}


void ExactMatcher::compute_AGBmatrices(const ArcMatch &arc_match){
	const Arc &arcA=arc_match.arcA();
	const Arc &arcB=arc_match.arcB();
	size_type idxA = arcA.idx();
	size_type idxB = arcB.idx();
	size_type number_of_posA = mappingA.number_of_valid_pos(idxA);
	size_type number_of_posB = mappingB.number_of_valid_pos(idxB);
	InftyArithInt score_A,score_B;
	for(size_type i=1;i<number_of_posA;i++){
		for(size_type j=1;j<number_of_posB;j++){
			//don't fill in A last column and last row except the entry (number_of_posA-1,number_of_posB-1)
			if((i==number_of_posA-1 && j!=number_of_posB-1) || (i!=number_of_posA-1 && j==number_of_posB-1)){
				 A(i,j)=infty_score_t::neg_infty;
			}
			else{
				A(i,j)=max(seq_str_matching(A,arc_match,i,j,false),infty_score_t::neg_infty);
			}
			G(i,j)=max(max(A(i,j),G(i-1,j-1)),max(G(i,j-1),G(i-1,j)));
			B(i,j)=max(seq_str_matching(B,arc_match,i,j,true),G(i,j));

		}
	}
}

infty_score_t ExactMatcher::seq_str_matching(ScoreMatrix &mat, const ArcMatch &arc_match, size_type i, size_type j, bool matrixB){
	infty_score_t score_seq = infty_score_t::neg_infty; infty_score_t score_str =infty_score_t::neg_infty;
	const Arc &arcA=arc_match.arcA();
	const Arc &arcB=arc_match.arcB();
	size_type idxA = arcA.idx();
	size_type idxB = arcB.idx();
	size_type PosSeqA = mappingA.get_pos_in_seq(arcA,i);
	size_type PosSeqB = mappingB.get_pos_in_seq(arcB,j);
	if(seqA[PosSeqA]==seqB[PosSeqB]){
		if(mappingA.seq_matching(idxA,i) && mappingB.seq_matching(idxB,j)){
			//sequential matching
			score_seq = mat(i-1,j-1)+alpha_1*100;
		}
		//special case if the previous entry in matrix B was taken from matrix G -> no check if sequential matching
		//is permitted
		if(matrixB){
			score_seq = max(G(i-1,j-1)+alpha_1*100,score_seq);
		}
		//structural matching
		for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(PosSeqA-1,PosSeqB-1).begin();
						    		    			arc_matches.common_right_end_list(PosSeqA-1,PosSeqB-1).end() != it; ++it ) {
			const ArcMatch &inner_am = arc_matches.arcmatch(*it);
			if(inner_am.arcA().left()>arcA.left() && inner_am.arcB().left()>arcB.left()){
				int pos_before_arcA = mappingA.get_pos_in_new_seq(arcA,inner_am.arcA().left()-1);
				int pos_before_arcB= mappingB.get_pos_in_new_seq(arcB,inner_am.arcB().left()-1);
				if(pos_before_arcA == -1 || pos_before_arcB == -1) continue; //no valid position before the arc
				infty_score_t tmp = mat(pos_before_arcA,pos_before_arcB)+alpha_1*100+score_for_arc_match(inner_am,true)+score_for_stacking(arc_match,inner_am);
				score_str=max(score_str,tmp);
			}
		}
	}
	return max(score_seq,score_str).normalized_neg();
}

void ExactMatcher::compute_F(){
	InftyArithInt score_seq,score_str;
	for(size_type i=1;i<=seqA.length();i++){
		for(size_type j=1;j<=seqB.length();j++){
			score_seq=InftyArithInt(0);score_str=InftyArithInt(0);
			//sequential matching
			if(seqA[i]==seqB[j]){
				score_seq = F(i-1,j-1)+alpha_1*100;
			}
			//structural matching
			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
			    			arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
			  
			    const ArcMatch &am = arc_matches.arcmatch(*it);
			    score_str = max(score_str,F(am.arcA().left()-1,am.arcB().left()-1)+score_for_arc_match(am,true));
			}
			F(i,j)=max(score_seq,score_str);
		}
	}
}

void ExactMatcher::compute_matrices() {
        // for all arc matches from inside to outside    
	for(ArcMatchVec::const_iterator it=arc_matches.begin();it!=arc_matches.end();it++){
		//check if valid arc_match
		if(seqA[it->arcA().left()]!=seqB[it->arcB().left()] || seqA[it->arcA().right()]!=seqB[it->arcB().right()]){
			arc_match_score.push_back(infty_score_t::neg_infty);
		}
		else{
			//compute matrices A,G and B for an arcMatch and store it in the vector arc_match_score
			compute_AGBmatrices(*it);
			size_type number_of_posA = mappingA.number_of_valid_pos(it->arcA().idx());
			size_type number_of_posB = mappingB.number_of_valid_pos(it->arcB().idx());
			infty_score_t temp = B(number_of_posA-1,number_of_posB-1)-alpha_1*100;
			arc_match_score.push_back(temp);
		}
	}
	compute_F_with_prob_external();
}

void ExactMatcher::compute_F_with_prob_external(){
	int count_seq = 0;
	int count_am=0;
	InftyArithInt score_seq,score_str;
	for(size_type i=1;i<=seqA.length();i++){
		for(size_type j=1;j<=seqB.length();j++){
			score_seq=InftyArithInt(0);score_str=InftyArithInt(0);
			//sequential matching
			if(seqA[i]==seqB[j] && !valid_external_pos(i,j)){
				count_seq++;
			}
			if(seqA[i]==seqB[j] && valid_external_pos(i,j)){
				score_seq = F(i-1,j-1)+alpha_1*100;
			}
			//structural matching
			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
			    			arc_matches.common_right_end_list(i,j).end() != it; ++it ) {

			    const ArcMatch &am = arc_matches.arcmatch(*it);
			    if(valid_external_arcmatch(am)){
			    score_str = max(score_str,F(am.arcA().left()-1,am.arcB().left()-1)+score_for_arc_match(am,true));
			    }
			    else{
			    	count_am++;
			    }
			}
			F(i,j)=max(score_seq,score_str);
		}
	}
	cout << "positions that weren't considered because of filtering " << count_seq << endl;
	cout << "arcmatches that weren't considered because of filtering " << count_am << endl;
}

void ExactMatcher::compute_EPMs_suboptimal(){
	//output maximum in F matrix
	int max = 0;
	for(int i=0;i<seqA.length()+1;i++){
		for(int j=0;j<seqB.length()+1;j++){
			if(F(i,j).finite_value()>max) max=F(i,j).finite_value();
		}
	}
	cout << "maximal score " << max << endl;
	//this->output_arc_match_score();
	vector<pair<int,int> > EPM_start_pos;
	find_start_pos_for_traceback(EPM_start_pos);
	/*for(int i=0;i<EPM_start_pos.size();i++){
		cout << "pos " << EPM_start_pos.at(i).first << "," << EPM_start_pos.at(i).second << endl;
	}*/
	//cout << "number of starting pos " << EPM_start_pos.size() << endl;
	//cout << "suboptimal score " << subopt_score << endl;
	for(int k=0;k<EPM_start_pos.size();k++){
		int i=EPM_start_pos.at(k).first;
		int j=EPM_start_pos.at(k).second;
		trace_in_F_suboptimal(i,j);
	}

	//trace_in_F_suboptimal(14,18);

	// chaining for suboptimal traceback case
//	const int& size1= (int)seqA.length();
//		const int& size2= (int)seqB.length();
//		int size= 0;
//		int score= 0;
//		cout << "#EPM: " << mcsPatterns.size() << endl;
//		time_t start_chaining = time (NULL);
//		//create LCSEPM object
//		//LCSEPM patterns(size1, size2 ,epm.get_patternPairMap(), myLCSEPM,size,score, EPM_min_size);
//		LCSEPM patterns(size1, size2 ,mcsPatterns, myLCSEPM,size,score, EPM_min_size);
//		//begin chaining algorithm
//		patterns.calculateLCSEPM();
//		time_t stop_chaining = time (NULL);
//	    cout << "time for chaining : " << stop_chaining - start_chaining << "sec " << endl;
//		//output patterns to PS files
//	    time_t start_ps = time (NULL);
//		patterns.MapToPS(sequenceA, sequenceB, size, myLCSEPM, file1, file2);
//		time_t stop_ps = time (NULL);
//		cout << "time for map to ps : " << stop_ps - start_ps << "sec " << endl;
}

void ExactMatcher::find_start_pos_for_traceback(vector<pair<int,int> > &EPM_start_pos){
	bool valid=true;
	cout << "compute EPMs suboptimal " << endl;
	cout << "suboptimal score " << subopt_score << endl;
	//store start points for maximally extended exact pattern matchings -> all patterns that cannot be extended
	for(int i=0;i<seqA.length()+1;i++){
		for(int j=0;j<seqB.length()+1;j++){
			if(F(i,j)>(infty_score_t)subopt_score){
				valid=true;
				if(i==seqA.length() || j==seqB.length()){
				}
				else if(seqA[i+1]==seqB[j+1] && valid_external_pos(i+1,j+1)){
					//don't consider (i,j)
					valid =false;
				}
				else{ //TODO not needed?
					for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
							arc_matches.common_left_end_list(i,j).end() != it; ++it){
						const ArcMatch &am = arc_matches.arcmatch(*it);
						if(!(arc_match_score.at(am.idx())==infty_score_t::neg_infty)
								&& valid_external_arcmatch(am)
						){
							//don't consider (i,j)
							valid=false;
							break;
						}
					}
				}
				if(valid){
					EPM_start_pos.push_back(pair<int,int>(i,j));
				}
			}
		}
	}
	cout << "Number of starting positions: " << EPM_start_pos.size() << endl;
}

bool ExactMatcher::valid_external_arcmatch(const ArcMatch &am){
	 return (mappingA.basepair_external(am.arcA().left(),am.arcA().right())
			&& mappingB.basepair_external(am.arcB().left(),am.arcB().right()));
}

bool ExactMatcher::valid_external_pos(size_type i, size_type j){
	return (mappingA.unpaired_external(i) && mappingB.unpaired_external(j));
}

void ExactMatcher::trace_in_F_suboptimal(int i,int j){
	//cout << "trace F " << i << "," << j << endl;
	static int count=0;

	static string seq1_id = seqA.names()[0];
	static string seq2_id = seqB.names()[0];

	list<info_for_trace_F > epms_to_proc;
	struct info_for_trace_F tmp = {epm,0,pair<int,int>(i,j)};
	epms_to_proc.push_back(tmp);
	int cur_i,cur_j,cur_score;

	while(epms_to_proc.begin()!=epms_to_proc.end()){
		epm = epms_to_proc.front().cur_epm;
		cur_i = epms_to_proc.front().curPos.first;
		cur_j = epms_to_proc.front().curPos.second;
		cur_score = epms_to_proc.front().score;
		if(!(F(cur_i,cur_j).finite_value()==0)){
			if(cur_score+alpha_1*100+F(cur_i-1,cur_j-1).finite_value()>subopt_score
			   && valid_external_pos(cur_i,cur_j)){
				epm.add(cur_i,cur_j,'.');
				struct info_for_trace_F tmp = {epm,cur_score+alpha_1*100,pair<int,int>(cur_i-1,cur_j-1)};
				epms_to_proc.push_back(tmp);
			}
			for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(cur_i,cur_j).begin();
							arc_matches.common_right_end_list(cur_i,cur_j).end() != it; ++it){
				const ArcMatch &am = arc_matches.arcmatch(*it);
				if(score_for_arc_match(am,true)==infty_score_t::neg_infty){
					//don't consider arcmatch
					break;
				}
				else if(valid_external_arcmatch(am) && cur_score + score_for_arc_match(am,true).finite_value() +
						F(am.arcA().left()-1,am.arcB().left()-1).finite_value() > subopt_score){
					// if score better than subopt score
					epm = epms_to_proc.front().cur_epm;
					epm.add_arcmatch(am);

					//suboptimal traceback
					list<info_for_trace_AGB> epms_to_proc_AGB;
					int posA = mappingA.get_pos_in_new_seq(am.arcA(),cur_i);
					int posB = mappingB.get_pos_in_new_seq(am.arcB(),cur_j);
					struct info_for_trace_AGB initial_info={epm,epms_to_proc.front().score+score_for_arc_match(am,false).finite_value(),in_B,pair<int,int>(posA,posB),epm.get_arcmatches_to_do()};
					epms_to_proc_AGB.push_back(initial_info);
					epms_to_proc_AGB.push_back(initial_info);
					//if(am.idx()==607){
					//trace_AGB_suboptimal_main(am,epm,epms_to_proc_AGB,epms_to_proc);
					//}

					//optimal traceback
					trace_AGB(am,epm);
					struct info_for_trace_F tmp = {epm,epms_to_proc.front().score+score_for_arc_match(am,true).finite_value(),pair<int,int>(am.arcA().left()-1,am.arcB().left()-1)};
					epms_to_proc.push_back(tmp);
				}
			}
		}
		else{
			//epm is computed
			count++;
			epm.sort_patVec();
			//epm.validate_epm();
			//epm.print_epm(cout,cur_score);
			
			stringstream ss;
			ss << "pat_" << count;
			string patId= ss.str();

			SinglePattern pattern1 = SinglePattern(patId,seq1_id,epm.getPat1Vec());
			SinglePattern pattern2 = SinglePattern(patId,seq2_id,epm.getPat2Vec());

			foundEPMs.add(patId, pattern1, pattern2, epm.getStructure(), cur_score);
		}
		epms_to_proc.pop_front();
	}
	epm.reset();
	//cout << "number of EPM " << count << endl;
}

void
ExactMatcher::trace_AGB_suboptimal_main(const ArcMatch &am, EPM &epm_to_store,list<info_for_trace_AGB> &epms_to_proc_AGB,list<info_for_trace_F > &epms_to_proc){
	//while(epms_to_proc_AGB.size()!=0){

		trace_AGB_suboptimal(am,epm_to_store,epms_to_proc_AGB, epms_to_proc);
		epms_to_proc_AGB.pop_front();
	//}
}


//TODO stacking probabilities
bool
ExactMatcher::trace_AGB_suboptimal(const ArcMatch &am, EPM &epm_to_store,list<info_for_trace_AGB> &epms_to_proc_AGB,list<info_for_trace_F > &epms_to_proc){
	compute_AGBmatrices(am);
	//this->print_matrices(am,20,20);
	cout << "trace AGB suboptimal " << endl;
	size_type number_of_posA = mappingA.number_of_valid_pos(am.arcA().idx());
	size_type number_of_posB = mappingB.number_of_valid_pos(am.arcB().idx());
	pair<int,int> str_pos;
	pair<int,int> curPos;
	int i,j;
	size_type posA,posB;
	struct info_for_trace_AGB cur_epm_AGB;
	int state=in_B;
	int cur_score;
	//Initialization
	epms_to_proc_AGB.pop_front();
	cur_epm_AGB = epms_to_proc_AGB.front();
	state = cur_epm_AGB.state;
	cur_score = cur_epm_AGB.score;
	pair<int,int> cur_pos_seq = pair<int,int>(cur_epm_AGB.cur_epm.getPat1Vec().back(),cur_epm_AGB.cur_epm.getPat2Vec().back());
	curPos=cur_epm_AGB.curPos;
	epm = cur_epm_AGB.cur_epm;

	while((curPos!=pair<int,int>(0,0) || state!=in_A)){
		i=curPos.first;
		j=curPos.second;
		posA = mappingA.get_pos_in_seq(am.arcA(),i);
		posB = mappingB.get_pos_in_seq(am.arcB(),j);
		switch(state){
			case in_B:
			{
				if(seqA[posA]==seqB[posB]){
					//structural matching
					for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(posA-1,posB-1).begin();
								arc_matches.common_right_end_list(posA-1,posB-1).end() != it; ++it ){
							const ArcMatch &inner_am = arc_matches.arcmatch(*it);
							if(inner_am.arcA().left()>am.arcA().left() && inner_am.arcB().left()>am.arcB().left()){
								int pos_before_arcA = mappingA.get_pos_in_new_seq(am.arcA(),inner_am.arcA().left()-1);
								int pos_before_arcB= mappingB.get_pos_in_new_seq(am.arcB(),inner_am.arcB().left()-1);
								if(pos_before_arcA == -1 || pos_before_arcB == -1) continue; //no valid position before the arc
								bool flag = (subopt_score < (score_for_arc_match(inner_am,true)+score_for_stacking(am,inner_am)+alpha_1*100+B(pos_before_arcA,pos_before_arcB)).finite_value());
								if(flag){
									cout << "structural matching " << endl;
									if(posA!=am.arcA().right() && posB!=am.arcB().right()){
										//add paired bases to structure and corresponding position
										curPos=pair<int,int>(pos_before_arcA,pos_before_arcB);
										epm.add(posA,posB,'.');
										epm.add_arcmatch(inner_am);
										struct info_for_trace_AGB tmp = {epm,cur_score+score_for_arc_match(inner_am,false).finite_value(),state,curPos,epm.get_arcmatches_to_do()};//;.push_back(inner_am.idx())}; //TODO cur score
										epms_to_proc_AGB.push_back(tmp);
										epm = cur_epm_AGB.cur_epm;
										epm.get_arcmatches_to_do().pop_back();
									}
								}
							}
						}
					//sequential matching
					bool flag = (subopt_score < cur_score+alpha_1*100+B(i-1,j-1).finite_value());
					if(flag && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						int add_score=0;
						if(posA!=am.arcA().right() && posB!=am.arcB().right()){
							epm.add(posA,posB,'.');
							add_score = alpha_1*100;
						}
						curPos=pair<int,int>(i-1,j-1);
						struct info_for_trace_AGB tmp = {epm,cur_score+add_score,state,curPos,epm.get_arcmatches_to_do()};
						epms_to_proc_AGB.push_back(tmp);
						epm = cur_epm_AGB.cur_epm;

					}
				}
				if(B(i,j)==G(i,j)){
					state=in_G; break;
				}
				bool flag2 = B(i,j)==G(i-1,j-1)+alpha_1*100 && !(mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j));
				//special case
				if(flag2){
					curPos=pair<int,int>(i-1,j-1);
					if(posA!=am.arcA().right() && posB!=am.arcB().right()){
						epm.add(posA,posB,'.');
					}
					state=in_G;
					break;
				}
			}
			case in_G:
				if(G(i,j)==A(i,j)){
						struct info_for_trace_AGB tmp = {epm,cur_score,in_A,curPos,epm.get_arcmatches_to_do()};
						epms_to_proc_AGB.push_back(tmp);
						epm = cur_epm_AGB.cur_epm;
					state=in_A; break;
				}
				if(i>=1 && j>=1 && G(i,j)==G(i-1,j-1)){
					curPos=pair<int,int>(i-1,j-1);
					break;
				}
				if(i>=1 && G(i,j)==G(i-1,j)){
					curPos=pair<int,int>(i-1,j);
					break;
				}
				if(j>=1 && G(i,j)==G(i,j-1)){
					curPos=pair<int,int>(i,j-1);
					break;
				}
			case in_A:
			{
				if(seqA[mappingA.get_pos_in_seq(am.arcA(),i)]==seqB[mappingB.get_pos_in_seq(am.arcB(),j)]){
					//structural matching
					for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(posA-1,posB-1).begin();
									arc_matches.common_right_end_list(posA-1,posB-1).end() != it; ++it ){
						const ArcMatch &inner_am = arc_matches.arcmatch(*it);
						if(inner_am.arcA().left()>am.arcA().left() && inner_am.arcB().left()>am.arcB().left()){
							int pos_before_arcA = mappingA.get_pos_in_new_seq(am.arcA(),inner_am.arcA().left()-1);
							int pos_before_arcB= mappingB.get_pos_in_new_seq(am.arcB(),inner_am.arcB().left()-1);
							if(pos_before_arcA == -1 || pos_before_arcB == -1) continue; //no valid position before the arc
							bool flag = (subopt_score < (score_for_arc_match(inner_am,true)+score_for_stacking(am,inner_am)+alpha_1*100+A(pos_before_arcA,pos_before_arcB)).finite_value());
							if(flag){
								cout << "structural matching " << endl;
								curPos=pair<int,int>(pos_before_arcA,pos_before_arcB);
								if(posA!=am.arcA().right() && posB!=am.arcB().right()){
									epm.add(posA,posB,'.');
									epm.add_arcmatch(inner_am);
									struct info_for_trace_AGB tmp = {epm,cur_score+score_for_arc_match(inner_am,false).finite_value(),state,curPos,epm.get_arcmatches_to_do()};//TODO
									epms_to_proc_AGB.push_back(tmp);
									epm = cur_epm_AGB.cur_epm;
								}
							}
						}
					}
					bool flag3 = (subopt_score < cur_score+alpha_1*100+A(i-1,j-1).finite_value());
					if(flag3 && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						curPos=pair<int,int>(i-1,j-1);
						if(posA!=am.arcA().right() && posB!=am.arcB().right()){
							epm.add(posA,posB,'.');
							struct info_for_trace_AGB tmp = {epm,cur_score+alpha_1*100,state,curPos,epm.get_arcmatches_to_do()};
							epms_to_proc_AGB.push_back(tmp);
						}
						break;
					}
				}
			}
		}
		if(state!=in_G){
			if(epms_to_proc_AGB.size()>1){
				epms_to_proc_AGB.pop_front();
				cur_epm_AGB = epms_to_proc_AGB.front();
			}
			else{
				break;
			}
			state = cur_epm_AGB.state;
			cur_score = cur_epm_AGB.score;
			pair<int,int> cur_pos_seq = pair<int,int>(cur_epm_AGB.cur_epm.getPat1Vec().back(),cur_epm_AGB.cur_epm.getPat2Vec().back());
			curPos=cur_epm_AGB.curPos;
			epm = cur_epm_AGB.cur_epm;
		}
	}
	struct info_for_trace_F tmp = {epm,cur_score,pair<int,int>(am.arcA().left()-1,am.arcB().left()-1)};
	epms_to_proc.push_back(tmp);
	return true;
}

void ExactMatcher::print_epms_to_proc_AGB(list<info_for_trace_AGB> &epms_to_proc_AGB){
	cout << "epms to proc " << endl;
	for(list<info_for_trace_AGB>::iterator it = epms_to_proc_AGB.begin();it!=epms_to_proc_AGB.end();it++){
		it->cur_epm.print_epm(cout,it->score);
		cout << "curPos " << it->curPos.first << "," << it->curPos.second << endl;
		cout << "curState " << it->state << endl;
		print_arcmatches_to_do(it->arcmatches_to_do_for_cur_epm);
		cout << endl;
	}
}

void ExactMatcher::print_arcmatches_to_do(std::vector<ArcMatch::idx_type> arcmatches_to_do){
	if(arcmatches_to_do.size()!=0){
		cout << "arcmatches to do " << endl;
		for(std::vector<ArcMatch::idx_type>::iterator it2 = arcmatches_to_do.begin();it2!=arcmatches_to_do.end();it2++){
			cout << *it2 << endl;
		}
	}
}

void ExactMatcher::print_epms_to_proc(list<info_for_trace_F > &epms_to_proc){
	cout << "printing epms to proc " << endl;
	int count =0;
	for(list<info_for_trace_F >::iterator it=epms_to_proc.begin();it!=epms_to_proc.end();it++){
		cout << count << ".EPM: ";
		EPM cur_epm = it->cur_epm;
		cur_epm.print_epm(cout,it->score);
		cout << "curPos " << it->curPos.first << "," << it->curPos.second << endl;
		count++;
	}
}

void ExactMatcher::compute_EPMs_heuristic() {
    trace_F();
    cout << "trace_F " << endl;
    trace_in_F();
}

//TODO: modifiy trace_F with prob_external!
void
ExactMatcher::trace_F(){
	for(int i=seqA.length();i>=1;i--){
		for(int j=seqB.length();j>=1;j--){
			//sequential matching
			bool flag = F(i,j)==F(i-1,j-1)+alpha_1*100;
			if(flag){
				bool flag2 = Trace(i,j).score>Trace(i+1,j+1).score+alpha_1*100;
				if(flag2){
					//don't consider traceback from position (i+1,j+1)
					Trace(i+1,j+1).score=infty_score_t::neg_infty;
				}
				else{
					pair<int,int> *next = Trace(i,j).next_pos;
					//set score of the element to which the old pointer (if it exists) points to to -inf
					//reset pointer
					if(next!=0){
						Trace(next->first,next->second).score=infty_score_t::neg_infty;
						*(Trace(i,j).next_pos)=pair<int,int>(i+1,j+1);
						
					}
					else{
						Trace(i,j).next_pos=new pair<int,int>(i+1,j+1);
					}
					//delete pointer to arc_match_idx if it exists
					if(Trace(i,j).arc_match_idx!=0){
						delete Trace(i,j).arc_match_idx;
						Trace(i,j).arc_match_idx=0;
					}
					
					//set score
					Trace(i,j).score= Trace(i+1,j+1).score+alpha_1*100;

				}
				
			}
			else{
				for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
										    			arc_matches.common_right_end_list(i,j).end() != it; ++it){
					const ArcMatch &am = arc_matches.arcmatch(*it);
					if(F(i,j)==F(am.arcA().left()-1,am.arcB().left()-1)+score_for_arc_match(am,true)){
						if(Trace(am.arcA().left(),am.arcB().left()).score>Trace(i+1,j+1).score+score_for_arc_match(am,true)){
							//don't consider traceback from position (i+1,j+1)
							Trace(i+1,j+1).score=infty_score_t::neg_infty;
						}
						else{
							pair<int,int> *next=Trace(am.arcA().left(),am.arcB().left()).next_pos;
							//set score of the element to which the old pointer (if it exists) points to to -inf
							//reset pointer
							if(next!=0){
								Trace(next->first,next->second).score=infty_score_t::neg_infty;
								*(Trace(am.arcA().left(),am.arcB().left()).next_pos)=pair<int,int>(i,j);
							}
							else{
								Trace(am.arcA().left(),am.arcB().left()).next_pos= new pair<int,int>(i,j);
							}
							//reset pointer to arc_match_idx
							if(Trace(am.arcA().left(),am.arcB().left()).arc_match_idx!=0){
								*(Trace(am.arcA().left(),am.arcB().left())).arc_match_idx=ArcMatch::idx_type(*it);
							}
							else{
								Trace(am.arcA().left(),am.arcB().left()).arc_match_idx = new ArcMatch::idx_type(*it);
							}
							//set score
							Trace(am.arcA().left(),am.arcB().left()).score = Trace(i+1,j+1).score+score_for_arc_match(am,true);
						}
					}
				}
			}
		}
	}
}

bool compare(pair<pair<int,int>,infty_score_t > entry1,pair<pair<int,int>,infty_score_t > entry2) {
	return entry1.second > entry2.second;
}

void
ExactMatcher::trace_in_F(){
	list<pair<pair<int,int>,infty_score_t> > EPM_start_pos;
	for(size_type i=0;i<seqA.length()+1;i++){
		for(size_type j=0;j<seqB.length()+1;j++){
			if(Trace(i,j).score>(infty_score_t)EPM_threshold){
				//consider element for EPM
				pair<int,int> pos = pair<int,int>(i,j);
				EPM_start_pos.push_back(pair<pair<int,int>,infty_score_t >(pos,Trace(i,j).score));
			}
		}
	}
	//sort EPMs according to their score
	EPM_start_pos.sort(compare);
	while(EPM_start_pos.size()!=0){
		get_matching(EPM_start_pos.front().first.first,EPM_start_pos.front().first.second);
		EPM_start_pos.pop_front();
	}


	//this->output_trace_matrix();
//	const int& size1= (int)seqA.length();
//	const int& size2= (int)seqB.length();
//	int size= 0;
//	int score= 0;
//	cout << "#EPM: " << mcsPatterns.size() << endl;
//	time_t start_chaining = time (NULL);
//	//create LCSEPM object
//	//LCSEPM patterns(size1, size2 ,epm.get_patternPairMap(), myLCSEPM,size,score, EPM_min_size);
//	LCSEPM patterns(size1, size2 ,mcsPatterns, myLCSEPM,size,score, EPM_min_size);
//	//begin chaining algorithm
//	patterns.calculateLCSEPM();
//	time_t stop_chaining = time (NULL);
//    cout << "time for chaining : " << stop_chaining - start_chaining << "sec " << endl;
//	//output patterns to PS files
//    time_t start_ps = time (NULL);
//	patterns.MapToPS(sequenceA, sequenceB, size, myLCSEPM, file1, file2);
//	time_t stop_ps = time (NULL);
//	cout << "time for map to ps : " << stop_ps - start_ps << "sec " << endl;
}

void
ExactMatcher::get_matching(size_type i, size_type j){
    bool valid=true;
	static int count= 0;
	static string seq1_id = seqA.names()[0];
	static string seq2_id = seqB.names()[0];
	pair<int,int> prevEl=pair<int,int>(i,j);
	pair<int,int> *curEl=Trace(i,j).next_pos;
	infty_score_t score1= Trace(i,j).score;
	if(!(Trace(i,j).score==infty_score_t::neg_infty)){
		while(curEl!=0){
			//structural pointer
			if(Trace(prevEl.first,prevEl.second).arc_match_idx!=0){
				ArcMatch::idx_type idx = *(Trace(prevEl.first,prevEl.second).arc_match_idx);
				const ArcMatch &am = arc_matches.arcmatch(idx);
				epm.add_arcmatch(am);
				if(!trace_AGB(arc_matches.arcmatch(idx),epm)) {
					valid=false;
					break;
				}
				//reset cur_it to the last element
				prevEl=pair<int,int>(am.arcA().right()+1,am.arcB().right()+1);
				curEl = Trace(prevEl.first,prevEl.second).next_pos;
			}
			else{
				epm.add(prevEl.first, prevEl.second, '.');
				prevEl=*curEl;
				curEl = Trace(curEl->first,curEl->second).next_pos;
			}
		}
		epm.sort_patVec();
		if(valid){
		  //epm.validate_epm();
		  set_el_to_neg_inf();
		  count++;
		  stringstream ss;
		  ss << "pat_" << count;
		  string patId= ss.str();
		  SinglePattern pattern1 = SinglePattern(patId,seq1_id,epm.getPat1Vec());
		  SinglePattern pattern2 = SinglePattern(patId,seq2_id,epm.getPat2Vec());
		  int score= score1.finite_value();
		  //epm.add_pattern(patId,pattern1,pattern2, score);
		  foundEPMs.add(patId, pattern1, pattern2, epm.getStructure(), score);
		}
		else{set_el_to_neg_inf();}
		epm.reset();
	}
	  
}

/*void ExactMatcher::set_el_to_inf(){
	for(int i=0;i<epm.getPat1Vec().size();i++){
		int pos1=epm.getPat1Vec().at(i);
		int pos2=epm.getPat2Vec().at(i);
		//Trace(pos1,pos2).score = infty_score_t::pos_infty;
		Trace(pos1,pos2).score= infty_score_t::neg_infty;
	}
}*/

void ExactMatcher::set_el_to_neg_inf(){
	if(!epm.isEmpty()){
		for(int i=0;i<epm.getPat1Vec().size();i++){
			int pos1=epm.getPat1Vec().at(i);
			int pos2=epm.getPat2Vec().at(i);
			if(!(Trace(pos1,pos2).score==infty_score_t::pos_infty)){
			Trace(pos1,pos2).score = infty_score_t::neg_infty;
			}
		}
		//reset also parts after the epm to -inf in the F matrix
		infty_score_t score_after = Trace(epm.getPat1Vec().at(epm.getPat1Vec().size()-1)+1,epm.getPat2Vec().at(epm.getPat2Vec().size()-1)+1).score;
		pair<int,int> curEl = pair<int,int>(epm.getPat1Vec().at(epm.getPat1Vec().size()-1)+1,epm.getPat2Vec().at(epm.getPat2Vec().size()-1)+1);
		pair<int,int> *next = Trace(curEl.first,curEl.second).next_pos;
		if(!(score_after == infty_score_t::pos_infty || score_after == infty_score_t::neg_infty || score_after == (infty_score_t)0)){
			while(next!=0){
				//arc match
				if(curEl.first != next->first-1 || curEl.second != next->second-1){
					if(!(Trace(curEl.first,curEl.second).score == infty_score_t::pos_infty)){
						Trace(curEl.first,curEl.second).score=infty_score_t::neg_infty;
					}
					if(!(Trace(next->first,next->second).score == infty_score_t::pos_infty)){
						Trace(next->first,next->second).score=infty_score_t::neg_infty;
					}
					curEl=pair<int,int>(next->first+1,next->second+1);
				}
				else{
					if(!(Trace(curEl.first,curEl.second).score == infty_score_t::pos_infty)){
						Trace(curEl.first,curEl.second).score=infty_score_t::neg_infty;
					}
					curEl = *next;
				}
				next = Trace(next->first,next->second).next_pos;
			}
		}
		//epm.reset();
	}
}


bool ExactMatcher::str_traceAGB(const ScoreMatrix &mat, const ArcMatch &am, size_type i, size_type j,pair<int,int> &curPos,EPM &epm_to_store){
	size_type posA = mappingA.get_pos_in_seq(am.arcA(),i);
	size_type posB = mappingB.get_pos_in_seq(am.arcB(),j);
	curPos=pair<int,int>(-1,-1);
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(posA-1,posB-1).begin();
			arc_matches.common_right_end_list(posA-1,posB-1).end() != it; ++it ){
		const ArcMatch &inner_am = arc_matches.arcmatch(*it);
		if(inner_am.arcA().left()>am.arcA().left() && inner_am.arcB().left()>am.arcB().left()){
			int pos_before_arcA = mappingA.get_pos_in_new_seq(am.arcA(),inner_am.arcA().left()-1);
			int pos_before_arcB= mappingB.get_pos_in_new_seq(am.arcB(),inner_am.arcB().left()-1);
			if(pos_before_arcA == -1 || pos_before_arcB == -1) continue; //no valid position before the arc
			bool flag = mat(i,j)==score_for_arc_match(inner_am,true)+score_for_stacking(am,inner_am)+alpha_1*100+mat(pos_before_arcA,pos_before_arcB);
			if(flag){
				curPos=pair<int,int>(pos_before_arcA,pos_before_arcB);
				if(posA!=am.arcA().right() && posB!=am.arcB().right()){
					//add paired bases to structure and corresponding position
					epm.add(posA,posB,'.');
				}
				return add_arcmatch(inner_am,epm_to_store);
				break;
			}
		}
	}
	return true;
}

bool
ExactMatcher::trace_AGB(const ArcMatch &am, EPM &epm_to_store){
	compute_AGBmatrices(am);
	size_type number_of_posA = mappingA.number_of_valid_pos(am.arcA().idx());
	size_type number_of_posB = mappingB.number_of_valid_pos(am.arcB().idx());
	int state=in_B;
	pair<int,int> str_pos;
	pair<int,int> curPos = pair<int,int>(number_of_posA-1,number_of_posB-1);
	while((curPos!=pair<int,int>(0,0) || state!=in_A)){
		int i=curPos.first;
		int j=curPos.second;
		switch(state){
			case in_B:
			{
				if(seqA[mappingA.get_pos_in_seq(am.arcA(),i)]==seqB[mappingB.get_pos_in_seq(am.arcB(),j)]){
					//structural matching
					if(!str_traceAGB(B,am,i,j,str_pos,epm_to_store)) return false;
					if(str_pos!=pair<int,int>(-1,-1)){
						curPos=str_pos; break;
					}
					//sequential matching
					bool flag = B(i,j)==B(i-1,j-1)+alpha_1*100;
					if(flag && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						curPos=pair<int,int>(i-1,j-1);
						if(!add(am,pair<int,int>(i,j),'.',epm_to_store)) return false;
						break;
					}
				}
				if(B(i,j)==G(i,j)){
					state=in_G; break;
				}
				bool flag2 = B(i,j)==G(i-1,j-1)+alpha_1*100;
				//special case
				if(flag2){
					curPos=pair<int,int>(i-1,j-1);
					if(!add(am,pair<int,int>(i,j),'.',epm_to_store)) return false;
					state=in_G;
					break;
				}
			}
			case in_G:
				if(G(i,j)==A(i,j)){
					state=in_A; break;
				}
				if(i>=1 && j>=1 && G(i,j)==G(i-1,j-1)){
					curPos=pair<int,int>(i-1,j-1);
					break;
				}
				if(i>=1 && G(i,j)==G(i-1,j)){
					curPos=pair<int,int>(i-1,j);
					break;
				}
				if(j>=1 && G(i,j)==G(i,j-1)){
					curPos=pair<int,int>(i,j-1);
					break;
				}
			case in_A:
			{
				if(seqA[mappingA.get_pos_in_seq(am.arcA(),i)]==seqB[mappingB.get_pos_in_seq(am.arcB(),j)]){
					//structural matching
					if(!str_traceAGB(A,am,i,j,str_pos,epm_to_store)) return false;
					if(str_pos!=pair<int,int>(-1,-1)){
						curPos=str_pos; break;
					}
					//sequential matching
					bool flag3 = A(i,j)==A(i-1,j-1)+alpha_1*100;
					if(flag3 && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						curPos=pair<int,int>(i-1,j-1);
						if(!add(am,pair<int,int>(i,j),'.',epm_to_store)) return false;
						break;
					}
				}
			}
		}
	}
	//if there are arcMatches left to process, the last arcMatch is processed next
	if(epm.arcmatch_to_process()){
		return trace_AGB(arc_matches.arcmatch(epm.next_arcmatch()),epm_to_store);
	}
	return true;
}


bool ExactMatcher::add_arcmatch(const ArcMatch &am,EPM &epm_to_store){
	//positions are already contained in an EPM
	if(Trace(am.arcA().right(),am.arcB().right()).score==infty_score_t::pos_infty
			|| Trace(am.arcA().left(),am.arcB().left()).score==infty_score_t::pos_infty){
		return false;
	}
	//add arcMatch to the EPM
	epm_to_store.add_arcmatch(am);
	epm_to_store.store_arcmatch(am.idx());
	//set the score of the corresponding positions to +inf
	//Trace(am.arcA().right(),am.arcB().right()).score=infty_score_t::pos_infty;
	//Trace(am.arcA().left(),am.arcB().left()).score=infty_score_t::pos_infty;
	return true;
}


bool ExactMatcher::add(const ArcMatch &am, pair<int,int> pos_, char c,EPM &epm_to_store){
	size_type posA = mappingA.get_pos_in_seq(am.arcA(),pos_.first);
	size_type posB = mappingB.get_pos_in_seq(am.arcB(),pos_.second);
	int number_of_posA = mappingA.number_of_valid_pos(am.arcA().idx());
	int number_of_posB = mappingB.number_of_valid_pos(am.arcB().idx());
	if(pos_.first==number_of_posA-1 && pos_.second==number_of_posB-1) return true;
	//positions are already contained in an EPM
	else if(Trace(posA,posB).score==infty_score_t::pos_infty) return false;
	else{
		//add position and corresponding structure
		epm_to_store.add(posA,posB,c);
		//set the score of the corresponding position to +inf
		//Trace(posA,posB).score=infty_score_t::pos_infty;
		return true;
	}
}



infty_score_t ExactMatcher::score_for_arc_match(const ArcMatch &am, bool with_part_under_am){
	if(easier_scoring_par){
		if(with_part_under_am){
			return arc_match_score.at(am.idx())+((2*alpha_1)+easier_scoring_par*2)*100;
		}
		else{
			return infty_score_t(0)+((2*alpha_1)+easier_scoring_par*2)*100;
		}
	}
	double probArcA = bpsA.get_arc_prob(am.arcA().left(),am.arcA().right());
	double probArcB = bpsB.get_arc_prob(am.arcB().left(),am.arcB().right());
	if(with_part_under_am){
	return arc_match_score.at(am.idx())+((2*alpha_1)+(probArcA+probArcB)*alpha_2)*100;
	}
	else{
		return infty_score_t(0)+((2*alpha_1)+(probArcA+probArcB)*alpha_2)*100;
	}
}


infty_score_t ExactMatcher::score_for_stacking(const ArcMatch &am, const ArcMatch &inner_am){
    double prob_stacking_arcA = 0;
	double prob_stacking_arcB = 0;
	//stacking arcA
	if(am.arcA().left()+1==inner_am.arcA().left() &&
	   am.arcA().right()-1==inner_am.arcA().right()){
	     prob_stacking_arcA = bpsA.get_arc_2_prob(am.arcA().left(),am.arcA().right());
	}
	//stacking arcB 
	if(am.arcB().left()+1==inner_am.arcB().left() &&
	   am.arcB().right()-1==inner_am.arcB().right()){
	     prob_stacking_arcB = bpsB.get_arc_2_prob(am.arcB().left(),am.arcB().right());
	}
	return infty_score_t(0)+(prob_stacking_arcA+prob_stacking_arcB)*100*alpha_3;
}

void ExactMatcher::print_matrices(const ArcMatch &am, size_type offsetA,size_type offsetB){
	size_type number_of_posA = mappingA.number_of_valid_pos(am.arcA().idx());
	size_type number_of_posB = mappingB.number_of_valid_pos(am.arcB().idx());
	if(offsetA>number_of_posA){offsetA=number_of_posA;}
	if(offsetB>number_of_posB){offsetB=number_of_posB;}
	cout << "number of pos A " << number_of_posA << endl;
	cout << "number of pos B " << number_of_posB << endl;

	cout << "A" << endl;
	for(int i=number_of_posA-offsetA;i<number_of_posA;i++){
		for(int j=number_of_posB-offsetB;j<number_of_posB;j++){
			cout << A(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "G" << endl;
	for(int i=number_of_posA-offsetA;i<number_of_posA;i++){
		for(int j=number_of_posB-offsetB;j<number_of_posB;j++){
			cout << G(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "B" << endl;
	for(int i=number_of_posA-offsetA;i<number_of_posA;i++){
		for(int j=number_of_posB-offsetB;j<number_of_posB;j++){
			cout << B(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void ExactMatcher::validate_epm(){
	//EPMs have the same size
	if(epm.getPat1Vec().size()!=epm.getPat2Vec().size()){
		cerr << "wrong pattern " << endl;
		return;
	}

	//two matched positions in the EPMs have the same nucleotide
	for(int i=0;i<epm.getPat1Vec().size();i++){
		if(seqA[epm.getPat1Vec().at(i)]!=seqB[epm.getPat2Vec().at(i)]){
			cerr << "no EPM " << endl;
		}
	}
	//test for Pat1Vec and Pat2Vec
	intVec vec;
	for(int i=0;i<2;i++){
		if(i==0) vec=epm.getPat1Vec();
		if(i==1) vec=epm.getPat2Vec();
		vector<pair<int,int> > arcmatches_to_validate;
		int balance=0;
		bool gap = false;
		arcmatches_to_validate.push_back(pair<int,int>(0,vec.size()-1));
		while(arcmatches_to_validate.size()!=0){
			pair<int,int> cur_arcmatch = arcmatches_to_validate.back();
			arcmatches_to_validate.pop_back();
			gap=false;
			for(int i=cur_arcmatch.first+1;i<=cur_arcmatch.second;i++){
				if(epm.getStructure().at(i)=='.'){
					if(i!=0 && !(vec.at(i-1)+1==vec.at(i))){
						if(gap) cerr << "ERROR";
						else gap=true;
					}
				}
				else if(epm.getStructure().at(i)=='('){
					int start = i;
					if(!(vec.at(i-1)+1==vec.at(i))){
						if(gap) cerr << "ERROR " << endl;
						gap=true;
					}
					int balance = 1; //count first open bracket at position i
					while(balance!=0){
						i++;
						if(epm.getStructure().at(i)=='(') balance++;
						if(epm.getStructure().at(i)==')') balance--;
					}
					arcmatches_to_validate.push_back(pair<int,int>(start,i));
				}
				else if(epm.getStructure().at(i)==')'){
					if(!(vec.at(i-1)+1==vec.at(i))){
						if(gap) cerr << "ERROR " << endl;
						gap=true;
					}
				}
			}
		}
	}
}


//for debugging
void ExactMatcher::print_EPM_start_pos(list<pair<pair<int,int>,infty_score_t> > &EPM_start_pos){
	cout << "EPM start pos" << endl;
	for(list<pair<pair<int,int>,infty_score_t> >::iterator it=EPM_start_pos.begin();it!=EPM_start_pos.end();it++){
		cout << it->first.first << "," << it->first.second << ": " << it->second << endl;
	}
}

void ExactMatcher::output_trace_matrix(){
	std::cout << "Trace Score" << std::endl;
	for(size_type i=0;i<=seqA.length()+1;i++){
		for(size_type j=0;j<=seqB.length()+1;j++){
			cout << Trace(i,j).score << "\t";
		}
		cout << endl;
	}
	cout << endl;
	std::cout << "Trace Pointer" << std::endl;
		for(size_type i=0;i<=seqA.length()+1;i++){
			for(size_type j=0;j<=seqB.length()+1;j++){
				if(Trace(i,j).next_pos!=0) cout << Trace(i,j).next_pos->first << "," << Trace(i,j).next_pos->second << "\t";
				else cout << "-"<< "\t";
			}
			cout << endl;
	}
	cout << endl;
	std::cout << "Trace ArcMatches " << endl;
	for(size_type i=0;i<=seqA.length()+1;i++){
			for(size_type j=0;j<=seqB.length()+1;j++){
				if(Trace(i,j).arc_match_idx!=0) cout << *(Trace(i,j).arc_match_idx) << "\t";
				else cout << "-"<< "\t";
			}
			cout << endl;
	}
	cout << endl;
}

void ExactMatcher::output_arc_match_score(){
	for(size_type i=0;i<arc_match_score.size();i++){
		if(!(arc_match_score.at(i)==infty_score_t::neg_infty)){
			const ArcMatch &am = arc_matches.arcmatch(i);
			std::cout << i << "):" << "(" << am.arcA().left() << "," << am.arcA().right() << "),("
			<< am.arcB().left() << "," << am.arcB().right() << ")" <<  ": " <<
			score_for_arc_match(am,true) << "," <<
			mappingA.basepair_external(am.arcA().left(),am.arcA().right()) << "," <<
			mappingB.basepair_external(am.arcB().left(),am.arcB().right()) << std::endl;
		}
	}
}




void Mapping::compute_mapping(){
	for(size_type k=0;k<bps.num_bps();k++){
		//for each base pair
		Arc arc = bps.arc(k);
		pos_vec pos_vec_;
		pos_vec new_pos_vec_;
		//add left endpoint of the arc
		pos_vec_.push_back(0);
		new_pos_vec_.push_back(0);
		//compute mapping
		size_type i=1;
		for(size_type j=arc.left()+1;j<arc.right();j++){
			if(is_valid_pos(arc,j)){
				pos_vec_.push_back(j-arc.left());
				new_pos_vec_.push_back(i);
				i++;
			}
			else{
				new_pos_vec_.push_back(-1);
			}
		}
		//add right endpoint of the arc
		pos_vec_.push_back(arc.right()-arc.left());
		new_pos_vec_.push_back(i);

		pos_vecs.push_back(pos_vec_);
		new_pos_vecs.push_back(new_pos_vec_);

	}
}


void Mapping::print_vec() const{
	for(size_type i=0;i<pos_vecs.size();i++){
			std::cout << "i " << i << std::endl;
			for(size_type j=0;j<pos_vecs[i].size();j++){
				std::cout << pos_vecs[i].at(j) << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "new pos vec " << std::endl;
		for(size_type i=0;i<new_pos_vecs.size();i++){
				std::cout << "i " << i << std::endl;
				for(size_type j=0;j<new_pos_vecs[i].size();j++){
					std::cout << new_pos_vecs[i].at(j) << " ";
				}
				std::cout << std::endl;
		}
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

string& PatternPair::get_struct()
{
  return structure;
};

//--------------------------------------------------------------------------
// class PatternPairMap
//    is able to manage a set of PatternPairs(EPMs), each with 2 SinglePatterns
//--------------------------------------------------------------------------
PatternPairMap::PatternPairMap()
{
   idMap.clear();
   patternList.clear();
   patternOrderedMap.clear();
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
   for(patListITER i = patternList.begin();i!=patternList.end();i++)
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
      for (orderedMapITER i=patternOrderedMap.begin();i!=patternOrderedMap.end();i++)
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
   for(patListITER i = patternList.begin();i!=patternList.end();i++)
   {
      bases += (*i)->getSize();
   }
   return bases;
}



LCSEPM::~LCSEPM()
{
	//cout << endl << " execute destructor..." << endl;

	EPM_Table2.clear();
	holeOrdering2.clear();
}

void    LCSEPM::calculateLCSEPM()
{
	cout << " LCSEPM preprocessing..."  <<endl;
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
				for (vector<PatternPairMap::SelfValuePTR>::iterator myIter = EPM_list.begin(); myIter < EPM_list.end(); myIter++){

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
	PatternPairMap::SelfValuePTR lastEPM 	= NULL;
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
			lastEPM = (*t).second;
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
				for (vector<PatternPairMap::SelfValuePTR>::iterator myIter = EPM_list.begin(); myIter < EPM_list.end(); myIter++){
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
  for (PatternPairMap::patListCITER i=myMap.getList().begin();i != myMap.getList().end();i++)
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
	x= 'X';
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

   for (PatternPairMap::patListCITER i=myMap.getList().begin();i != myMap.getList().end();i++)
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

	for (PatternPairMap::patListCITER i=matchedEPMs.getList().begin();i != matchedEPMs.getList().end();i++)
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
	for (int i=last_edge_seq1+1;i<=seqA.length() && seq1_1.size()<seqA.length();++i)
	{
		seq1_1.push_back('.');
		seq1_2.push_back('.');
		seq1_3.push_back('.');
	}

	for (int j=last_edge_seq2+1;j<=seqB.length() && seq2_1.size()<seqB.length();++j)
	{
		seq2_1.push_back('.');
		seq2_2.push_back('.');
		seq2_3.push_back('.');
	}

	seq1_1 += "#";
	seq1_2 += "#";
	seq1_3 += "#";

	seq2_1 += "#";
	seq2_2 += "#";
	seq2_3 += "#";

	outLocARNAfile << ">"<< seqA.names()[0] << endl << upperCase(sequenceA) << endl;
	outLocARNAfile << seq1_1 << endl << seq1_2 << endl << seq1_3 << endl;
	outLocARNAfile << ">"<<seqB.names()[0] << endl << upperCase(sequenceB) << endl;
	outLocARNAfile << seq2_1 << endl << seq2_2 << endl << seq2_3 << endl << endl;

	outLocARNAfile.close();


}


} //end namespace
