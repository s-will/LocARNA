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
			   const string& sequenceA_,
			   const string& sequenceB_,
			   const string& file1_,
			   const string& file2_
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
      sequenceA(sequenceA_),
      sequenceB(sequenceB_),
      file1(file1_),
      file2(file2_)
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
	myLCSEPM= PatternPairMap();
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
	//for(int i=0;i<EPM_start_pos.size();i++){
		//cout << "pos " << EPM_start_pos.at(i).first << "," << EPM_start_pos.at(i).second << endl;
	//}
	cout << "number of starting pos " << EPM_start_pos.size() << endl;
	//cout << "suboptimal score " << subopt_score << endl;
	/*for(int k=0;k<EPM_start_pos.size();k++){
		int i=EPM_start_pos.at(k).first;
		int j=EPM_start_pos.at(k).second;
		trace_in_F_suboptimal(i,j);
	}*/
	trace_in_F_suboptimal(226,363);
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
	cout << "trace in F subopt " << i << "," << j << endl;
	//cout << "score " << F(i,j) << endl;
	//cout << "f" << endl;
	//cout << F << endl;
	static int count=0;
	list<EPM> epms_to_proc_AM;
	list<EPM> epms_to_proc;
	EPM trace_F_cur;
	vector<list<EPM>::iterator> epms_to_proc_delete;
	epm.set_score(0);
	epm.set_curPos(pair<int,int>(i,j));
	epms_to_proc.push_back(epm);
	int cur_i,cur_j,cur_score;
	while(epms_to_proc.begin()!=epms_to_proc.end()){
		//cout << "while " << endl;
		epm = epms_to_proc.front();
		cur_i = epms_to_proc.front().get_curPos().first;
		cur_j = epms_to_proc.front().get_curPos().second;
		cur_score = epms_to_proc.front().get_score();
		//cout << "hier hier " << endl;
		//cout << "cur epm " << endl;
		//epm.print_epm(cout,cur_score);
		if(!(F(cur_i,cur_j).finite_value()==0)){
			//cout << "F matrix " << endl;
			EPM trace_F=epm;
			//cout << "trace F " << endl;
			//trace_F.print_epm(cout,0);
			cout << "cur_i " << cur_i << endl;
			cout << "cur_j " << cur_j << endl;
			if(cur_score+alpha_1*100+F(cur_i-1,cur_j-1).finite_value()>subopt_score
			   && valid_external_pos(cur_i,cur_j)){
				cout << "sequential matching " << endl;
				epm.add(cur_i,cur_j,'.');
				//struct info_for_trace_F tmp = {epm,cur_score+alpha_1*100,pair<int,int>(cur_i-1,cur_j-1)};
				epm.set_score(cur_score+alpha_1*100);
				epm.set_curPos(pair<int,int>(cur_i-1,cur_j-1));
				epms_to_proc.push_back(epm);
				//this->print_epms_to_proc_AGB(epms_to_proc);
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
					//cout << "cur score " << cur_score << endl;
					//cout << "score for arcmatch " << score_for_arc_match(am,true) << endl;
					///cout << "F " << F(am.arcA().left()-1,am.arcB().left()-1).finite_value() << endl;

					//cout << "consider am " << endl;
					//if(am.idx()==607){

					//EPM trace_F =epm;
					//epm.reset();
					//epm.set_score(epms_to_proc.front().get_score()+score_for_arc_match(am,true).finite_value());
					//epm.add_arcmatch(am);
					//list<EPM> epms_to_proc_AM;
					//trace_AGB_suboptimal(am,epm,epm,epms_to_proc,epms_to_proc_AM);
					//copy_epm_trace(epms_to_proc_AM,trace_F,epms_to_proc);
					//epms_to_proc_AM.clear();
					EPM initial_epm = trace_F;
					//cout << "trae F " << endl;
					//initial_epm.print_epm(cout,0);
					//initial_epm.set_arcmatches_to_do(initial_epm.get_arcmatches_to_do().push_back(am.idx()));
					//cout << "store am " << am.idx() << endl;
					initial_epm.store_arcmatch(am.idx());
					initial_epm.set_score(0);
					//initial_epm.set_arcmatches_to_do(initial_epm.get_arcmatches_to_do().push_back(am.idx()));
					epms_to_proc.push_front(initial_epm);
					//cout <<" epms to proc hier " << endl;
					//this->print_epms_to_proc_AGB(epms_to_proc);
					bool finished = false;
					//static int count=0;
					//count++;
					//if(count>100) return;
					cout << "consider am arcA " << am.arcA().left() << "," << am.arcA().right() << endl;
					cout << "consider am arcB " << am.arcB().left() << "," << am.arcB().right() << endl;
					while(!finished){
						//static int count2=0;
						//count2++;
						//if(count2>100) return;
						finished = true;
						for(list<EPM>::iterator it=epms_to_proc.begin();it!=epms_to_proc.end();it++){
							//cout << "iterate over epms to proc " << endl;
							if(it->get_arcmatches_to_do().size()!=0){
								cout << "still am to do " << endl;
								it->print_epm(cout,it->get_score());
								//cout << "iterate over epms to proc " << endl;
								//this->print_epms_to_proc_AGB(epms_to_proc);
								//cout << "trace " << endl;
								//it->print_epm(cout,0);
								//this->print_arcmatches_to_do(it->get_arcmatches_to_do());
								//cout << "cur pos trace " << it->get_curPos().first << "," << it->get_curPos().second << endl;
								//static int count3=0;
								//count3++;
								//if(count3>20){cout << "hi34 " << endl; return;}
								//cout << "count " <<count3 << endl;
								trace_F_cur = *it;

								//cout << "trace_F_cur" << endl;
								//trace_F_cur.print_epm(cout,0);
								//cout << "trace F pos " << trace_F_cur.get_curPos().first << "," << trace_F_cur.get_curPos().second << endl;
								//cout << "tracing " << endl;
								//it->print_epm(cout,it->get_score());
								//cout << it->get_curPos().first << "," << it->get_curPos().second << endl;
								ArcMatch::idx_type next_am_idx = it->next_arcmatch();
								//cout << "trace 2" << endl;
								//it->print_epm(cout,0);
								//this->print_arcmatches_to_do(it->get_arcmatches_to_do());
								const ArcMatch &inner_am = arc_matches.arcmatch(next_am_idx);
								//cout << "TRAC AGB SUBOPTIMAL with " << endl;
								epm.reset();
								epm.set_score(cur_score+score_for_arc_match(am,true).finite_value()
										+F(am.arcA().left()-1,am.arcB().left()-1).finite_value());
								epm.add_arcmatch(inner_am);
								//cout << "before trace AGB subopt " << endl;
								cout << "inner arc A " << inner_am.arcA().left() << "," << inner_am.arcA().right() << endl;
								cout << "inner arc B " << inner_am.arcB().left() << "," << inner_am.arcB().right() << endl;
								//cout << "score " << epm.get_score() << endl;
								//cout << "score vorher " << epms_to_proc.front().get_score() << endl;
								this->trace_AGB_suboptimal(inner_am,epm,epm,epms_to_proc_AM);
								//cout << "after trace AGB subopt " << endl;
								//this->print_epms_to_proc_AGB(epms_to_proc_AM);
								///cout << "arc B " << inner_am.arcB().left() << ","<< inner_am.arcB().right() << endl;
								//cout << "hier " << endl;
								if(epms_to_proc_AM.size()>1){
									//cout << "hier " << endl;
									//cout << "in trace F " << endl;
									//this->print_epms_to_proc_AGB(epms_to_proc_AM);
								}
								//this->print_epms_to_proc_AGB(epms_to_proc_AM);
								//bool traced_arcmatch = false;
								//cout << "trace F cur " << endl;
								//trace_F_cur.print_epm(cout,0);
								//cout << "hier hier " << endl;
								//cout << trace_F_cur.get_curPos().first <<"," << trace_F_cur.get_curPos().second << endl;
								//cout << "inner am " << inner_am.idx() << endl;
								//cout << "inner am " << inner_am.arcA().left() << "," << inner_am.arcB().left() << endl;
								/*for(int i=0;i<trace_F_cur.getStructure().size();i++){
									if(trace_F_cur.getStructure().at(i)!='.'){
										//already traced an arcmatch
										//cout << "already traced arcmatch " << endl;
										//cout << "structure " << trace_F_cur.getStructure() << endl;
										traced_arcmatch = true;
										break;
									}
								}*/
								//if(!traced_arcmatch){
								if(inner_am.arcA().left() < trace_F_cur.get_curPos().first &&
									inner_am.arcB().left() < trace_F_cur.get_curPos().second){
									//cout << "new pos " << inner_am.arcA().left()-1 << "," << inner_am.arcB().left()-1 << endl;
									trace_F_cur.set_curPos(pair<int,int>(inner_am.arcA().left()-1,inner_am.arcB().left()-1));
									trace_F_cur.set_score_to_substract_in_F(F(trace_F_cur.get_curPos().first,trace_F_cur.get_curPos().second).finite_value());

									//cout << "new score substract " << trace_F_cur.get_score_to_substract_in_F() << endl;
									//cout << "new pos " <<  trace_F_cur.get_curPos().first << "," << trace_F_cur.get_curPos().second << endl;
								}
								//cout << "cur pos F " << trace_F_cur.get_curPos().first << "," << trace_F_cur.get_curPos().second << endl;
								copy_epm_trace(epms_to_proc_AM,trace_F_cur,epms_to_proc);
								//cout << "after epm trace " << endl;
								//this->print_epms_to_proc_AGB(epms_to_proc);
								epms_to_proc_AM.clear();
								epms_to_proc_delete.push_back(it);
								finished=false;
							}
						}
						//cout << "delet hier " << endl;
						//this->print_epms_to_proc_AGB(epms_to_proc);
						for(int i=0;i<epms_to_proc_delete.size();i++){
							epms_to_proc.erase(epms_to_proc_delete.at(i));
						}
						epms_to_proc_delete.clear();
						//cout << "after deleting " << endl;
						//this->print_epms_to_proc_AGB(epms_to_proc);
						//if(!test){ finished=true; cout << "finished " << endl;}

					} //while !finished
					cout << "finished " << endl;
				} //valid arcmatch
			} //for iterate over arcmatches
			cout << "done iterating " << endl;
			for(list<EPM>::iterator it=epms_to_proc.begin();it!=epms_to_proc.end();it++){
				if(it->get_score_to_substract_in_F()!=0 &&it->get_arcmatches_to_do().size()==0){
					//cout << "HIER " << it->get_score_to_substract_in_F() << endl;
					//it->print_epm(cout,it->get_score());
					it->set_score(it->get_score()-it->get_score_to_substract_in_F());
					//cout << "new score " << it->get_score() << endl;
					it->set_score_to_substract_in_F(0);
				}
			}
		}
		else{
			//cout << "epm found " << endl;
			//epm is computed
			count++;
			//epm.print_epm(cout,0);
			epm.sort_patVec();
			cout << "final epm " << endl;
			epm.print_epm(cout,cur_score);
			if(!validate_epm()) { cout << "falsches EPM "<< endl; return;}
		}
		//cout << "pop front " << endl;
		epms_to_proc.pop_front();
		//cout << "after pop front " << endl;
		//this->print_epms_to_proc_AGB(epms_to_proc);
	}
	epm.reset();
	cout << "number of EPM " << count << endl;
}

//TODO stacking probabilities
void
ExactMatcher::trace_AGB_suboptimal(const ArcMatch &am, EPM &epm_to_store,EPM cur_epm_AM, list<EPM> &epms_to_proc_AM){
	compute_AGBmatrices(am);
	if(am.arcA().left()==4 && am.arcA().right()==17 && am.arcB().left()==4 && am.arcB().right()==17){
		//cout << " trace AGB suboptimal " << endl;
		//cout << " arc A " << am.arcA().left() << "," << am.arcA().right() << endl;
		//cout << " arc B " << am.arcB().left() << "," << am.arcB().right() << endl;
		//this->print_matrices(am,20,20);
	}
	//cout << " arc A " << am.arcA().left() << "," << am.arcA().right() << endl;
	//cout << " arc B " << am.arcB().left() << "," << am.arcB().right() << endl;
	//this->print_matrices(am,20,20);
	size_type number_of_posA = mappingA.number_of_valid_pos(am.arcA().idx());
	size_type number_of_posB = mappingB.number_of_valid_pos(am.arcB().idx());
	pair<int,int> str_pos;
	pair<int,int> curPos;
	int i,j;
	size_type posA,posB;
	//EPM cur_epm_AM;
	int state=in_B;
	int cur_score;
	//bool seq_matching;
	//Initialization
	//cur_epm_AM = epm_to_proc_AM;//epms_to_proc_AM.front();
	//state = cur_epm_AM.get_state();
	cur_score = cur_epm_AM.get_score() - arc_match_score.at(am.idx()).finite_value();
	if(am.arcA().left()==1 && am.arcA().right()==16 && am.arcB().left()==1 && am.arcB().right()==16 ||
	am.arcA().left()==5 && am.arcA().right()==14 && am.arcB().left()==5 && am.arcB().right()==11)
	{
		//cout << "cur score hier " << cur_score << endl;
		//cur_epm_AM.print_epm(cout,cur_score);
	}
	//cout << "cur score " << cur_score << endl;
	//cout << "arc match score under " << arc_match_score.at(am.idx()).finite_value() << endl;
	//cout << "cur score " << cur_score << endl;
	//pair<int,int> cur_pos_seq = pair<int,int>(cur_epm_AM.getPat1Vec().back(),cur_epm_AM.getPat2Vec().back());
	curPos=pair<int,int>(mappingA.number_of_valid_pos(am.arcA().idx())-1,mappingB.number_of_valid_pos(am.arcB().idx())-1);
	bool not_finished=true;
	//list<EPM> epms_to_proc_AM;
	epms_to_proc_AM.push_back(cur_epm_AM);
	//cout << "curPos " << curPos.first << "," << curPos.second << endl;
	//int count =0;
	while(not_finished){
		//static int count=0;
		//count++;
		//if(count>100) return;
		//cout << "epms to proc AM " << endl;
		//this->print_epms_to_proc_AGB(epms_to_proc_AM);
		not_finished=false;
	while((curPos!=pair<int,int>(0,0) || state!=in_A)){
		//static int count=0;
		//count++;
		//if(count>100) return;
		//cout << "cur epm " << endl;
		//epm.print_epm(cout,0);
		//this->print_arcmatches_to_do(epm.get_arcmatches_to_do());
		//cout << "trace AGB sub epms to store AM " << endl;
		//this->print_epms_to_proc_AGB(epms_to_proc_AM);
		//count++;
		//if(count>20) return;
		//cout << "print epms to proc AM " << endl;
		//this->print_epms_to_proc(epms_to_proc_AM);
		//count++;
		i=curPos.first;
		j=curPos.second;
		posA = mappingA.get_pos_in_seq(am.arcA(),i);
		posB = mappingB.get_pos_in_seq(am.arcB(),j);
		//cout << "posA " << posA << endl;
		//cout << "posB " << posB << endl;
		//cout << "bool " << (seqA[posA]==seqB[posB]) << endl;
		//cout << " epms to proc AM " << endl;
		//this->print_epms_to_proc_AGB(epms_to_proc_AM);
		//cout << " cur pos begin " << curPos.first << "," << curPos.second << endl;
		//cout << " state " << state << endl;
		//cur_epm_AM.print_epm(cout,-3);
		//this->print_arcmatches_to_do(cur_epm_AM.get_arcmatches_to_do());
		//seq_matching = false;
		cur_epm_AM.set_curPos(curPos);
		epm = cur_epm_AM;
		//epm.print_epm(cout,epm.get_score());
		//cout << "pos epm "<<epm.get_curPos().first << "," << epm.get_curPos().second << endl;

		//cout << "seq matching " << seq_matching << endl;
		switch(state){
		//static int count=0;
				//count++;
				//if(count>100) return;
			case in_B:
			{
				//cur_epm_AM.print_epm(cout,-8);
				if(seqA[posA]==seqB[posB]){
					//structural matching
					str_trace_AGB_suboptimal(B,posA,posB,am,cur_score,state,curPos,epms_to_proc_AM,cur_epm_AM);
					//sequential matching
					//cout << "i-1 " << i-1 << endl;
					//cout << "j-1 " << j-1 << endl;
					//cout << "B " << B(i-1,j-1) << endl;
					//cout << "seq matching " << seq_matching << endl;
					if(i>=1&&j>=1){
					bool flag = (subopt_score < cur_score+alpha_1*100+B(i-1,j-1).finite_value());
					//cout << "subopt score " << subopt_score << endl;
					//cout << "cur score " << cur_score << endl;
					//cout << "B " << B(i-1,j-1) << endl;
					//cout << "falg " << flag << endl;
					//cout << "Mapping A " << mappingA.seq_matching(am.arcA().idx(),i) << endl;
					//cout << "MappingB " << mappingB.seq_matching(am.arcB().idx(),j) << endl;
					if(flag && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						//cout << " sequential matching B " << endl;
						int add_score=0;
						if(posA!=am.arcA().right() && posB!=am.arcB().right()){
							epm.add(posA,posB,'.');
							add_score = alpha_1*100;
						}
						//c//out << "Hi§$" << endl;
							//cur_epm_AM.print_epm(cout,-10);
							//this->print_arcmatches_to_do(cur_epm_AM.get_arcmatches_to_do());
							//epm.print_epm(cout,-5);
							//this->print_arcmatches_to_do(epm.get_arcmatches_to_do());
							curPos=pair<int,int>(i-1,j-1);
							epm.set_score(cur_score+add_score);
							epm.set_state(state);
							epm.set_curPos(curPos);
							epms_to_proc_AM.push_back(epm);
							//epm.print_epm(cout,epm.get_score());
							epm = cur_epm_AM;
							//epm.print_epm(cout,-5);
							//cur_epm_AM.print_epm(cout,-6);
							//cout << "hier " << endl;
							//cout << "i,j " << i << "," << j << endl;
							//cout << "state " << state << endl;
							//seq_matching = true;
					}
					}
					//cout << "seq matching " << seq_matching << endl;
				}
				//cout << "i,j " << i << ", " << j << endl;
				//cout << "seq_matching " << seq_matching << endl;
				//cout << "B " << B(i,j) << endl;
				//cout << "G " << G(i,j) << endl;
				if(B(i,j)==G(i,j)){
					//epm.print_epm(cout,-2);
					//cout << "curPos epm " <<epm.get_curPos().first << "," << epm.get_curPos().second << endl;
					//cout << " B -> G " << endl;
					epm.set_state(in_G);
					epms_to_proc_AM.push_back(epm);
					//cout << "new epm " << endl;
					//epm.print_epm(cout,epm.get_score());
					//state=in_G;
					epm = cur_epm_AM;
					break;
				}
				bool flag2 = B(i,j)==G(i-1,j-1)+alpha_1*100 && !(mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j));
				//special case
				if(flag2){
					//cout << " special case " << endl;
					curPos=pair<int,int>(i-1,j-1);
					if(posA!=am.arcA().right() && posB!=am.arcB().right()){
						epm.add(posA,posB,'.');
					}
					//cout << "special case " << endl;
					epm.set_state(in_G);
					epm.set_curPos(curPos);
					epms_to_proc_AM.push_back(epm);
					//epm.print_epm(cout,epm.get_score());
					//state=in_G;
					epm = cur_epm_AM;
					break;
				}
			//this->print_epms_to_proc_AGB(epms_to_proc_AM);
			//epm.print_epm(cout,epm.get_score());
			}
			break;
			case in_G:
			{
				if(G(i,j)==A(i,j)){
					//cout << " G->A "<< endl;
					epm.set_score(cur_score);
					epm.set_state(in_A);
					epm.set_curPos(curPos);
					epms_to_proc_AM.push_back(epm);
					epm = cur_epm_AM;
					state=in_A; break;
				}
				if(i>=1 && j>=1 && G(i,j)==G(i-1,j-1)){
					//cout << " G diagonal " << endl;
					curPos=pair<int,int>(i-1,j-1);
					break;
				}
				if(i>=1 && G(i,j)==G(i-1,j)){
					//cout << " G up " << endl;
					curPos=pair<int,int>(i-1,j);
					break;
				}
				if(j>=1 && G(i,j)==G(i,j-1)){
					//cout << " G left " << endl;
					curPos=pair<int,int>(i,j-1);
					break;
				}
			}
			case in_A:
			{
				if(seqA[posA]==seqB[posB]){
					//structural matching
					str_trace_AGB_suboptimal(A,posA,posB,am,cur_score,state,curPos,epms_to_proc_AM,cur_epm_AM);
					//cout << "after structural matching " << endl;
					//sequential matching
					if(i>=1&&j>=1){
					bool flag3 = (subopt_score < cur_score+alpha_1*100+A(i-1,j-1).finite_value());
					//cout << "i,j " << i << "," << j << endl;
					//cout << "flag3 " << (subopt_score < cur_score+alpha_1*100+A(i-1,j-1).finite_value())<<endl;;
					//cout << "Mapping A " << mappingA.seq_matching(am.arcA().idx(),i) << endl;
					//cout << "gemappte pos " << mappingA.get_pos_in_seq(am.arcA(),i) << endl;
					//cout << "Mapping B " << mappingB.seq_matching(am.arcB().idx(),j) << endl;
					//cout << "gemappte pos " << mappingB.get_pos_in_seq(am.arcB(),j) << endl;
					if(flag3 && mappingA.seq_matching(am.arcA().idx(),i) && mappingB.seq_matching(am.arcB().idx(),j)){
						//cout << " sequential matching A" << endl;
						curPos=pair<int,int>(i-1,j-1);
						if(posA!=am.arcA().right() && posB!=am.arcB().right()){
							epm.add(posA,posB,'.');
							epm.set_score(cur_score+alpha_1*100);
							epm.set_state(state);
							epm.set_curPos(curPos);
							epms_to_proc_AM.push_back(epm);
							//cout << "sequential matching A " << endl;
							//TODO reset cur epm
							epm = cur_epm_AM;
						}
						break;
					}
					}
				}
			}
		}
		if(state!=in_G){
			//cout << " pop fonrt AM " << endl;
			if(epms_to_proc_AM.size()>1){
				//cout << "pop front " << endl;
				epms_to_proc_AM.pop_front();
				cur_epm_AM = epms_to_proc_AM.front();
				//this->print_epms_to_proc_AGB(epms_to_proc_AM);
				//cout << endl;
			}
			else{
				cout << "break " << endl;
				this->print_epms_to_proc_AGB(epms_to_proc_AM);
				break;
			}
			state = cur_epm_AM.get_state();
			cur_score = cur_epm_AM.get_score();
			//pair<int,int> cur_pos_seq = pair<int,int>(cur_epm_AM.getPat1Vec().back(),cur_epm_AM.getPat2Vec().back());
			curPos=cur_epm_AM.get_curPos();
			//cout << "cur pos hier " << curPos.first << "," << curPos.second << endl;
			//epm = cur_epm_AM;
			//epm.print_epm(cout,4);
		}
		//cout << "curPos " << curPos.first << "," << curPos.second << endl;
		//cout << "state " << state << endl;
		}//first while
		//cout << "store  epm " << endl;
		//epms_to_proc.push_back(epm); //store only part under the am
		//cout << "store epm in epms to proc " << endl;
		//epm.print_epm(cout,0);
		//epms_to_proc_AM.pop_front(); //remove element from epms_to_proc_AM
		//if(epms_to_proc_AM.size()!=0){
			//cur_epm_AM = epms_to_proc_AM.front();
		//}
		//check if finished
		for(list<EPM>::iterator it = epms_to_proc_AM.begin();it!=epms_to_proc_AM.end();it++){
			if(it->get_curPos()!=pair<int,int>(0,0)){
				//static int count2=0;
				//count2++;
				//if(count2>10){cout << "abbruch " << endl; return;}
				//cout << " not finished yet " << endl;
				//it->print_epm(cout,it->get_score());
				not_finished=true;
				cur_epm_AM = *it;
				//Initialize
				state = cur_epm_AM.get_state();
				cur_score = cur_epm_AM.get_score();
				//pair<int,int> cur_pos_seq = pair<int,int>(cur_epm_AM.getPat1Vec().back(),cur_epm_AM.getPat2Vec().back());
				curPos=cur_epm_AM.get_curPos();
				epm = cur_epm_AM;
				break;
			}
		}
		//if(not_finished){this->print_epms_to_proc_AGB(epms_to_proc_AM);}
	}//second while
	if(am.arcA().left()==4 && am.arcA().right()==17 && am.arcB().left()==4 && am.arcB().right()==17){
		cout << " finished!" << endl;
		return;

	//	cout << "epms to proc for am (4,17),(4,17) " << endl;
	//	this->print_epms_to_proc_AGB(epms_to_proc_AM);
	}
}

void
ExactMatcher::copy_epm_trace(list<EPM> &epms_to_proc_AM,EPM &trace_F, list<EPM> &epms_to_proc){
	//cout << "copy epm trace " << endl;
	//cout << "epms to proc AM " << endl;
	//this->print_epms_to_proc_AGB(epms_to_proc_AM);
	//cout << "epms to proc " << endl;
	//this->print_epms_to_proc_AGB(epms_to_proc);
	//cout << "trace f " << endl;
	//trace_F.print_epm(cout,0);
	//cout << trace_F.get_curPos().first << ","<< trace_F.get_curPos().second << endl;
	EPM combined_epm;
	for(list<EPM>::iterator it = epms_to_proc_AM.begin();it!=epms_to_proc_AM.end();it++){
		combined_epm = trace_F;
		//copy information from trace AGB
		for(int i=0;i<it->getPat1Vec().size();i++){
			combined_epm.add(it->getPat1Vec().at(i),it->getPat2Vec().at(i),it->getStructure().at(i));
		}
		//copy arcmatches_to_do
		//it->get_arcmatches_to_do().pop_back();
		combined_epm.set_arcmatches_to_do(it->get_arcmatches_to_do());
		//cout << "combined epm " << endl;
		//combined_epm.print_epm(cout,0);
		//cout << "curPos " << combined_epm.get_curPos().first << "," << combined_epm.get_curPos().second <<endl;
		combined_epm.set_score(it->get_score());
		//combined_epm.set_score_to_substract_in_F(0);
		//cout << "score hier " << combined_epm.get_score() << endl;
		//store information in epms_to_proc
		//cout << "combined epm " << endl;
		//combined_epm.print_epm(cout,combined_epm.get_score());
		epms_to_proc.push_back(combined_epm);
	}
	//cout << "after copy epm trace " << endl;
	//this->print_epms_to_proc_AGB(epms_to_proc);
}

/*void
ExactMatcher::trace_AGB_suboptimal_main(const ArcMatch &am, EPM &epm_to_store,list<EPM> &epms_to_proc_AGB,list<EPM> &epms_to_proc){
	//while(epms_to_proc_AGB.size()!=0){
		list<EPM> tmp;
		tmp.push_back(epms_to_proc_AGB.front());
		trace_AGB_suboptimal(am,epm_to_store,tmp,epms_to_proc_AGB, epms_to_proc);
		epms_to_proc_AGB.pop_front();
	//}
}*/


void
ExactMatcher::str_trace_AGB_suboptimal(ScoreMatrix &mat,int posA, int posB,const ArcMatch &am,int cur_score,int state, pair<int,int> curPos,list<EPM> &epms_to_proc_AM,EPM cur_epm_AM){
	//cout << "posA-1 " << posA-1 << endl;
	//cout << "posB-1 " << posB-1 << endl;
	for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(posA-1,posB-1).begin();
			arc_matches.common_right_end_list(posA-1,posB-1).end() != it; ++it ){
		const ArcMatch &inner_am = arc_matches.arcmatch(*it);
		if(arc_match_score.at(inner_am.idx())==infty_score_t::neg_infty){
			continue;
		}
		if(inner_am.arcA().left()>am.arcA().left() && inner_am.arcB().left()>am.arcB().left()){
			//cout << "am Arc A " << inner_am.arcA().left() << "," << inner_am.arcA().right() << endl;
			//cout << "am Arc B " << inner_am.arcB().left() << "," << inner_am.arcB().right() << endl;
			int pos_before_arcA = mappingA.get_pos_in_new_seq(am.arcA(),inner_am.arcA().left()-1);
			int pos_before_arcB= mappingB.get_pos_in_new_seq(am.arcB(),inner_am.arcB().left()-1);
			//cout << "pos before arcA " << pos_before_arcA << endl;
			//cout << "pos before arcB " << pos_before_arcB << endl;
			if(pos_before_arcA == -1 || pos_before_arcB == -1) continue; //no valid position before the arc
			//cout << "arcA " << inner_am.arcA().left() << "," << inner_am.arcA().right() << endl;
			//cout << "arcB " << inner_am.arcB().left() << "," << inner_am.arcB().right() << endl;
			//cout << "mat " << mat(pos_before_arcA,pos_before_arcB) << endl;
			//cout << "mat " << mat(pos_before_arcA,pos_before_arcB).finite_value() << endl;
			bool flag = (subopt_score < (cur_score + (score_for_arc_match(inner_am,true)+
						score_for_stacking(am,inner_am)+alpha_1*100+
						mat(pos_before_arcA,pos_before_arcB)).finite_value()));
			//cout << "cur score " << cur_score << endl;
			//cout << "score for inner am " << score_for_arc_match(inner_am,true) << endl;
			//cout << "score for stacking " << score_for_stacking(am,inner_am) << endl;
			//cout << "sequential score " << alpha_1*100 << endl;
			//cout << "pos before arc A " << pos_before_arcA << endl;
			//cout << "pos before arc B " << pos_before_arcB << endl;
			//cout << "score before " << mat(pos_before_arcA,pos_before_arcB) << endl;
			//cout << "score before " << A(pos_before_arcA,pos_before_arcB) << endl;
			//cout << "gesamt " << (cur_score + (score_for_arc_match(inner_am,true)+
			//		score_for_stacking(am,inner_am)+alpha_1*100+
			//		mat(pos_before_arcA,pos_before_arcB)).finite_value());
			//cout << "flag " << flag << endl;
			if(flag){
				//cout << "structural matching " << endl;
				if(posA!=am.arcA().right() && posB!=am.arcB().right()){
					//cout << "arcA " << inner_am.arcA().left() << "," << inner_am.arcA().right() << endl;
					//cout << "arcB " << inner_am.arcB().left() << "," << inner_am.arcB().right() << endl;
					//add paired bases to structure and corresponding position
					//int newposA = mappingA.number_of_valid_pos(inner_am.arcA().idx())-1;
					//int newposB = mappingB.number_of_valid_pos(inner_am.arcB().idx())-1;
					//curPos=pair<int,int>(newposA,newposB);
					//cout << "new cur pos " << newposA << "," << newposB << endl;
					epm.add(posA,posB,'.');
				}
				//cout << " str_trace AGB " << endl;
				//cout << " arcA " << inner_am.arcA().left() << "," << inner_am.arcA().right() << endl;
				//cout << " arcB " << inner_am.arcB().left() << "," << inner_am.arcB().right() << endl;
				//epm.add_arcmatch(inner_am);
				//cout << "new score " << cur_score+(score_for_arc_match(inner_am,true)+score_for_stacking(am,inner_am)+alpha_1*100).finite_value() << endl;
				epm.store_arcmatch(inner_am.idx());
				epm.set_score(cur_score+(score_for_arc_match(inner_am,true)+score_for_stacking(am,inner_am)+alpha_1*100).finite_value());
				epm.set_state(state);
				epm.set_curPos(pair<int,int>(pos_before_arcA,pos_before_arcB));
				epms_to_proc_AM.push_back(epm);
				//cout << "new epm " << endl;
				//epm.print_epm(cout,epm.get_score());
				//this->print_arcmatches_to_do(epm.get_arcmatches_to_do());
				epm = cur_epm_AM;
				//epm.print_epm(cout,epm.get_score());
				//this->print_arcmatches_to_do(epm.get_arcmatches_to_do());
				//epm.get_arcmatches_to_do().pop_back();
			}
		}
	}
	//cout << "after str trace subopt " << endl;
}


/*void ExactMatcher::swap(list<EPM> &epms_to_proc_AM,list<EPM>::iterator &next_el){
	//struct info_for_trace_AGB tmp = first_el;
	list<info_for_trace_AGB>::iterator i = epms_to_proc_AM.begin();
	//cout << (*i) << endl;
	//cout << (*j) << endl;
	info_for_trace_AGB temp = *i;
	*i = *next_el;
	*next_el = temp;
	//cout << (*i) << endl;
	//cout << (*j) << endl;
}*/

void ExactMatcher::print_epms_to_proc_AGB(list<EPM> &epms_to_proc_AGB){
	cout << "_________________________________________________" << endl;
	cout << " epms to proc AGB " << endl;
	for(list<EPM>::iterator it = epms_to_proc_AGB.begin();it!=epms_to_proc_AGB.end();it++){
		it->print_epm(cout,it->get_score());
		cout << " curPos " << it->get_curPos().first << "," << it->get_curPos().second << endl;
		cout << " curState " << it->get_state() << endl;
		print_arcmatches_to_do(it->get_arcmatches_to_do());
	}
	cout << endl;
	cout << "_____________________________________________________" << endl;
	cout << endl;
}

void ExactMatcher::print_arcmatches_to_do(std::vector<ArcMatch::idx_type> arcmatches_to_do){
	if(arcmatches_to_do.size()!=0){
		cout << " arcmatches to do " << endl;
		for(std::vector<ArcMatch::idx_type>::iterator it2 = arcmatches_to_do.begin();it2!=arcmatches_to_do.end();it2++){
			cout << " " << *it2 << endl;
			const ArcMatch &am = arc_matches.arcmatch(*it2);
			cout << " arcA " << am.arcA().left() << "," << am.arcA().right() << endl;
			cout << " arcB " << am.arcB().left() << "," << am.arcB().right() << endl;
		}
	}
}

void ExactMatcher::print_epms_to_proc(list<EPM> &epms_to_proc){
	cout << "printing epms to proc " << endl;
	cout << "_____________________________________________________" << endl;
	int count =0;
	for(list<EPM>::iterator it=epms_to_proc.begin();it!=epms_to_proc.end();it++){
		cout << count << ".EPM: ";
		//EPM cur_epm = *it;
		it->print_epm(cout,it->get_score());
		cout << "curPos " << it->get_curPos().first << "," << it->get_curPos().second << endl;
		count++;
	}
	cout << "_____________________________________________________" << endl;
	cout << endl;
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
	const int& size1= (int)seqA.length();
	const int& size2= (int)seqB.length();
	int size= 0;
	int score= 0;
	cout << "#EPM: " << mcsPatterns.size() << endl;
	time_t start_chaining = time (NULL);
	//create LCSEPM object
	//LCSEPM patterns(size1, size2 ,epm.get_patternPairMap(), myLCSEPM,size,score, EPM_min_size);
	LCSEPM patterns(size1, size2 ,mcsPatterns, myLCSEPM,size,score, EPM_min_size);
	//begin chaining algorithm
	patterns.calculateLCSEPM();
	time_t stop_chaining = time (NULL);
    cout << "time for chaining : " << stop_chaining - start_chaining << "sec " << endl;
	//output patterns to PS files
    time_t start_ps = time (NULL);
	patterns.MapToPS(sequenceA, sequenceB, size, myLCSEPM, file1, file2);
	time_t stop_ps = time (NULL);
	cout << "time for map to ps : " << stop_ps - start_ps << "sec " << endl;
}

void
ExactMatcher::get_matching(size_type i, size_type j){
    bool valid=true;
	static int count= 0;
	static int count2=0;
	static int count_matching = 0;
	count_matching++;
	const std::string &seq1_id= "sequence 1";
	const std::string &seq2_id= "sequence 2";
	SinglePattern pattern1,pattern2;
	string patId;
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
		  patId= ss.str();
		  pattern1 = SinglePattern(patId,seq1_id,epm.getPat1Vec());
		  pattern2 = SinglePattern(patId,seq2_id,epm.getPat2Vec());
		  int score= score1.finite_value();
		  //epm.add_pattern(patId,pattern1,pattern2, score);
		  mcsPatterns.add(patId, epm.getPat1Vec().size(), pattern1, pattern2, epm.getStructure(), score);
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

bool ExactMatcher::validate_epm(){
	//EPMs have the same size
	if(epm.getPat1Vec().size()!=epm.getPat2Vec().size()){
		cerr << "wrong pattern " << endl;
		return false;
	}

	//two matched positions in the EPMs have the same nucleotide
	for(int i=0;i<epm.getPat1Vec().size();i++){
		if(seqA[epm.getPat1Vec().at(i)]!=seqB[epm.getPat2Vec().at(i)]){
			cerr << "no EPM " << endl;
			return false;
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
						if(gap){ cerr << "ERROR " << endl; return false;}
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
						if(gap){ cerr << "ERROR " << endl;return false;}
						gap=true;
					}
				}
			}
		}
	}
	return true;
}

void ExactMatcher::output_locarna(){
  
	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	for (PatternPairMap::patListCITER i=myLCSEPM.getList().begin();i != myLCSEPM.getList().end();i++)
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
	string outname = "locarna_constraints_input.txt"; //"/home/radwan/Exparna_P/LocARNA/src/locarna_constraints_input.txt";
	ofstream outLocARNAfile (outname.c_str());

	int last_edge_seq1,last_edge_seq2;
	last_edge_seq1=0;
	last_edge_seq2=0;

	string seq1_1,seq1_2,seq1_3,seq2_1,seq2_2,seq2_3;
	int edge = 100;

	for (vector<intPair>::iterator i_edge = matchingsLCSEPM.begin(); (i_edge != matchingsLCSEPM.end() && seq1_1.size()<sequenceA.length() && seq2_2.size()<sequenceB.length());++i_edge)
	{
		//cout << "first: " << (*i_edge).first << " second: " << (*i_edge).second << endl;

		for (int i=last_edge_seq1+1;(i<(int)((*i_edge).first)&&seq1_1.size()<sequenceA.length());++i)
		{
			seq1_1.push_back('.');
			seq1_2.push_back('.');
			seq1_3.push_back('.');
		}

		for (int j=last_edge_seq2+1;(j<(int)((*i_edge).second)&&seq2_2.size()<sequenceB.length());++j)
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
	for (int i=last_edge_seq1+1;i<=(int)sequenceA.length()&&seq1_1.size()<sequenceA.length();++i)
	{
		seq1_1.push_back('.');
		seq1_2.push_back('.');
		seq1_3.push_back('.');
	}

	for (int j=last_edge_seq2+1;j<=(int)sequenceB.length()&&seq2_1.size()<sequenceB.length();++j)
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

	outLocARNAfile << "> Sequence1" << endl << upperCase(sequenceA) << endl;
	outLocARNAfile << seq1_1 << endl << seq1_2 << endl << seq1_3 << endl;
	outLocARNAfile << "> Sequence2" << endl << upperCase(sequenceB) << endl;
	outLocARNAfile << seq2_1 << endl << seq2_2 << endl << seq2_3 << endl << endl;

	outLocARNAfile.close();
	
	
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
// a pattern consist of elements
//--------------------------------------------------------------------------


SinglePattern::SinglePattern(const string& myId_,const string& seqId_,const intVec& mySinglePattern_)
   :myId(myId_),seqId(seqId_),pattern(mySinglePattern_)
   {
   }

SinglePattern::~SinglePattern()
{
	pattern.clear();
}

const string& SinglePattern::getmyId() const
   {
      return myId;
   }

const string& SinglePattern::getseqId() const
   {
      return seqId;
   }

const intVec& SinglePattern::getPat() const
   {
      return pattern;
   }


//--------------------------------------------------------------------------
// class PatternPair
//    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
//--------------------------------------------------------------------------
const string& PatternPair::getId() const
{
   return id;
};

const int& PatternPair::getSize() const
{
   return size;
};

const SinglePattern& PatternPair::getFirstPat() const
{
   return first;
};

const SinglePattern& PatternPair::getSecPat() const
{
   return second;
};

void PatternPair::resetBounds()
{
	  insideBounds.clear();
}

void PatternPair::setOutsideBounds(intPPair myPPair)
{
	  outsideBounds = myPPair;
};

intPPair PatternPair::getOutsideBounds()
{
	  return outsideBounds;
};

void PatternPair::addInsideBounds(intPPair myPPair)
{
	  insideBounds.push_back(myPPair);
};

const vector<intPPair>& PatternPair::getInsideBounds()
{
	  return insideBounds;
};

void PatternPair::setEPMScore(int myScore)
{
	  score = myScore;
};

string& PatternPair::get_struct()
{
  return structure;
};

int PatternPair::getScore()
      {
    	  return score;
      };

int PatternPair::getEPMScore()
      {
    	  return EPMscore;
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
                         const int& mysize,
                         const SinglePattern& first,
                         const SinglePattern& second,
			 const string& structure,
			 int score
			)
{
   PatternPair* p= new PatternPair(id,mysize,first,second,structure,score);
   SelfValuePTR myP = SelfValuePTR(p);
   patternList.push_back(myP);
   idMap.insert(make_pair(id,myP));
}

void PatternPairMap::add(const SelfValuePTR value)
{
   SelfValuePTR myP = SelfValuePTR(new PatternPair(*value));
   patternList.push_back(myP);
   idMap.insert(make_pair(value->getId(),myP));
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
	EPM_Table.clear();
	for(HoleMapCITER it= holeOrdering.begin(); it!=holeOrdering.end(); it++)
		delete (*it).second;
	holeOrdering.clear();
	
}

void    LCSEPM::calculateLCSEPM()
{
 	preProcessing();
	calculateHoles();
	vector < vector<int> > last_vec;
	int i,k;
	i = 1;
	k = 1;
	LCSEPMscore = D_rec(i,size1,k,size2,last_vec,false);
	cout << "    Score LCS-EPM: "<< LCSEPMscore <<endl;
	calculateTraceback(i,size1,k,size2,last_vec);
	LCSEPMsize = matchedEPMs.getMapBases();
	cout << "    #EPMs: "<< matchedEPMs.size() << " / matched Bases: "<< LCSEPMsize <<endl;
}

void    LCSEPM::calculatePatternBoundaries(PatternPair*   myPair)
{
   const vector<unsigned int>& myPatStr1 = myPair->getFirstPat().getPat();
   const vector<unsigned int>& myPatStr2 = myPair->getSecPat().getPat();

   myPair->resetBounds();

	   for (unsigned int k=1;k < (myPatStr1.size());++k)
	   {
 	   if ((myPatStr1[k]-EPM_min_size > myPatStr1[k-1])
     		  &&(myPatStr2[k]-EPM_min_size > myPatStr2[k-1]))
 		  {
    	  myPair->addInsideBounds(std::make_pair(make_pair(myPatStr1[k-1],myPatStr1[k]),make_pair(myPatStr2[k-1],myPatStr2[k])));
 		  }
	   }
	   // insert global min/max of the pattern
	   myPair->setOutsideBounds(make_pair(make_pair(myPatStr1.front(),myPatStr1.back()),make_pair(myPatStr2.front(),myPatStr2.back())));
 
 }

void LCSEPM::preProcessing()
{
    EPM_Table.resize(size1+1);
    for (unsigned int i = 0; i < EPM_Table.size();++i)
    	EPM_Table[i].resize(size2+1);

    for (PatternPairMap::patListCITER myPair = patterns.getList().begin(); myPair != patterns.getList().end(); ++myPair)
    {
        calculatePatternBoundaries(*myPair);


        EPM_Table[(*myPair)->getOutsideBounds().first.second][(*myPair)->getOutsideBounds().second.second] = (*myPair);

        for(IntPPairCITER h = (*myPair)->getInsideBounds().begin(); h != (*myPair)->getInsideBounds().end(); ++h)
        {
            HoleKeyPTR myHoleKey = HoleKeyPTR(new HoleKey());
            int holeSize = (*h).first.second - (*h).first.first - 1;
            myHoleKey->bounds = (*h);
            myHoleKey->pattern = (*myPair);
            holeOrdering.insert(make_pair(holeSize,myHoleKey));
        }
    }
}


int LCSEPM::max3(int a, int b, int c)
{
	int tmp = a>b? a:b;
	return (tmp>c? tmp:c);
}

int LCSEPM::D_rec(const int& i,const  int& j,const int& k,const int& l,vector < vector<int> >& D_h,const bool debug)
{

	int 					score_EPM;
	int 					pos_before_EPM_Str1;
	int 					pos_before_EPM_Str2;

	D_h.clear();
	D_h.resize(j - i + 2);
	for (unsigned int a = 0; a < D_h.size();++a)
		D_h[a].resize(l - k + 2,0);
	
	for(unsigned int j_1 = 1; j_1 < (j-i+2); ++j_1)
		for (unsigned int l_2 = 1; l_2 < (l-k+2); ++l_2)
		{
			if (EPM_Table[i + j_1-1][k + l_2-1] == NULL)
			{
				D_h[j_1][l_2] = (D_h[j_1-1][l_2]>D_h[j_1][l_2-1])? D_h[j_1-1][l_2]:D_h[j_1][l_2-1] ;
			}
			else
			{
				pos_before_EPM_Str1 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().first.first ) - i;
				pos_before_EPM_Str2 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().second.first ) - k;
				if ((pos_before_EPM_Str1 < 0)||(pos_before_EPM_Str2 <0))
					score_EPM = 0;
				else
					score_EPM = D_h[pos_before_EPM_Str1][pos_before_EPM_Str2] + EPM_Table[i + j_1-1][k + l_2-1]->getScore();
				D_h[j_1][l_2] = max3(score_EPM,D_h[j_1-1][l_2],D_h[j_1][l_2-1]);
			}
		}
	return (D_h[j - i + 1][l - k + 1]);
}

void LCSEPM::calculateHoles()
{
	vector < vector<int> > vec;
	for (HoleMapCITER t = holeOrdering.begin();t != holeOrdering.end();++t)
    {
		bool deb=false;
		(*t).second->pattern->setEPMScore(	(*t).second->pattern->getScore() +
											D_rec((*t).second->bounds.first.first+1,(*t).second->bounds.first.second-1,\
											(*t).second->bounds.second.first+1,(*t).second->bounds.second.second-1,vec,deb));
    }
}

void LCSEPM::calculateTraceback(const int i,const  int j,const int k,const int l,vector < vector<int> > holeVec)
{
	int j_1 = holeVec.size()-1;
	int l_2 = holeVec[0].size()-1;
	while ((j_1 >= 1)&&(l_2 >= 1))
	{
		if (holeVec[j_1 - 1][l_2] == holeVec[j_1][l_2])
			--j_1;
		else
			if (holeVec[j_1][l_2 - 1] == holeVec[j_1][l_2])
				--l_2;
			else
			{
				vector < vector<int> > tmpHoleVec;
				matchedEPMs.add(EPM_Table[i + j_1 - 1][k + l_2 - 1]);
				for(IntPPairCITER h = EPM_Table[i + j_1 - 1][k + l_2 - 1]->getInsideBounds().begin(); h != EPM_Table[i + j_1-1][k + l_2-1]->getInsideBounds().end(); ++h)
				{
					int sc = D_rec((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec,false);
					if (sc != 0)
						calculateTraceback((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec);
				}
				int s1 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().first.first ) - i;
				int s2 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().second.first) - k;
				j_1 = s1;
				l_2 = s2;
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
void LCSEPM::MapToPS(const string& sequenceA, const string& sequenceB, int& mySize, PatternPairMap& myMap, const string& file1, const string& file2)
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

   
   string psfilename = file1+".ps";
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
   
   psfilename = file2+".ps";
   free(structure);
   structure= NULL;
   if(base_pair){ free(base_pair); base_pair= NULL;}
   free_arrays();
   structure= getStructure(myMap,false,sequenceB.length());
   fold(upperCase(sequenceB).c_str(), structure);
   
   PS_rna_plot_a(const_cast<char*>(sequenceB.c_str()),
                 const_cast<char*>(structure),
                 const_cast<char*>(psfilename.c_str()),
                 const_cast<char*>(clus2_str.c_str()),
                 const_cast<char*>(pos1.c_str()));
    free(structure);
    structure= NULL;
    if(base_pair){ free(base_pair); base_pair= NULL;}
    free_arrays();
}

} //end namespace
