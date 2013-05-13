#include "sparsification_mapper.hh"
#include "sequence.hh"


namespace LocARNA {

void SparsificationMapper::compute_mapping_idx_arcs(){
	info_for_pos struct_pos;
	left_adj_vec.resize(bps.num_bps());
	info_valid_seq_pos_vecs.resize(bps.num_bps());
	valid_mat_pos_vecs_before_eq.resize(bps.num_bps());
	for(size_type k=0;k<bps.num_bps();k++){
		//pos_type max_size = 0;
		struct_pos.reset();
		const Arc &arc = bps.arc(k);
		//add initialization
		struct_pos.unpaired=true;
		struct_pos.seq_pos=arc.left();
		info_valid_seq_pos_vecs.at(k).push_back(struct_pos);
		valid_mat_pos_vecs_before_eq.at(k).push_back(0);
		left_adj_vec.at(k).resize(arc.right()-arc.left());
		//compute mapping
		for(size_type j=arc.left()+1;j<arc.right();j++){
			struct_pos.reset();
			if(is_valid_pos(arc,j)){
				struct_pos.seq_pos=j;
				struct_pos.unpaired=true;
			}
			for(BasePairs::RightAdjList::const_iterator inner_arc=bps.right_adjlist(j).begin();
					inner_arc!=bps.right_adjlist(j).end();inner_arc++){
				if(inner_arc->left()<=arc.left()) break;
				if(!is_valid_arc(*inner_arc,arc)) continue;
				left_adj_vec.at(arc.idx()).at(inner_arc->left()-arc.left()).push_back(inner_arc->idx());
				struct_pos.seq_pos=j;//-arc.left();
				struct_pos.valid_arcs.push_back(inner_arc->idx());
			}
			if(struct_pos.seq_pos==j){
				info_valid_seq_pos_vecs.at(k).push_back(struct_pos);
				//max_size++;
			}
			valid_mat_pos_vecs_before_eq.at(k).push_back(info_valid_seq_pos_vecs.at(k).size()-1);

		}
		size_type max_size = info_valid_seq_pos_vecs.at(k).size();
		if (max_info_vec_size < max_size) max_info_vec_size = max_size;
		//if (max_info_vec_size > max_size )
		//    max_info_vec_size = max_size;
	}
	//cout << "valid positions for indices " << info_valid_seq_pos_vecs << endl;
}

void SparsificationMapper::compute_mapping_idx_left_ends(){

	info_for_pos struct_pos;
	size_type seq_length = rnadata.get_sequence().length();
//	std::cout << "compute_mapping_idx_left_ends: seq_length=" << seq_length << std::endl;
	info_valid_seq_pos_vecs.resize(seq_length+1);
	valid_mat_pos_vecs_before_eq.resize(seq_length+1);
	//go over all left ends
	for(pos_type cur_left_end=0;cur_left_end<=seq_length;cur_left_end++){
		size_type max_size = 0;
		struct_pos.reset();
		//add initialization
		struct_pos.unpaired=true;
		struct_pos.seq_pos=cur_left_end;
		info_valid_seq_pos_vecs.at(cur_left_end).push_back(struct_pos);
		valid_mat_pos_vecs_before_eq.at(cur_left_end).push_back(0);
		pos_type max_right_end= (bps.left_adjlist(cur_left_end).begin()==bps.left_adjlist(cur_left_end).end()) ? 0
				:(--bps.left_adjlist(cur_left_end).end())->right();
		if (cur_left_end == 0)
		    max_right_end = seq_length+1;
		for(pos_type cur_pos=cur_left_end+1;cur_pos<max_right_end;cur_pos++){
			struct_pos.reset();
			if (cur_left_end == 0)
			    valid_pos_external(cur_pos,0,struct_pos);
			else
			    iterate_left_adj_list(cur_left_end,cur_pos,0,struct_pos);
			for(BasePairs::RightAdjList::const_iterator inner_arc=bps.right_adjlist(cur_pos).begin();
				inner_arc!=bps.right_adjlist(cur_pos).end();inner_arc++){
			    if(inner_arc->left()<=cur_left_end) break;

			    if (cur_left_end == 0)
				valid_pos_external(cur_pos,&(*inner_arc),struct_pos);
			    else
				iterate_left_adj_list(cur_left_end,cur_pos,&(*inner_arc),struct_pos);
			}
			if(struct_pos.seq_pos==cur_pos){
				info_valid_seq_pos_vecs.at(cur_left_end).push_back(struct_pos);
				max_size++;
			}
			valid_mat_pos_vecs_before_eq.at(cur_left_end).push_back(info_valid_seq_pos_vecs.at(cur_left_end).size()-1);

		}
//		if (max_right_end != 0)		    valid_mat_pos_vecs_before_eq.at(cur_left_end).push_back(info_valid_seq_pos_vecs.at(cur_left_end).size()-1); //toask: ask Christina for max_right

		if (max_info_vec_size < max_size )
		    max_info_vec_size = max_size;
	}
//	cout << "max_info_vec_size " << max_info_vec_size << endl;
//	cout << "valid positions for indices " << info_valid_seq_pos_vecs << endl;
}

void SparsificationMapper::valid_pos_external(pos_type cur_pos,const Arc *inner_arc, info_for_pos &struct_pos){

    if(!inner_arc){
	if(is_valid_pos_external(cur_pos))
	{
	    struct_pos.unpaired=true;
	    struct_pos.seq_pos=cur_pos;
	}
    }
    else if(is_valid_arc_external(*inner_arc))
    {
	struct_pos.valid_arcs.push_back(inner_arc->idx());
	struct_pos.seq_pos=cur_pos;
    }
}

void SparsificationMapper::iterate_left_adj_list(pos_type cur_left_end, pos_type cur_pos,const Arc *inner_arc, info_for_pos &struct_pos){
	for(BasePairs::LeftAdjList::const_iterator arc=bps.left_adjlist(cur_left_end).begin();
			arc!=bps.left_adjlist(cur_left_end).end();arc++){
		if(cur_pos>=arc->right()) continue;
		if(!inner_arc){
			if(!is_valid_pos(*arc,cur_pos)) continue;
			struct_pos.unpaired=true;
			struct_pos.seq_pos=cur_pos;
			break;
		}
		else if(!is_valid_arc(*inner_arc,*arc)) continue;
		struct_pos.valid_arcs.push_back(inner_arc->idx());
		struct_pos.seq_pos=cur_pos;
		break;
	}
}

std::ostream &operator << (std::ostream &out, const vector<SparsificationMapper::InfoForPosVec> &pos_vecs_) {
	size_type idx = 0;
	for (vector<SparsificationMapper::InfoForPosVec>::const_iterator it = pos_vecs_.begin();it!=pos_vecs_.end();it++){
		out << "Idx " << idx << endl;
		out << (*it) << endl;
		idx++;
	}
	return out;
}

std::ostream &operator << (std::ostream &out, const SparsificationMapper::InfoForPosVec &pos_vec_){
	for(SparsificationMapper::InfoForPosVec::const_iterator it_bp = pos_vec_.begin();
			it_bp!=pos_vec_.end();it_bp++){
		out << "pos " << it_bp->seq_pos;
		//int type = it_bp->type_of_pos;
		bool unpaired = it_bp->unpaired;
		if(unpaired) out << " unpaired" ;
		if(!it_bp->valid_arcs.empty()) out << " ArcIdxVec ";
		cout << it_bp->valid_arcs << endl;
	}
	return out;
}

} //end namespace
