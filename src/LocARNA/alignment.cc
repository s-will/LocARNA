#include <math.h>

#include "alignment.hh"
#include "alignment_impl.hh"

#include "sequence.hh"
#include "basepairs.hh"
#include "rna_data.hh"
#include "anchor_constraints.hh"
#include "multiple_alignment.hh"
#include "string1.hh"
#include "plusvector.hh"

#ifdef HAVE_LIBRNA
extern "C" {
#  include <string.h>
#  include <ViennaRNA/data_structures.h>
#  include <ViennaRNA/fold_vars.h>
#  include <ViennaRNA/alifold.h>
}
#endif

namespace LocARNA {

    Alignment::Alignment(const Sequence &seqA, 
			 const Sequence &seqB)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	clear();
    }
    

    Alignment::Alignment(const Sequence &seqA, 
			 const Sequence &seqB, 
			 const string1 &alistrA,
			 const string1 &alistrB)
	: pimpl_(new AlignmentImpl(this,seqA,seqB)) {
	throw failure("Alignment::Alignment(...) constructor from alignment strings not implemented.");
	
	assert(alistrA.length()==alistrB.length());

	//size_t posA;
	//size_t posB;
	for(size_t i=0; i<alistrA.length(); ++i) {
	    //...
	}
    }

    Alignment::Alignment(const Alignment &alignment)
	: pimpl_(new AlignmentImpl(*alignment.pimpl_))
    {
    }
    
    Alignment::~Alignment() { delete pimpl_; }

    Alignment &
    Alignment::operator =(const Alignment &alignment)
    {
	pimpl_ = new AlignmentImpl(*alignment.pimpl_);
	return *this;
    }
    
    void
    Alignment::set_consensus_structure(const RnaStructure &structure) {
	throw("Alignment::set_consensus_structure(...) not implemented");
    }
    
    void
    Alignment::set_structures(const RnaStructure &structureA,const RnaStructure &structureB) {
	throw("Alignment::set_structures(...) not implemented");
    }

    void Alignment::clear() {
	pimpl_->strA_.resize(pimpl_->seqA_.length()+1);
	pimpl_->strB_.resize(pimpl_->seqB_.length()+1);
	for (std::vector<char>::iterator it=pimpl_->strA_.begin(); it!=pimpl_->strA_.end(); ++it) *it='.';
	for (std::vector<char>::iterator it=pimpl_->strB_.begin(); it!=pimpl_->strB_.end(); ++it) *it='.';
	
	pimpl_->a_.clear();
	pimpl_->b_.clear();
    }

    void 
    Alignment::append(int i, int j) {
	pimpl_->a_.push_back(i);
	pimpl_->b_.push_back(j);
    }

    void
    Alignment::add_basepairA(int i, int j) {
	pimpl_->strA_[i]='(';
	pimpl_->strA_[j]=')';
    }
    
    void
    Alignment::add_basepairB(int i, int j) {
	pimpl_->strB_[i]='(';
	pimpl_->strB_[j]=')';
    }
    
    // void
    // invert_position_vector(const std::vector<int> &x,std::vector<int> &inv_a) const {
    // 	throw failure("AlignmentImpl::invert_position_vector(...) not implemented.");

    // }

    void Alignment::write(std::ostream &out, int width, infty_score_t score,
			  bool local_out, bool pos_out, bool write_structure) const {
	write_clustal(out, width, score, local_out, pos_out,false,write_structure);
    }




    void Alignment::write_clustal(std::ostream &out,
				  int width,
				  infty_score_t score,
				  bool local_out,
				  bool pos_out,
				  bool clustal_format,
				  bool write_structure) const {
	if (!pos_out) std::cout<<std::endl;

	const char loc_blank='~';

	size_type alisize=pimpl_->a_.size();

	assert(pimpl_->b_.size()==pimpl_->a_.size());

	/* DEBUGGING output
	   std::cout <<"a: " << a_.size()<<std::endl;
	   for (size_type i=0; i<a_.size(); i++) {
	   std::cout << a_[i]<< " ";
	   }
	   std::cout << std::endl;

	   std::cout <<"b: " << b_.size()<<std::endl;
	   for (size_type i=0; i<b_.size(); i++) {
	   std::cout << b_[i]<< " ";
	   }
	   std::cout << std::endl;
	*/

	Sequence aliStrA;
	Sequence aliSeqA;
	Sequence aliSeqB;
	Sequence aliStrB;
    
	// ----------------------------------------
	// for new sequence objects, generate appropriately named sequence entries with empty sequences 
	aliStrA.append(Sequence::SeqEntry("",""));
	aliStrB.append(Sequence::SeqEntry("",""));
	//
	for(Sequence::const_iterator it=pimpl_->seqA_.begin(); pimpl_->seqA_.end()!=it; ++it)
	    aliSeqA.append(Sequence::SeqEntry(it->name(),""));
	//
	for(Sequence::const_iterator it=pimpl_->seqB_.begin(); pimpl_->seqB_.end()!=it; ++it)
	    aliSeqB.append(Sequence::SeqEntry(it->name(),""));
	// ----------------------------------------
    
    
	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B


	for (size_type i=0; i<alisize; i++) {
	    //out << "("<<a_[i]<<","<<b_[i]<<") "; out.flush();

	    if (!local_out) {
		for (int j=lastA; j<pimpl_->a_[i]; j++) {
		    aliStrA += pimpl_->strA_[j];
		    aliSeqA += pimpl_->seqA_[j]; lastA++;
		    aliSeqB += loc_blank;
		    aliStrB += loc_blank;
		}
		for (int j=lastB; j<pimpl_->b_[i]; j++) {
		    aliStrA += loc_blank;
		    aliSeqA += loc_blank;
		    aliSeqB += pimpl_->seqB_[j]; lastB++;
		    aliStrB += pimpl_->strB_[j];
		}
	    }

	    if ( pimpl_->a_[i]==-1 ) {
		aliStrA += '-';
		aliSeqA += '-';
	    } else if ( pimpl_->a_[i]==-2 ) {
		aliStrA += '_';
		aliSeqA += '_';
	    } else {
		aliStrA += pimpl_->strA_[pimpl_->a_[i]]; lastA++;
		aliSeqA += pimpl_->seqA_[pimpl_->a_[i]];
	    }
	    if ( pimpl_->b_[i]==-1 ) {
		aliStrB += '-';
		aliSeqB += '-';
	    }
	    else if ( pimpl_->b_[i]==-2 ) {
		aliStrB += '_';
		aliSeqB += '_';
	    } else {
		aliStrB += pimpl_->strB_[pimpl_->b_[i]]; lastB++;
		aliSeqB += pimpl_->seqB_[pimpl_->b_[i]];
	    }
	}


	if (!local_out) {
	    for (size_type j=lastA; j<=pimpl_->seqA_.length(); j++) {
		aliStrA += pimpl_->strA_[j];
		aliSeqA += pimpl_->seqA_[j]; lastA++;
		aliSeqB += loc_blank;
		aliStrB += loc_blank;
	    }
	    for (size_type j=lastB; j<=pimpl_->seqB_.length(); j++) {
		aliStrA += loc_blank;
		aliSeqA += loc_blank;
		aliSeqB += pimpl_->seqB_[j]; lastB++;
		aliStrB += pimpl_->strB_[j];
	    }
	}

	if (clustal_format) {
	    out << "CLUSTAL W --- "<<PACKAGE_STRING; // <<" - Local Alignment of RNA"
	    if (pimpl_->seqA_.row_number()==1 && pimpl_->seqB_.row_number()==1)
		out  <<" --- Score: " <<score;
	    out  <<std::endl<<std::endl;
	}
	//else if (!pos_out)
	//	if (pimpl_->seqA_.row_number()==1 && pimpl_->seqB_.row_number()==1)
	//	    out << "SCORE: "<<score<<std::endl;


	if (alisize>0) {
	    size_type local_start_A=pimpl_->a_[0]; //!< pos where local alignment starts for A
	    size_type local_start_B=pimpl_->b_[0]; //!< pos where local alignment starts for B

	    size_type local_end_A=pimpl_->a_[alisize-1]; //!< pos where local alignment ends for A
	    size_type local_end_B=pimpl_->b_[alisize-1]; //!< pos where local alignment ends for B

	    size_type length=aliSeqA.length();

	    if (pos_out) {
		out << "HIT " << score
		    << " " << local_start_A
		    << " " << local_start_B
		    << " " << local_end_A
		    << " " << local_end_B
		    <<std::endl;
	    }


	    if (local_out) {
		out << std::endl;
		out << "\t+" << local_start_A << std::endl;
		out << "\t+" << local_start_B << std::endl;
	    }
	
	    if (!pos_out || local_out) {
		out<<std::endl;

		size_type k=1;
		while (k <= length) {

		    if (write_structure)
			aliStrA.write( out, k, std::min(length,k+width-1) );
		    aliSeqA.write( out, k, std::min(length,k+width-1) );
		    aliSeqB.write( out, k, std::min(length,k+width-1) );

		    if (write_structure)
			aliStrB.write( out, k, std::min(length,k+width-1) );

		    out<<std::endl;
		    k+=width;
		}
	    }

	    if (local_out) {
		out << "\t+" << local_end_A << std::endl;
		out << "\t+" << local_end_B << std::endl;
	    }
	} // end if (alisize>0)
    }


    /*
      write pp output, which can be reread for progressive alignment!
    */

    void Alignment::write_pp(std::ostream &out,
			     const RnaData &rna_dataA,
			     const RnaData &rna_dataB,
			     const AnchorConstraints &seqConstraints,
			     int width,
			     bool use_alifold,
			     double expA,
			     double expB,
			     bool stacking
			     ) const {


	size_type alisize = pimpl_->a_.size();

	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B

	plusvector<int> aliA;
	plusvector<int> aliB;

	for (size_type i=0; i<alisize; i++) {
	    // out << "("<<a_[i]<<","<<b_[i]<<") "; out.flush();
	    for (int j=lastA; j<pimpl_->a_[i]; j++) {
		aliA += j; lastA++;
		aliB += -2;
	    }
	    for (int j=lastB; j<pimpl_->b_[i]; j++) {
		aliA += -2;
		aliB += j; lastB++;
	    }
	    if ( pimpl_->a_[i] < 0 ) {
		aliA += -1;
	    } else {
		aliA += pimpl_->a_[i]; lastA++;
	    }
	    if ( pimpl_->b_[i] < 0 ) {
		aliB += -1;
	    } else {
		aliB += pimpl_->b_[i]; lastB++;
	    }
	}

	for (size_type j=lastA; j<=pimpl_->seqA_.length(); j++) {
	    aliA += j; lastA++;
	    aliB += -2;
	}
	for (size_type j=lastB; j<=pimpl_->seqB_.length(); j++) {
	    aliA += -2;
	    aliB += j; lastB++;
	}
	/*
	  for (int i=0; i<aliA.size(); i++) {
	  std::cerr <<i<<":"<< aliA[i]<<";"<<aliB[i]<<" ";
	  }
	  std::cerr <<std::endl;
	*/

	write(out, width, (infty_score_t)0);

	// ------------------------------------------------------------
	// write consensus constraint string
	//
	if (!seqConstraints.empty()) {
	    out << "#C ";

	    std::vector<std::string> seqCStrings(seqConstraints.name_size());

	    std::string noname="";
	    for (size_type j=0; j<seqConstraints.name_size(); ++j) noname+=".";

	    for (size_type i=0; i<aliA.size(); i++) {
		const std::string &nameA = (aliA[i]>0)?seqConstraints.get_name_a(aliA[i]):"";
		const std::string &nameB = (aliB[i]>0)?seqConstraints.get_name_b(aliB[i]):"";

		std::string name=noname;

		if (nameA!="") name =  nameA;
		else if (nameB!="") name = nameB;

		for (size_type j=0; j<seqConstraints.name_size(); ++j) {
		    seqCStrings[j] += name[j];
		}
	    }

	    for (size_type j=0; j<seqConstraints.name_size(); ++j) {
		if (j!=0) out<<"#";
		out <<seqCStrings[j];
	    }
	}

	out << "\n\n";

	// ------------------------------------------------------------
	// write probability dot plot

	out <<"#"<<std::endl;

#    ifdef HAVE_LIBRNA
	if (use_alifold) {
	    double p_minA = rna_dataA.arc_cutoff_prob();
	    double p_minB = rna_dataB.arc_cutoff_prob();
	    double p_minMean =
		exp(
		    (log(p_minA)*pimpl_->seqA_.row_number()
		     + log(p_minB)*pimpl_->seqB_.row_number())
		    / (pimpl_->seqA_.row_number() + pimpl_->seqB_.row_number())
		    );

	    pimpl_->write_alifold_consensus_dot_plot(out,p_minMean);
	    return;
	}
#    endif

	assert(use_alifold==false /*HAVE_LIBRNA undefined*/);
	pimpl_->write_consensus_dot_plot(out,aliA,aliB,rna_dataA,rna_dataB,expA,expB,stacking);

    }


    size_type
    Alignment::get_local_startA() const {return pimpl_->a_[0];}
    
    size_type
    Alignment::get_local_endA() const {return pimpl_->a_[pimpl_->a_.size()-1];}
    
    size_type
    Alignment::get_local_startB()  const {return pimpl_->b_[0];}
	
    size_type
    Alignment::get_local_endB() const {return pimpl_->b_[pimpl_->b_.size()-1];}

    const Sequence &
    Alignment::get_seqA() const {return pimpl_->seqA_;} 

    const Sequence &
    Alignment::get_seqB() const {return pimpl_->seqB_;} 

    const std::vector<int> &
    Alignment::get_a() const {return pimpl_->a_;} 

    const std::vector<int> &
    Alignment::get_b() const {return pimpl_->b_;} 


#ifdef HAVE_LIBRNA
    void
    AlignmentImpl::write_alifold_consensus_dot_plot(std::ostream &out, double cutoff) const {
	char **sequences;

	plist *pl;

	// set some global variables that control alifold behavior
	RibosumFile = NULL;
	ribo = 1; //activate ribosum scoring
	cv_fact=0.6; // best values for ribosum scoring according to RNAalifold man page
	nc_fact=0.5;
	dangles=2;
	oldAliEn = 0;

	// construct sequences array from sequences in alignment
	MultipleAlignment ma(*self_); // generate multiple alignment from alignment object
	sequences = new char *[ma.row_number()+1];
	for (size_t i=0; i<ma.row_number(); i++) {
	    sequences[i] = new char[ma.length()+1];
	    // copy ma row to sequences[i]
	    strcpy(sequences[i],ma.seqentry(i).seq().to_string().c_str());
	}
	sequences[ma.row_number()]=NULL;

	int length = strlen(sequences[0]);
	char *structure = new char[length+1];

	// estimate specific scaling factor only for large instances
	if (length > 1000) {
	    // estimate nice scaling factor
	    double min_en = alifold((const char **)sequences, structure);

	    double sfact         = 1.07; // from RNAalifold.c code

	    double kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
	    pf_scale = exp(-(sfact*min_en)/kT/length); // set global variable
	} else {
	    pf_scale = 1.0;
	}

	// call ali pf fold
	// double energy = // we don't need the energy
	alipf_fold((const char **)sequences,structure,&pl);

	// read the pair_info structures from array pl
	// and write the base pair probabilities to out
	//
	for(plist *pair=pl; pair->i!=0 && pair->j!=0; pair++) {
	    if ( pair->p  >  cutoff ) {
		out << pair->i << " "
		    << pair->j << " "
		    << pair->p << std::endl;
	    }
	}

	// free heap space
	delete [] structure;
	for (size_t i=0; i<ma.row_number(); i++) {
	    delete [] sequences[i];
	}
	delete [] sequences;

	free(pl);
	free_alifold_arrays();
    }
#endif

    void
    AlignmentImpl::write_consensus_dot_plot(std::ostream &out,
					    const std::vector<int> &aliA,
					    const std::vector<int> &aliB,
					    const RnaData &rna_dataA,
					    const RnaData &rna_dataB,
					    double p_expA,
					    double p_expB,
					    bool stacking
					    ) const {
	
	double p_minA = rna_dataA.arc_cutoff_prob();
	double p_minB = rna_dataB.arc_cutoff_prob();
	
	double p_minMean =
	    exp(
		(log(p_minA)*seqA_.row_number()
		 + log(p_minB)*seqB_.row_number())
		/ (seqA_.row_number() + seqB_.row_number())
		);

	//std::cout << "avg p_min: " << p_minMean << std::endl;
	//std::cout << "p_expA: " << p_expA << std::endl;
	//std::cout << "p_expB: " << p_expB << std::endl;

	for (size_type i=0; i<aliA.size(); i++) {
	    for (size_type j=i+3; j<aliB.size(); j++) { // min loop size=3
		// here we compute consensus pair probabilities

		double pA =
		    (aliA[i]<0 || aliA[j]<0)
		    ? 0
		    : rna_dataA.arc_prob(aliA[i], aliA[j]);

		double pB =
		    (aliB[i]<0 || aliB[j]<0)
		    ? 0
		    : rna_dataB.arc_prob(aliB[i], aliB[j]);

		double p = average_probs(pA,pB,p_minMean,p_expA,p_expB);

		if (stacking) {

		    double st_pA =
			(aliA[i]<0 || aliA[j]<0)
			? 0
			: rna_dataA.joint_arc_prob(aliA[i], aliA[j]);

		    double st_pB =
			(aliB[i]<0 || aliB[j]<0)
			? 0
			: rna_dataB.joint_arc_prob(aliB[i], aliB[j]);

		    double st_p = average_probs(st_pA,st_pB,p_minMean,p_expA,p_expB);

		    if (p > p_minMean || st_p > p_minMean) {
			out << i+1 << " " << j+1 << " " << p << " " << st_p <<std::endl;
			// std::cerr << p <<" <- "<<pA<<","<<pB<<" : "<<p_minA<<","<<p_minB<<";"<<p_minMean<<std::endl;
		    }

		} else {
		    if (p > p_minMean) {
			out << i+1 << " " << j+1 << " " << p <<std::endl;
			// std::cerr << p <<" <- "<<pA<<","<<pB<<" : "<<p_minA<<","<<p_minB<<";"<<p_minMean<<std::endl;
		    }
		}
	    }
	}
    }

    double
    AlignmentImpl::average_probs(double pA, double pB, double p_min,
				 double p_expA, double p_expB) const {

	/*
	  probably better, something like

	  if (pA<p_min*1.05) {
	  pA = std::min(p_expA,p_min*0.75);
	  }


	*/

	pA = std::max(std::min(p_expA,p_min*0.75), pA);
	pB = std::max(std::min(p_expB,p_min*0.75), pB);

	// weighted geometric mean
	double p = exp(
		       (log(pA)*seqA_.row_number() + log(pB)*seqB_.row_number()) /
		       //---------------------------------------------------
		       (seqA_.row_number() + seqB_.row_number())
		       );

	// weighted arithmetic mean
	// p = (seqA_.row_number()*pA + seqB_.row_number()*pB)/(seqA_.row_number()+seqB_.row_number());
	return p;
    }

    void AlignmentImpl::write_debug(std::ostream &out) const {

	for (size_type i=0; i<a_.size(); i++) {
	    out << a_[i] << " ";
	}
	out<<std::endl;
	for (size_type i=0; i<b_.size(); i++) {
	    out << b_[i] << " ";
	}
	out<<std::endl;

    }


}
