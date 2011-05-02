#include <math.h>
#include <string.h>

#include "alignment.hh"
#include "basepairs.hh"
#include "rna_data.hh"
#include "scoring.hh"
#include "anchor_constraints.hh"
#include "multiple_alignment.hh"

#ifdef HAVE_LIBRNA
extern "C" {
#  include "ViennaRNA/data_structures.h"
#  include "ViennaRNA/fold_vars.h"
#  include "ViennaRNA/alifold.h"
}
#endif

namespace LocARNA {

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
	
	size_type alisize=a_.size();
    

	std::vector<Sequence> aliString;
    
	aliString.resize(4,Sequence());

	Sequence &aliStrA=aliString[0];
	Sequence &aliSeqA=aliString[1];
	Sequence &aliSeqB=aliString[2];
	Sequence &aliStrB=aliString[3];

	aliStrA.init_buffer("");
	aliSeqA.init_buffer(seqA_);
	aliSeqB.init_buffer(seqB_);
	aliStrB.init_buffer("");
    
    
	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B
    
	
	for (size_type i=0; i<alisize; i++) {
	    //out << "("<<a_[i]<<","<<b_[i]<<") "; out.flush();
		
	    if (!local_out) {
		for (int j=lastA; j<a_[i]; j++) {
		    aliStrA += strA_[j];
		    aliSeqA += seqA_[j]; lastA++;
		    aliSeqB += loc_blank;
		    aliStrB += loc_blank;
		}
		for (int j=lastB; j<b_[i]; j++) {
		    aliStrA += loc_blank;
		    aliSeqA += loc_blank;
		    aliSeqB += seqB_[j]; lastB++;
		    aliStrB += strB_[j];
		}
	    }
		
	    if ( a_[i]==-1 ) {
		aliStrA += '-';
		aliSeqA += '-';
	    } else {
		aliStrA += strA_[a_[i]]; lastA++;
		aliSeqA += seqA_[a_[i]];
	    }
	    if ( b_[i]==-1 ) {
		aliStrB += '-';
		aliSeqB += '-';
	    } else {
		aliStrB += strB_[b_[i]]; lastB++;
		aliSeqB += seqB_[b_[i]];
	    }
	}
       

	if (!local_out) {
	    for (size_type j=lastA; j<=seqA_.length(); j++) {
		aliStrA += strA_[j];
		aliSeqA += seqA_[j]; lastA++;
		aliSeqB += loc_blank;
		aliStrB += loc_blank;
	    }
	    for (size_type j=lastB; j<=seqB_.length(); j++) {
		aliStrA += loc_blank;
		aliSeqA += loc_blank;
		aliSeqB += seqB_[j]; lastB++;
		aliStrB += strB_[j];
	    }
	}

	if (clustal_format) {
	    out << "CLUSTAL W --- "<<PACKAGE_STRING; // <<" - Local Alignment of RNA"
	    if (seqA_.row_number()==1 && seqB_.row_number()==1)
		out  <<" --- Score: " <<score;
	    out  <<std::endl<<std::endl;
	}
	//else if (!pos_out)
	//	if (seqA_.row_number()==1 && seqB_.row_number()==1)
	//	    out << "SCORE: "<<score<<std::endl;
    
	size_type local_start_A=a_[0]; //!< pos where local alignment starts for A
	size_type local_start_B=b_[0]; //!< pos where local alignment starts for B
		
	size_type local_end_A=a_[alisize-1]; //!< pos where local alignment ends for A
	size_type local_end_B=b_[alisize-1]; //!< pos where local alignment ends for B

	size_type length=aliSeqA.length();
	size_type k=1;
 
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
	
	    while (k <= length) {
	    
		if (write_structure)
		    aliString[0].write( out, k, std::min(length,k+width-1) );
		aliString[1].write( out, k, std::min(length,k+width-1) );
		aliString[2].write( out, k, std::min(length,k+width-1) );

		if (write_structure)
		    aliString[3].write( out, k, std::min(length,k+width-1) );

		out<<std::endl;
		k+=width;
	    }
	}

	if (local_out) {
	    out << "\t+" << local_end_A << std::endl;
	    out << "\t+" << local_end_B << std::endl;
	}
    }


    /*
      write pp output, which can be reread for progressive alignment!
    */

    void Alignment::write_pp(std::ostream &out,
			     const BasePairs &bpsA,
			     const BasePairs &bpsB,
			     const Scoring &scoring, 
			     const AnchorConstraints &seqConstraints, 
			     int width,
			     bool use_alifold) const { 


	size_type alisize = a_.size();

	int lastA=1; // bases consumed in sequence A
	int lastB=1; // ---------- "" ------------ B

	plusvector<int> aliA;
	plusvector<int> aliB;

	for (size_type i=0; i<alisize; i++) {
	    // out << "("<<a_[i]<<","<<b_[i]<<") "; out.flush();
	    for (int j=lastA; j<a_[i]; j++) {
		aliA += j; lastA++;
		aliB += -2;
	    }
	    for (int j=lastB; j<b_[i]; j++) {
		aliA += -2;
		aliB += j; lastB++;
	    }
	    if ( a_[i]==-1 ) {
		aliA += -1;
	    } else {
		aliA += a_[i]; lastA++;
	    }
	    if ( b_[i]==-1 ) {
		aliB += -1;
	    } else {
		aliB += b_[i]; lastB++;
	    }
	}
  
	for (size_type j=lastA; j<=seqA_.length(); j++) {
	    aliA += j; lastA++;
	    aliB += -2;
	}
	for (size_type j=lastB; j<=seqB_.length(); j++) {
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
	    double p_minA = bpsA.prob_min();
	    double p_minB = bpsB.prob_min();    
	    double p_minMean =
		exp(
		    (log(p_minA)*seqA_.row_number()
		     + log(p_minB)*seqB_.row_number())
		    / (seqA_.row_number() + seqB_.row_number())
		    );
	    
	    write_alifold_consensus_dot_plot(out,p_minMean);
	    return;
	}
#    endif

	assert(use_alifold==false /*HAVE_LIBRNA undefined*/);
	write_consensus_dot_plot(out,aliA,aliB,bpsA,bpsB,scoring);

    }

#ifdef HAVE_LIBRNA
    void
    Alignment::write_alifold_consensus_dot_plot(std::ostream &out, double cutoff) const {
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
	MultipleAlignment ma(*this); // generate multiple alignment from alignment object
	sequences = new char *[ma.size()+1];
	for (size_t i=0; i<ma.size(); i++) {
	    sequences[i] = new char[ma.length()+1];
	    // copy ma row to sequences[i]
	    strcpy(sequences[i],ma.seqentry(i).seq().to_string().c_str());
	}
	sequences[ma.size()]=NULL;
	
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
	double energy = alipf_fold((const char **)sequences,structure,&pl);
	
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
	delete structure;
	for (size_t i=0; i<ma.size(); i++) {
	    delete sequences[i];
	}
	delete sequences;
	
	free(pl);
	free_alifold_arrays();
    }
#endif
    
    void
    Alignment::write_consensus_dot_plot(std::ostream &out,
					const plusvector<int> &aliA,
					const plusvector<int> &aliB,
					const BasePairs &bpsA,
					const BasePairs &bpsB,
					const Scoring &scoring 
					) const {
	int lenA=seqA_.length();
	int lenB=seqB_.length();
    
	double p_minA = bpsA.prob_min();
	double p_minB = bpsB.prob_min();
	double p_expA = scoring.prob_exp(lenA);
	double p_expB = scoring.prob_exp(lenB);
    
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
		    : bpsA.get_arc_prob(aliA[i], aliA[j]);

		double pB = 
		    (aliB[i]<0 || aliB[j]<0)
		    ? 0
		    : bpsB.get_arc_prob(aliB[i], aliB[j]);
	    
		double p = average_probs(pA,pB,p_minMean,p_expA,p_expB);	    
	    
		if (scoring.stacking()) {
		
		    double st_pA =
			(aliA[i]<0 || aliA[j]<0)
			? 0
			: bpsA.get_arc_2_prob(aliA[i], aliA[j]);
		    
		    double st_pB = 
			(aliB[i]<0 || aliB[j]<0)
			? 0
			: bpsB.get_arc_2_prob(aliB[i], aliB[j]);
		    
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
    Alignment::average_probs(double pA, double pB, double p_min,
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

    void Alignment::write_debug(std::ostream &out) const {
    
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
