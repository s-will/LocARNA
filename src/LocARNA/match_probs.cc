#include "match_probs.hh"

/* 
   Nachteile von Probcons(=pairHMM)-Ansatz

   * nicht local oder einfach erweiterbar
   * Parameter von ProbconsRNA nicht einfach nachvollziehbar
   * Implementierung kritisch (nochmal checken)
   * muss noch in log-space transformiert werden für numerische Stabilität
   
   => probalign / proba Ansatz
   
*/

/*
  ACHTUNG: die Algorithmen müssen letzlich auf Sequenzen von Alignmentspalten laufen!
  Momentan wird in pairHMM_probs nur die erste Zeile benutzt! (seqX[i][0]!)

  Der pf-Algorithmus sollte damit klarkommen. Das ist aber sehr einfach
  implementiert (d.h. langsam).
*/



#include <iostream>
#include <fstream>
#include <sstream>

#include <cmath>

#include "sequence.hh"
#include "alphabet.hh"
#include "ribosum.hh"
#include "stral_score.hh"

namespace  LocARNA {
    
    MatchProbs::ProbConsParameter::ProbConsParameter(const std::string &filename) {
	std::ifstream in(filename.c_str());
	if (!in.good()) {
	    std::cerr << "Cannot open file "<<filename<<" for reading."<<std::endl;
	    exit(-1);
	}
	
	try {
	    in >> initM >> initX >> initY;
	    
	    in >> startX >> startY;
	    
	    extendM = 1-startX-startY;
	    
	    in >> extendX >> extendY;
	    startMFromX = 1 - extendX;
	    startMFromY = 1 - extendY;   
	    
	    getline(in,basenames); // eat line before base names
	    getline(in,basenames); // read base names
	    
	    if (basenames!="ACGUTN") {
		throw(std::ifstream::failure("Expected base names ACGUTN. Found line: "+basenames));
	    }
	    
	    // read emmission probs
	    emmission.resize(6,6);
	    
	    std::string line;
	    int i=0;
	    while (i<6 && getline(in,line)) {
		std::istringstream sline(line);
		for (int j=0; j<=i; j++) {
		    double p;
		    sline >> p;
		    emmission(i,j)=p;
		    emmission(j,i)=p;
		}
		i++;
	    }
	    if (i!=6) {
		throw(std::ifstream::failure("Cannot read enough emmission probabilities."));
	    }

	    // read background probabilities

	    background.resize(6);
	    if (getline(in,line)) {
		std::istringstream sline(line);
		for (int j=0; j<6; j++) {
		    double p;
		    sline >> p;
		    background[j]=p;
		}
	    } else {
		throw(std::ifstream::failure("Cannot read background probabilities."));
	    }
	} catch (std::ifstream::failure e) {
	    std::cerr << "Cannot parse "<<filename<<". " <<e.what()<< std::endl
		      << "File not in probcons parameter format. Now exiting." << std::endl;
	    exit(-1);
	}
    }

    MatchProbs::MatchProbs()
	: probs() {
    }

    /*
      MatchProbs::MatchProbs(const std::string &filename)
      : probs() {
      read(filename);
      }
    */

    void
    MatchProbs::pairHMM_probs(const Sequence &seqA,
			      const Sequence &seqB,
			      const std::string &filename) {

	typedef size_t size_type;
    
	ProbConsParameter p(filename);    
    
	// use a translation table to convert nucleotide characters to indices
    
	std::vector<int> trans_tab(256,5);
	for(size_type i=0; i<p.basenames.length(); i++) {
	    char c=p.basenames[i];
	    trans_tab[c]=i;
	}
    
	size_t lenA=seqA.length();
	size_t lenB=seqB.length();
    
	if (seqA.get_rows()!=1 || seqB.get_rows()!=1) {
	    std::cerr
		<< "WARNING: the base match probabilities are currently computed only on the first sequence" << std::endl
		<< "of a multiple alignment. I.e., this does not work correctly for multiple alignment yet." << std::endl;
	}
    
	// ------------------------------------------------------------
	// calculate the forward variables
    
	Matrix<double> fwdM(lenA+1,lenB+1);
	Matrix<double> fwdX(lenA+1,lenB+1);
	Matrix<double> fwdY(lenA+1,lenB+1);
	// due to the constructors all entries in fwd* are initialized to 0
    
	// init
	fwdM(0,0)=p.initM;
	fwdX(0,0)=p.initX;
	fwdY(0,0)=p.initY;

    
	// init first column of fwdX
	for (size_type i=1; i<=seqA.length(); i++) {
	    int Ai=trans_tab[seqA[i][0]];
	    fwdX(i,0) = fwdX(i-1,0) * p.extendX * p.background[Ai];
	}

	// init first row of fwdY
	for (size_type j=1; j<=seqB.length(); j++) {
	    int Bj=trans_tab[seqB[j][0]];
	    fwdY(0,j) = fwdY(0,j-1) * p.extendY * p.background[Bj];
	}

    
	// recursion
	for (size_type i=1; i<=seqA.length(); i++) {
	    for (size_type j=1; j<=seqB.length(); j++) {
		size_type Ai=trans_tab[seqA[i][0]];
		size_type Bj=trans_tab[seqB[j][0]];
	    
		fwdM(i,j) = ( fwdM(i-1,j-1) * p.extendM 
			      + fwdX(i-1,j-1) * p.startMFromX
			      + fwdY(i-1,j-1) * p.startMFromY
			      )
		    * p.emmission(Ai,Bj);
	    
		fwdX(i,j) = ( fwdM(i-1,j) * p.startX
			      + fwdX(i-1,j) * p.extendX ) * p.background[Ai];
	    
		fwdY(i,j) = ( fwdM(i,j-1) * p.startY
			      + fwdY(i,j-1) * p.extendY ) * p.background[Bj];
	    }
	}
    
	// ------------------------------------------------------------
	// calculate the backward variables


	Matrix<double> bckM(lenA+1,lenB+1);
	Matrix<double> bckX(lenA+1,lenB+1);
	Matrix<double> bckY(lenA+1,lenB+1);
	// due to the constructors all entries in bck* are initialized to 0
    
	// init
	bckM(lenA,lenB)=1;
	bckX(lenA,lenB)=1;
	bckY(lenA,lenB)=1;
    
	// init first column of bckX
	for (int i=seqA.length()-1; i>0; i--) {
	    int Ai=trans_tab[seqA[i+1][0]];
	
	    bckM(i,lenB) = p.startX * p.background[Ai] * bckX(i+1,lenB);
	    bckX(i,lenB) = p.extendX * p.background[Ai] * bckX(i+1,lenB);
	    bckY(i,lenB) = 0;
	}
    
	// init first row of bckY
	for (int j=seqB.length()-1; j>0; j--) {
	    int Bj=trans_tab[seqB[j+1][0]];
	    bckM(lenA,j) = p.startY * p.background[Bj] * bckY(lenA,j+1);
	    bckX(lenA,j) = 0;
	    bckY(lenA,j) = p.extendY * p.background[Bj] * bckY(lenA,j+1);
	}
    
	// recursion
	for (int i=seqA.length()-1; i>0; i--) {
	    for (int j=seqB.length()-1; j>0; j--) {
		int Ai=trans_tab[seqA[i][0]];
		int Bj=trans_tab[seqB[j][0]];
	    
		bckM(i,j) =
		    p.extendM * p.emmission(Ai,Bj) * bckM(i+1,j+1)
		    + p.startX * p.background[Ai] * bckX(i+1,j)
		    + p.startY * p.background[Bj] * bckY(i,j+1)
		    ;
	
		bckX(i,j) = 
		    p.startMFromX * p.emmission(Ai,Bj) * bckM(i+1,j+1)
		    + p.extendX * p.background[Ai] * bckX(i+1,j);

		bckY(i,j) = 
		    p.startMFromY * p.emmission(Ai,Bj) * bckM(i+1,j+1)
		    + p.extendY * p.background[Bj] * bckY(i,j+1);
	    }
	}
    

	double pAB = fwdM(lenA,lenB)+fwdX(lenA,lenB)+fwdY(lenA,lenB);

	// now compute base match probabilities
	probs.resize(lenA+1,lenB+1);
	for (size_type i=1; i<=seqA.length(); i++) {
	    for (size_type j=1; j<=seqB.length(); j++) {
		probs(i,j) = fwdM (i,j)*bckM(i,j) / pAB;
		//std::cout << i<<","<<j<<": "<<probs(i,j)<<std::endl;
	    }
	}
    }



    //! perform the partition variant of Gotoh on seqA and seqB
    //! and fill the three matrices
    void MatchProbs::pf_gotoh(size_type lenA,
			      size_type lenB,
			      Matrix<double> &zM,
			      Matrix<double> &zA,
			      Matrix<double> &zB,
			      const StralScore &score,
			      double temp,
			      bool local
			      ) {

	/*
	  Gotoh partition version
      
	  In contrast to the standard Gotoh algorithm,
	  we need to make the algorithm non-ambiguous.
	  This means that ZM is only for alignments that
	  do not end in gaps. 
	  We achieve this by reducing ZM_ij to
	  ZM_i-1,j-1, ZA_i-1,j-1, and ZB_i-1,j-1
      
	  The matrices ZA and ZB represent alignments
	  that end in a gap in a or b, resp.
      
	*/

	// Boltzman-weights for gap opening and extension
	double g_open = exp( score.alpha() / temp );
	double g_ext  = exp( score.beta() / temp );
    
	// std::cout << "g_open: "<<g_open<<std::endl;
	// std::cout << "g_ext: "<<g_ext<<std::endl;

	zM.resize(lenA+1,lenB+1);
	zA.resize(lenA+1,lenB+1);
	zB.resize(lenA+1,lenB+1);
    
	// initialization
	//
	// we start with entries that are equal for
	// global and local alignment

	zM(0,0) = local?0:1;
	zA(0,0) = 0;
	zB(0,0) = 0;
    
	for (size_type i=1; i<=lenA; i++) { zM(i,0) = 0; }
	for (size_type j=1; j<=lenB; j++) { zM(0,j) = 0; }
    
	for (size_type j=1; j<=lenB; j++) { zA(0,j) = 0; }
	for (size_type i=1; i<=lenA; i++) { zB(i,0) = 0; }
    
	zA(1,0)=g_open * g_ext;
	zB(0,1)=g_open * g_ext;
	
	// init that differs for global and local
	for (size_type i=2; i<=lenA; i++) {
	    zA(i,0) = ((local?g_open:0) + zA(i-1,0)) * g_ext;
	}
    
	for (size_type j=2; j<=lenB; j++) {
	    zB(0,j) = ((local?g_open:0) + zB(0,j-1)) * g_ext;
	}
    
	// recursion
	for (size_type i=1; i<=lenA; i++) {
	    for (size_type j=1; j<=lenB; j++) {
	    
		double match_score_ij = score.sigma(i,j);
	    
	    
		// Boltzman-weight for match of i and j
		double match_ij = exp( match_score_ij / temp );
	    
		zM(i,j) = 
		    + zM(i-1,j-1) * match_ij
		    + zA(i-1,j-1) * match_ij
		    + zB(i-1,j-1) * match_ij
		    + (local?match_ij:0)
		    ;
	    
		zA(i,j) = 
		    + zA(i-1,j) * g_ext
		    + zM(i-1,j) * g_open * g_ext
		    + zB(i-1,j) * g_open * g_ext
		    + (local?(g_open * g_ext):0)
		    ;
		   
		zB(i,j) =
		    + zB(i,j-1) * g_ext
		    + zM(i,j-1) * g_open * g_ext
		    + zA(i,j-1) * g_open * g_ext
		    + (local?(g_open * g_ext):0)
		    ;
	    }
	}
    }


    void
    MatchProbs::pf_probs(const RnaData &rnaA,
			 const RnaData &rnaB, 
			 const Matrix<double> &sim_mat,
			 const Alphabet<char> &alphabet,
			 double gap_opening,
			 double gap_extension,
			 double pf_struct_weight,
			 double temp,
			 bool flag_local)
    {
	
	size_type lenA=rnaA.get_sequence().length();
	size_type lenB=rnaB.get_sequence().length();

	Matrix<double> zM;
	Matrix<double> zA;
	Matrix<double> zB;
    
	Matrix<double> zMr;
	Matrix<double> zAr;
	Matrix<double> zBr;
    
	StralScore score( rnaA, rnaB,
			  sim_mat, alphabet, 
			  pf_struct_weight,
			  gap_opening, 
			  gap_extension
			  );
    
	pf_gotoh(lenA,lenB,
		 zM,zA,zB,
		 score,
		 temp,
		 flag_local
		 );
    
	score.reverse();
    
	pf_gotoh(lenA,lenB,
		 zMr,zAr,zBr,
		 score, // reversed !
		 temp,
		 flag_local
		 );
    
	double z; // total partition function

	if (flag_local) {
	    // for the local pf we need to sum over all matrix entries
	    z = 1; // weight of the empty alignment
	    for (size_type i=0; i<=lenA; i++) {
		for (size_type j=0; j<=lenB; j++) {
		    z += zM(i,j)+zA(i,j)+zB(i,j);
		}
	    }
	    for (size_type i=0; i<=lenA; i++) {
		z -= lenB * zA(i,0);
	    }
	    for (size_type j=0; j<=lenB; j++) {
		z -= lenA * zB(0,j);
	    }
	} else { // global
	    z = zM(lenA,lenB)+zA(lenA,lenB)+zB(lenA,lenB);
	}
    
	// std::cout << "Z=" << z << std::endl;
	// std::cout << " ?= " << zMr(lenA,lenB)+zAr(lenA,lenB)+zBr(lenA,lenB) << std::endl;
    
	/*
	  std::cout << "ZM:" << std::endl << zM << std::endl;
	  std::cout << "ZA:" << std::endl << zA << std::endl;
	  std::cout << "ZB:" << std::endl << zB << std::endl;
	  std::cout << "ZMr:" << std::endl << zMr << std::endl;
	  std::cout << "ZAr:" << std::endl << zAr << std::endl;
	  std::cout << "ZBr:" << std::endl << zBr << std::endl;
	*/
    
	probs.resize(lenA+1,lenB+1);
    
	// in local alignment, we need to add 1 for the empty alignment
	// for avoiding redundancy the weight of the empty alignment is
	// not included in either matrix entry
	// 
	double locality_add = (flag_local?1:0);
    
	for (size_type i=1; i<=lenA; i++) {
	    for (size_type j=1; j<=lenB; j++) {
		probs(i,j) = 
		    (zM(i,j) * (zMr(lenA-i,lenB-j) + zAr(lenA-i,lenB-j) + zBr(lenA-i,lenB-j) + locality_add ) )
		    / z ;
		// std::cout <<i<<" "<<j<<": "<<probs(i,j)<<std::endl;
	    }
	}
    
	// std::cout << "Probs:" << std::endl << probs << std::endl;
    }

    /*
      void
      MatchProbs::read(const std::string &filename) {
      std::ifstream in(filename.c_str());
      read(in);
      }


      std::istream &
      MatchProbs::read(std::istream &in) {
      size_type x,y;
      try {
	
      in >> x >> y;
	
      if (x<=0 || y<=0) throw std::istream::failure("Found invalid dimensions.");
	
      probs.resize(x,y);
	
      in >> probs;
	
      } catch (std::istream::failure e) {
      std::cerr << "Parse error "<<e.what()<<std::endl;
      exit(-1);
      }
    
      return in;
      }
    */


    void
    MatchProbs::read_sparse(const std::string &filename, size_type lenA, size_type lenB) {    
	std::ifstream in(filename.c_str());
	read_sparse(in,lenA,lenB);
    }

    std::istream &
    MatchProbs::read_sparse(std::istream &in, size_type lenA, size_type lenB) {
	probs.resize(lenA+1,lenB+1);

	probs.fill(0);

	size_type i,j;
	double p;
    
	while (in >> i >> j >> p) {
	    probs(i,j)=p;
	}
    
	return in;
    }


    /*
      void
      MatchProbs::write(const std::string &filename) const {
      std::ofstream out(filename.c_str());
      write(out);
      }

      std::ostream &
      MatchProbs::write(std::ostream &out) const {
      out << probs.sizes().first << " " << probs.sizes().second << std::endl;
      out << probs;
      return out;
      }
    */

    void
    MatchProbs::write_sparse(const std::string &filename, double threshold) const {
	std::ofstream out(filename.c_str());
	write_sparse(out, threshold);
    }

    std::ostream &
    MatchProbs::write_sparse(std::ostream &out, double threshold) const {
       
	size_type lenA=probs.sizes().first-1;
	size_type lenB=probs.sizes().second-1;

	for (size_type i=1; i<=lenA; i++) {
	    for (size_type j=1; j<=lenB; j++) {
		if (probs(i,j) >= threshold) {
		    out << i<<" "<<j<<" "<<probs(i,j)<<std::endl;
		}
	    }
	}
	return out;
    }

} // end namespace LocARNA
