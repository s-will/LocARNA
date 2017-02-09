#include "ribosum.hh"
#include "aux.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

namespace LocARNA {

    Ribosum::Ribosum(const std::string &filename)
        : name(),
          bm(),
          am(),
          basename_alphabet(),
          arcname_alphabet(),
          char_basename_alphabet()
    {
        std::ifstream in(filename.c_str());
        if (!in) {
            std::ostringstream err;
            err << "Cannot open file "<<filename<<" for reading ribosum data.";
            throw failure(err.str());
        }

        read_ribosum(in);

        in.close();
    }

    Ribosum::Ribosum()
        : name(),
          bm(),
          am(),
          basename_alphabet(),
          arcname_alphabet(),
          char_basename_alphabet()
    {}


    Ribosum::~Ribosum() {
    }

    void
    Ribosum::read_ribosum(std::istream &in) {

        try {
            std::string line;

            if (! std::getline(in,line) )
                throw(std::ifstream::failure("Expecting name."));
            name = line;

            // ------------------------------------------------------------
            // read the 4x4 matrix

            std::vector<std::string> basenames;
            basenames.push_back("A");
            basenames.push_back("C");
            basenames.push_back("G");
            basenames.push_back("U");

            basename_alphabet = Alphabet<std::string>(basenames);
            char_basename_alphabet = make_char_alphabet();

            read_matrix(in,bm,basename_alphabet);

            // ignore two lines -- not clean, but we don't need the info currently
            getline(in,line);
            getline(in,line);

            // ------------------------------------------------------------
            // read the 16x16 matrix

            std::vector<std::string> arcnames;
            for (std::vector<std::string>::const_iterator i=basenames.begin();
                 i!=basenames.end();
                 ++i)
                for (std::vector<std::string>::const_iterator j=basenames.begin();
                     j!=basenames.end();
                     ++j)
                    arcnames.push_back((*i)+(*j));

            arcname_alphabet = Alphabet<std::string>(arcnames);

            read_matrix(in,am,arcname_alphabet);

            // ignore two lines
            getline(in,line);
            getline(in,line);

        } catch (std::ifstream::failure &e) {
            std::ostringstream err;
            err << "Cannot parse ribosum input. " <<e.what()<< std::endl
                << "File not in ribosum format.";
            throw failure(err.str());
        }
    }

    Alphabet<char>
    Ribosum::make_char_alphabet() const {
        // transform alphabet over singleton strings into one over characters
        std::vector<char> alphchars;
        for (Ribosum::alphabet_type::const_iterator it=basename_alphabet.begin();
             it!=basename_alphabet.end();
             ++it)
            alphchars.push_back((*it)[0]);

        return Alphabet<char>(alphchars);
    }

    std::istream &
    Ribosum::read_matrix(std::istream &in, Matrix<double> &mat, const Alphabet<std::string> &alph) const {
        Alphabet<std::string>::size_type siz=alph.size();

        std::string line;

        while (std::getline(in,line) && line == "")
            ;

        {
            std::istringstream linestream(line);

            for (size_t i=0; i<siz; i++) {
                std::string name;
                linestream >> name;

                if (name!=alph.elem(i))
                    throw(std::ifstream::failure("Expecting correct table header. Found: "+line));
            }
        }

        mat.resize(siz,siz);
        for(size_t i=0; i<siz; i++) {
            std::getline(in,line);
            std::istringstream linestream(line);
            std::string base;
            linestream >> base;
            if (base != alph.elem(i))
                throw(std::ifstream::failure("Expecting base name "+ alph.elem(i) +" as row header"));

            for (size_t j=0; j<=i; j++) {
                double number;
                linestream >> number;
                mat(i,j) = mat(j,i) = number;
            }
        }
        return in;
    }

    std::ostream &
    Ribosum::write_matrix(std::ostream &out, const Matrix<double> &mat, const Alphabet<std::string> &alph) const {
        out << alph << std::endl;
        out << mat  << std::endl;
        return out;
    }

    /**
     * Output operator
     *
     * @param out output stream
     * @param ribosum Ribosum to be written to stream
     *
     * @return output stream after writing
     */
    std::ostream & operator << (std::ostream &out, const Ribosum &ribosum) {
        out << ribosum.name << std::endl << std::endl;

        ribosum.write_matrix(out,ribosum.bm,ribosum.basename_alphabet);
        ribosum.write_matrix(out,ribosum.am,ribosum.arcname_alphabet);

        return out;
    }


    // ================================================================================
    // class RibosumFreq


    RibosumFreq::RibosumFreq(const std::string &filename)
        : Ribosum(),
          base_probs_(),
          base_nonstruct_probs_(),
          basepair_probs_(),
          basematch_probs_(),
          arcmatch_probs_()
    {
        std::ifstream in(filename.c_str());
        if (!in.is_open()) {
            std::cerr << "Cannot open file "<<filename<<" for reading ribosum data."<<std::endl;
            exit(-1);
        }
        read_ribosum(in);

        read_frequencies(in);

        in.close();
    }

    RibosumFreq::RibosumFreq()
        : Ribosum(),
          base_probs_(),
          base_nonstruct_probs_(),
          basepair_probs_(),
          basematch_probs_(),
          arcmatch_probs_() {
    }

    /**
     * Test for only blank characters
     *
     * @param s string to be tested
     *
     * @return whether s consists of only blank characters
     */
    bool is_blank(std::string &s) {
        for (size_t i=0; i<s.length(); i++)
            if (s[i]!=' ') return false;
        return true;
    }


    void
    RibosumFreq::read_frequencies(std::istream &in)  {
        read_matrix(in,"BASE FREQUENCIES",base_probs_,4,1);
        read_matrix(in,"BASE NONSTRUCTURAL FREQUENCIES",base_nonstruct_probs_,4,1);
        read_matrix(in,"BASE PAIR FREQUENCIES",basepair_probs_,4,4);
        read_matrix(in,"BASE MATCH FREQUENCIES",basematch_probs_,4,4);
        read_matrix(in,"BASE PAIR MATCH FREQUENCIES",arcmatch_probs_,16,16);}

    void
    RibosumFreq::read_matrix(std::istream &in,
                             const std::string &header,
                             matrix_t &mat,
                             size_t xdim,
                             size_t ydim) {
        try {
            std::string line;
            while (std::getline(in,line) && is_blank(line))
                ;

            if (line != header) {
                throw std::ifstream::failure("Expected header "+header+"."+" Read instead '"+line+"'.");
            }

            mat.resize(xdim,ydim);

            for (size_t i=0; i<xdim; i++)
                for (size_t j=0; j<ydim; j++)
                    in >> mat(i,j);

        } catch(std::ifstream::failure &e) {
            std::ostringstream err;
            err << "Cannot parse ribosum frequency input. " <<e.what()<< std::endl
                << "File not in extended ribosum format.";
            throw failure(err.str());
        }
    }


    std::ostream &
    RibosumFreq::write_matrix(std::ostream &out, const std::string &name, const Matrix<double> &mat) const {
        out << name << std::endl;
        out << mat  << std::endl;
        return out;
    }


    /**
     * Output operator
     *
     * @param out output stream
     * @param ribosum Ribosum to be written to stream
     *
     * @return output stream after writing
     */
    std::ostream & operator << (std::ostream &out, const RibosumFreq &ribosum) {
        out << (Ribosum)(ribosum);
        out << std::endl;

        ribosum.write_matrix(out,"BASE FREQUENCIES",ribosum.get_base_probs());
        ribosum.write_matrix(out,"BASE NONSTRUCTURAL FREQUENCIES",ribosum.get_base_nonstruct_probs());
        ribosum.write_matrix(out,"BASE PAIR FREQUENCIES",ribosum.get_basepair_probs());
        ribosum.write_matrix(out,"BASE MATCH FREQUENCIES",ribosum.get_basematch_probs());
        ribosum.write_matrix(out,"ARC MATCH FREQUENCIES",ribosum.get_arcmatch_probs());

        return out;
    }


    void RibosumFreq::write_CC_matrix(std::ostream &out,
                                      const std::string &ribname,
                                      const std::string &matname,
                                      int x, int y, const Ribosum::matrix_t &m) const{
        out << "const double "<<ribname<<"::"<<matname<<"[] = {" <<std::endl;
        for (int i=0;i<x;i++) {
            out << "    ";
            for (int j=0;j<y;j++) {
                out << m(i,j);
                if (i<(x-1) || j<(y-1)) out << ", "; else out << " ";
            }
            out << std::endl;
        }
        out << "};"<<std::endl<<std::endl;
    }

    void
    RibosumFreq::write_ICC_code(std::ostream &out, const std::string &ribname) const {
        std::string macro_headerfile=ribname;
        transform(macro_headerfile.begin(),
                  macro_headerfile.end(),
                  macro_headerfile.begin(),
                  ::toupper);

        macro_headerfile="LOCARNA_"+macro_headerfile+"_HH";

        out << "#ifndef "<<macro_headerfile<<std::endl
            << "#define "<<macro_headerfile<<std::endl
            << "#include \"LocARNA/ribosum.hh\""<<std::endl
            << "namespace LocARNA {"<<std::endl
            << "class "<<ribname<<": public RibosumFreq {" << std::endl
            << "    static const double bm_init[];" << std::endl
            << "    static const double am_init[];" << std::endl
            << "    static const double base_probs_init[];" << std::endl
            << "    static const double base_nonstruct_probs_init[];" << std::endl
            << "    static const double basepair_probs_init[];" << std::endl
            << "    static const double basematch_probs_init[];" << std::endl
            << "    static const double arcmatch_probs_init[];" << std::endl
            << "    static const std::string basename_alphabet_init[];"<< std::endl
            << "    static const std::string arcname_alphabet_init[];" << std::endl
            << std::endl
            << "  public:"<<std::endl
            << "    "<<ribname<<"(): RibosumFreq() {" << std::endl
            << "        bm = matrix_t(4,4,bm_init);" << std::endl
            << "        am = matrix_t(16,16,am_init);" << std::endl

            << "        base_probs_ = matrix_t(4,1,base_probs_init);" << std::endl
            << "        base_nonstruct_probs_ = matrix_t(4,1,base_nonstruct_probs_init);" << std::endl
            << "        basepair_probs_ = matrix_t(4,4,basepair_probs_init);" << std::endl
            << "        basematch_probs_ = matrix_t(4,4,basematch_probs_init);" << std::endl
            << "        arcmatch_probs_ = matrix_t(16,16,arcmatch_probs_init);" << std::endl
            << "        set_basename_alphabet(basename_alphabet_init);" << std::endl
            << "        set_arcname_alphabet(arcname_alphabet_init);" << std::endl
            << "    }" << std::endl
            << "};" << std::endl;

        write_CC_matrix(out, ribname,"bm_init",4,4,get_basematch_scores());
        write_CC_matrix(out, ribname,"am_init",16,16,get_arcmatch_scores());

        write_CC_matrix(out, ribname,"base_probs_init",4,1,get_base_probs());
        write_CC_matrix(out, ribname,"base_nonstruct_probs_init",
                        4,1,get_base_nonstruct_probs());

        write_CC_matrix(out, ribname,"basepair_probs_init",4,4,get_basepair_probs());
        write_CC_matrix(out, ribname,"basematch_probs_init",4,4,get_basematch_probs());

        write_CC_matrix(out, ribname,"arcmatch_probs_init",16,16,get_arcmatch_probs());

        out << "const std::string "<<ribname<<"::basename_alphabet_init[] = "
            << "{\"A\",\"C\",\"G\",\"U\"};" << std::endl << std::endl

            << "const std::string "<<ribname<<"::arcname_alphabet_init[] = "
            << "{\"AA\",\"AC\",\"AG\",\"AU\",\"CA\",\"CC\",\"CG\",\"CU\",\"GA\",\"GC\",\"GG\",\"GU\",\"UA\",\"UC\",\"UG\",\"UU\"};" << std::endl << std::endl

            << "}" << std::endl // end of namespace LocARNA
            << "#endif" << std::endl;
    }


    double
    RibosumFreq::basematch_score_corrected(char i,char j) const {

        double background_i=base_nonstruct_prob(i);
        double background_j=base_nonstruct_prob(j);

        return log( basematch_prob(i,j) / (background_i * background_j) ) / log(2);
    }

    void
    RibosumFreq::print_basematch_scores_corrected(std::ostream &out) const {
        for (Alphabet<char>::const_iterator it=char_basename_alphabet.begin(); char_basename_alphabet.end()!=it; ++it) {
            out << (*it) <<" ";
            for (Alphabet<char>::const_iterator it2=char_basename_alphabet.begin(); char_basename_alphabet.end()!=it2; ++it2) {
                out << basematch_score_corrected(*it,*it2)<<" ";
            }
            out << std::endl;
        }
        out << std::endl;
    }

}
