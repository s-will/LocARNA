#include "ribosum.hh"
#include "aux.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

namespace LocARNA {

    Ribosum::Ribosum(const std::string &filename)
        : name_(),
          bm_(),
          am_(),
          basename_alphabet_(),
          arcname_alphabet_(),
          char_basename_alphabet_() {
        std::ifstream in(filename.c_str());
        if (!in) {
            std::ostringstream err;
            err << "Cannot open file " << filename
                << " for reading ribosum data.";
            throw failure(err.str());
        }

        read_ribosum(in);

        in.close();
    }

    Ribosum::Ribosum()
        : name_(),
          bm_(),
          am_(),
          basename_alphabet_(),
          arcname_alphabet_(),
          char_basename_alphabet_() {}

    Ribosum::~Ribosum() {}

    void
    Ribosum::read_ribosum(std::istream &in) {
        try {
            std::string line;

            if (!std::getline(in, line))
                throw(std::ifstream::failure("Expecting name."));
            name_ = line;

            // ------------------------------------------------------------
            // read the 4x4 matrix

            std::array<std::string, 4> basenames = {"A", "C", "G", "U"};
            basename_alphabet_ = base_alphabet_type(basenames);
            char_basename_alphabet_ = make_char_alphabet();

            read_matrix(in, bm_, basename_alphabet_);

            // ignore two lines -- not clean, but we don't need the info
            // currently
            getline(in, line);
            getline(in, line);

            // ------------------------------------------------------------
            // read the 16x16 matrix

            std::vector<std::string> arcnames;
            for (const auto &i : basenames)
                for (const auto &j : basenames)
                    arcnames.push_back(i + j);
            arcname_alphabet_ = arc_alphabet_type(arcnames);

            read_matrix(in, am_, arcname_alphabet_);

            // ignore two lines
            getline(in, line);
            getline(in, line);

        } catch (std::ifstream::failure &e) {
            std::ostringstream err;
            err << "Cannot parse ribosum input. " << e.what() << std::endl
                << "File not in ribosum format.";
            throw failure(err.str());
        }
    }

    Ribosum::char_alphabet_type
    Ribosum::make_char_alphabet() const {
        // transform alphabet over singleton strings into one over characters
        std::array<char, 4> alphchars;
        for (size_type i=0; i<4; ++i) {
            alphchars[i] = basename_alphabet_[i][0];
        }
        return char_alphabet_type(alphchars);
    }

    /**
     * Output operator
     *
     * @param out output stream
     * @param ribosum Ribosum to be written to stream
     *
     * @return output stream after writing
     */
    std::ostream &
    operator<<(std::ostream &out, const Ribosum &ribosum) {
        out << ribosum.name_ << std::endl << std::endl;

        ribosum.write_matrix(out, ribosum.bm_, ribosum.basename_alphabet_);
        ribosum.write_matrix(out, ribosum.am_, ribosum.arcname_alphabet_);

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
          arcmatch_probs_() {
        std::ifstream in(filename.c_str());
        if (!in.is_open()) {
            std::cerr << "Cannot open file " << filename
                      << " for reading ribosum data." << std::endl;
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
          arcmatch_probs_() {}

    /**
     * Test for only blank characters
     *
     * @param s string to be tested
     *
     * @return whether s consists of only blank characters
     */
    bool
    is_blank(std::string &s) {
        for (size_t i = 0; i < s.length(); i++)
            if (s[i] != ' ')
                return false;
        return true;
    }

    void
    RibosumFreq::read_frequencies(std::istream &in) {
        read_matrix(in, "BASE FREQUENCIES", base_probs_, 4, 1);
        read_matrix(in, "BASE NONSTRUCTURAL FREQUENCIES", base_nonstruct_probs_,
                    4, 1);
        read_matrix(in, "BASE PAIR FREQUENCIES", basepair_probs_, 4, 4);
        read_matrix(in, "BASE MATCH FREQUENCIES", basematch_probs_, 4, 4);
        read_matrix(in, "BASE PAIR MATCH FREQUENCIES", arcmatch_probs_, 16, 16);
    }

    void
    RibosumFreq::read_matrix(std::istream &in,
                             const std::string &header,
                             matrix_t &mat,
                             size_t xdim,
                             size_t ydim) {
        try {
            std::string line;
            while (std::getline(in, line) && is_blank(line))
                ;

            if (line != header) {
                throw std::ifstream::failure("Expected header " + header + "." +
                                             " Read instead '" + line + "'.");
            }

            mat.resize(xdim, ydim);

            for (size_t i = 0; i < xdim; i++)
                for (size_t j = 0; j < ydim; j++)
                    in >> mat(i, j);

        } catch (std::ifstream::failure &e) {
            std::ostringstream err;
            err << "Cannot parse ribosum frequency input. " << e.what()
                << std::endl
                << "File not in extended ribosum format.";
            throw failure(err.str());
        }
    }

    std::ostream &
    RibosumFreq::write_matrix(std::ostream &out,
                              const std::string &name,
                              const Matrix<double> &mat) const {
        out << name << std::endl;
        out << mat << std::endl;
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
    std::ostream &
    operator<<(std::ostream &out, const RibosumFreq &ribosum) {
        out << (Ribosum)(ribosum);
        out << std::endl;

        ribosum.write_matrix(out, "BASE FREQUENCIES", ribosum.get_base_probs());
        ribosum.write_matrix(out, "BASE NONSTRUCTURAL FREQUENCIES",
                             ribosum.get_base_nonstruct_probs());
        ribosum.write_matrix(out, "BASE PAIR FREQUENCIES",
                             ribosum.get_basepair_probs());
        ribosum.write_matrix(out, "BASE MATCH FREQUENCIES",
                             ribosum.get_basematch_probs());
        ribosum.write_matrix(out, "ARC MATCH FREQUENCIES",
                             ribosum.get_arcmatch_probs());

        return out;
    }

    void
    RibosumFreq::write_CC_matrix(std::ostream &out,
                                 const std::string &ribname,
                                 const std::string &matname,
                                 int x,
                                 int y,
                                 const Ribosum::matrix_t &m) const {
        out << "const double " << ribname << "::" << matname << "[] = {"
            << std::endl;
        for (int i = 0; i < x; i++) {
            out << "    ";
            for (int j = 0; j < y; j++) {
                out << m(i, j);
                if (i < (x - 1) || j < (y - 1))
                    out << ", ";
                else
                    out << " ";
            }
            out << std::endl;
        }
        out << "};" << std::endl << std::endl;
    }

    void
    RibosumFreq::write_ICC_code(std::ostream &out,
                                const std::string &ribname) const {
        std::string macro_headerfile = ribname;
        transform(macro_headerfile.begin(), macro_headerfile.end(),
                  macro_headerfile.begin(), ::toupper);

        macro_headerfile = "LOCARNA_" + macro_headerfile + "_HH";

        out << "#ifndef " << macro_headerfile << std::endl
            << "#define " << macro_headerfile << std::endl
            << "#include \"LocARNA/ribosum.hh\"" << std::endl
            << "namespace LocARNA {" << std::endl
            << "class " << ribname << ": public RibosumFreq {" << std::endl
            << "    static const double bm_init[];" << std::endl
            << "    static const double am_init[];" << std::endl
            << "    static const double base_probs_init[];" << std::endl
            << "    static const double base_nonstruct_probs_init[];"
            << std::endl
            << "    static const double basepair_probs_init[];" << std::endl
            << "    static const double basematch_probs_init[];" << std::endl
            << "    static const double arcmatch_probs_init[];" << std::endl
            << std::endl
            << "  public:" << std::endl
            << "    " << ribname << "(): RibosumFreq() {" << std::endl
            << "        bm_ = matrix_t(4,4,bm_init);" << std::endl
            << "        am_ = matrix_t(16,16,am_init);" << std::endl

            << "        base_probs_ = matrix_t(4,1,base_probs_init);"
            << std::endl
            << "        base_nonstruct_probs_ = "
               "matrix_t(4,1,base_nonstruct_probs_init);"
            << std::endl
            << "        basepair_probs_ = matrix_t(4,4,basepair_probs_init);"
            << std::endl
            << "        basematch_probs_ = matrix_t(4,4,basematch_probs_init);"
            << std::endl
            << "        arcmatch_probs_ = matrix_t(16,16,arcmatch_probs_init);"
            << std::endl
            << "        set_basename_alphabet({\"A\",\"C\",\"G\",\"U\"});"
            << std::endl
            << "        set_arcname_alphabet({\"AA\",\"AC\",\"AG\",\"AU\",\"CA\",\"CC\",\"CG\",\"CU\",\"GA\",\"GC\",\"GG\",\"GU\",\"UA\",\"UC\",\"UG\",\"UU\"});"
            << std::endl
            << "    }" << std::endl
            << "};" << std::endl;

        write_CC_matrix(out, ribname, "bm_init", 4, 4, get_basematch_scores());
        write_CC_matrix(out, ribname, "am_init", 16, 16, get_arcmatch_scores());

        write_CC_matrix(out, ribname, "base_probs_init", 4, 1,
                        get_base_probs());
        write_CC_matrix(out, ribname, "base_nonstruct_probs_init", 4, 1,
                        get_base_nonstruct_probs());

        write_CC_matrix(out, ribname, "basepair_probs_init", 4, 4,
                        get_basepair_probs());
        write_CC_matrix(out, ribname, "basematch_probs_init", 4, 4,
                        get_basematch_probs());

        write_CC_matrix(out, ribname, "arcmatch_probs_init", 16, 16,
                        get_arcmatch_probs());

        out << "}" << std::endl // end of namespace LocARNA
            << "#endif" << std::endl;
    }

    double
    RibosumFreq::basematch_score_corrected(char i, char j) const {
        double background_i = base_nonstruct_prob(i);
        double background_j = base_nonstruct_prob(j);

        return log(basematch_prob(i, j) / (background_i * background_j)) /
            log(2);
    }

    void
    RibosumFreq::print_basematch_scores_corrected(std::ostream &out) const {
        for (const auto &x: char_basename_alphabet_) {
            out << x << " ";
            for (const auto &y : char_basename_alphabet_) {
                out << basematch_score_corrected(x,y) << " ";
            }
            out << std::endl;
        }
        out << std::endl;
    }
}
