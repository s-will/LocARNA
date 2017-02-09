#include "aux.hh"
#include "string1.hh"
#include "matrix.hh"

#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

namespace LocARNA {
    failure::~failure() throw() {}

    const char *
    failure::what() const throw() {
        //exception::what();
        return msg_.c_str();
    }

    // ------------------------------------------------------------
    // gaps

    // gap symbols
    const std::string the_gap_symbols = "-_~.";
    // simplified gap symbols
    const std::string simplified_gap_symbols = "---.";


    size_t Gap::size=0;
    const Gap Gap::regular  = Gap(Gap::size++);
    const Gap Gap::loop     = Gap(Gap::size++);
    const Gap Gap::locality = Gap(Gap::size++);
    const Gap Gap::other    = Gap(Gap::size++);

    bool
    is_gap_symbol(char c) {
        return the_gap_symbols.find(c)!=std::string::npos;
    }

    char
    gap_symbol(Gap gap) {
        return simplified_gap_symbols[gap.idx()];
    }

    char
    special_gap_symbol(Gap gap) {
        return the_gap_symbols[gap.idx()];
    }

    Gap gap_code(char symbol) {
        assert(is_gap_symbol(symbol));
        return Gap(the_gap_symbols.find(symbol));
    }

    // ------------------------------------------------------------

    /**
     * \brief Converts char to upper case
     *
     * Function class used by transform_toupper().
     */
    class ToUpper {
    public:
        //! convert to upper case
        char operator() (char c) const  { return std::toupper(c); }
    };

    void transform_toupper(std::string &s) {
        std::transform(s.begin(),s.end(),s.begin(),ToUpper());
    }

    void
    normalize_rna_sequence(std::string &seq) {
        transform_toupper(seq);
        for (size_type i=0; i<seq.length(); i++) {
            if (seq[i]=='T') seq[i]='U';
        }
    }

    bool
    has_prefix(const std::string &s, const std::string &p, size_t start) {
        if (s.length()<p.length()-start) {
            return false;
        }
        return s.substr(start,p.length())==p;
    }


    void
    split_at_separator(const std::string &s, char sep, std::vector<std::string> &v) {
        std::string str=s;
        v.clear();
        size_t pos;
        while ((pos=str.find(sep))!=std::string::npos) {
            if (pos>0) {
                v.push_back(str.substr(0,pos));
            } else {
                v.push_back("");
            }
            str = str.substr(pos+1); // note: if pos+1 == length of str, substr yields empty
        }
        v.push_back(str);
    }

    std::vector<std::string>
    split_at_separator(const std::string &s, char sep) {
        std::vector<std::string> v;
        split_at_separator(s,sep,v);
        return v;
    }


    std::string
    concat_with_separator(const std::vector<std::string> &v, char sep) {
        if (v.size()==0) return "";
        std::string s=v[0];
        for (std::vector<std::string>::const_iterator it=v.begin()+1; v.end()!=it; ++it) {
            s += sep + *it;
        }
        return s;
    }

    bool
    get_nonempty_line(std::istream &in,
                               std::string &line) {
        const std::string whitespace = " \t";
        while(getline(in,line)) {
            if (line.length()>0 &&
                !isspace(line[0])) {

                //  collect lines in case of quoted newline
                while (line[line.length()-1]=='\\') {
                    line=line.substr(0,line.length()-1);
                    std::string line1;
                    if (getline(in,line1)) {
                        line+=" "+line1;
                    } else {
                        break;
                    }
                }

                // remove trailing white space
                const size_t strEnd = line.find_last_not_of(whitespace);
                line = line.substr(0, strEnd+1);
                return true;
            }
        }
        line="";
        return false;
    }

    double
    sequence_identity(const string1 &seqA, const string1 &seqB) {
        size_t n=seqA.length();
        size_t m=seqB.length();

        Matrix<size_t> D(n+1,m+1,0);
        for (size_t i=1; i<=n;i++) {
            for (size_t j=1; j<=m;j++) {
                D(i,j) = std::max(D(i-1,j),D(i,j-1));
                D(i,j) = std::max(D(i,j),
                                  D(i-1,j-1) + ((seqA[i]==seqB[j])?1:0)
                                  );
            }
        }
        return 100 * D(n,m) / (double)std::min(n,m);
    }
}

