#include "aux.hh"
#include <iomanip>

namespace LocARNA {
    

    failure::~failure() throw() {};
    
    const char *
    failure::what() const throw() {
	//exception::what();
	return msg_.c_str();
    }

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
	for (size_type i=0; i<=seq.length(); i++) {
	    if (seq[i]=='T') seq[i]='U';
	}
    }

    // global stop watch
    StopWatch stopwatch;
    // ------------------------------------------------------------
    // implement StopWatch
    
    StopWatch::StopWatch(bool print_on_exit_): print_on_exit(print_on_exit_) {
    }

    StopWatch::~StopWatch() {
	if (print_on_exit) {
	    print_info(std::cerr);
	}
    }

    void
    StopWatch::set_print_on_exit(bool print_on_exit_) {
	print_on_exit=print_on_exit_;
    }
    
    bool
    StopWatch::start(const std::string &name) {
	timer_t &t=timers[name];
	
	if (t.running) return false;
	
	t.last_start=current_time();
	t.running=true;
	
	return true;
    }

    bool
    StopWatch::stop(const std::string &name) {
	assert(timers.find(name)!=timers.end());
	
	timer_t &t=timers[name];
	
	if (!t.running) return false; //allow stop without start
	
	t.cycles++;
	t.total += current_time() - t.last_start;
	t.running=false;
	
	return true;
    }

    bool
    StopWatch::is_running(const std::string &name) const {
	map_t::const_iterator it = timers.find(name);
	assert(it!=timers.end());
	const timer_t &t=it->second;
	return t.running;
    }

    double
    StopWatch::current_total(const std::string &name) const {
	map_t::const_iterator it = timers.find(name);
	assert(it!=timers.end());
	const timer_t &t=it->second;
	return
	    t.total
	    +
	    (t.running?current_time()-t.last_start:0.0);
    }
	
    size_t StopWatch::current_cycles(const std::string &name) const {
	map_t::const_iterator it = timers.find(name);
	assert(it!=timers.end());
	const timer_t &t=it->second;
	return
	    t.cycles
	    +
	    (t.running?1:0);
    }

    std::ostream &
    StopWatch::print_info(std::ostream &out,const std::string &name) const {
	return
	    out << " " 
		<< std::setw(14) << std::left <<name 
		<< " "
		<< std::setw(8) << std::right << std::fixed << std::setprecision(2) << current_total(name) 
		<< "s ("
		<<current_cycles(name)
		<<" cycles)" << std::left << std::endl;
    }

    std::ostream &
    StopWatch::print_info(std::ostream &out) const {
	
	if (timers.size()==0) return out;
	
	out << "------------------------------"<<std::endl;
	out << "Stopped Times"<<std::endl;
	
	for (map_t::const_iterator it=timers.begin(); timers.end()!=it; ++it) {
	    print_info(out,it->first);
	}
	return out;
    }

    double 
    StopWatch::current_time () const { //returns current time in seconds
	timeval tv;
	gettimeofday(&tv, NULL);
	double rtn_value = (double) tv.tv_usec;
	rtn_value /= 1e6;
	rtn_value += (double) tv.tv_sec;
	return rtn_value;
    }	
}

