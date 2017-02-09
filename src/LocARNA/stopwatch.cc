#include "stopwatch.hh"
#include <sys/time.h>
#include <assert.h>
#include <iostream>
#include <iomanip>



namespace LocARNA {
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
        size_t oldwidth=out.width();
        std::ios_base::fmtflags oldfmt = out.flags();

        out << " "
            << std::setw(14) << std::left << name
            << " ";

        std::streamsize oldprecision = out.precision(3);

        out
            << std::setw(8) << std::right << std::fixed << current_total(name)
            << "s ("
            <<current_cycles(name)
            <<" cycles)"  << std::endl;

        out.width(oldwidth);
        out.precision(oldprecision);
        out.setf(oldfmt);
        out.unsetf(std::ios_base::fixed);

        return out;
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
