#include "fitonoff.hh"

#include <iostream>
#include <limits>

namespace LocARNA {

    //! square difference
    double sqdiff(double x,double y) {
	return (x-y)*(x-y);
    }


    double
    FitOnOff::viterbi(double c0, double c1, bool traceback) {
	v[0][0]=0;
	v[1][0]=0;
    
	for (size_type i=1;i<=x.size(); ++i) {
	    v[0][i] = sqdiff(x[i-1],c0) + std::min(v[0][i-1],v[1][i-1]+delta_10);
	    v[1][i] = sqdiff(x[i-1],c1) + std::min(v[1][i-1],v[0][i-1]+delta_01);
	
	    if (traceback) {
		t[0][i] = (v[0][i-1] > v[1][i-1]+delta_10); // 0 for going to 0, 1 for going to 1
		t[1][i] = !(v[1][i-1] > v[0][i-1]+delta_01);
	    }
	}

	if (traceback) {
	    bool state = v[1][x.size()]<v[0][x.size()]; // 1 if v[1] end smaller, 0 otherwise
	
	    for (size_type i=x.size() ;i>0; --i) {
		trace[i]=state;
		state=t[state][i];
	    }
	}
    
	return std::min(v[0][x.size()],v[1][x.size()]);
    }

    double FitOnOff::best_once_on(double c0, double c1) {

	size_type min_size=std::min((size_t)20,x.size()); // on signal at least of size 20!
    
	std::vector<double> sum_dist_on;
	std::vector<double> sum_dist_off;
	sum_dist_on.resize(x.size()+1);
	sum_dist_off.resize(x.size()+1);

	sum_dist_on[0]=0;
	sum_dist_off[0]=0;
	for (size_type i=1; i<=x.size(); ++i) {
	    sum_dist_on[i]  = sum_dist_on[i-1] +sqdiff(x[i-1],c1);
	    sum_dist_off[i] = sum_dist_off[i-1]+sqdiff(x[i-1],c0);
	}
    
	size_type best_i=1;
	size_type best_j=1;
	double best_dist=std::numeric_limits<double>::max();
    
	for (size_type i=1; i<=x.size(); ++i) {
	    for (size_type j=i+min_size-1; j<=x.size(); ++j) {
		// compute distance for signal that is on in interval [i..j] and off otherwise.
		// This distance is \sum_{k=1..i-1} sqdiff(x[k],c0) + \sum_{k=i..j} sqdiff(x[k],c1) + \sum_{k=i+1..x.size()} sqdiff(x[k],c0)
		// and computed from sum_dist_off and sum_dist_on in constant time.
		double dist = sum_dist_off[i-1] + ( sum_dist_on[j]-sum_dist_on[i-1] )  + sum_dist_off[x.size()]-sum_dist_off[j+1];
	    
		if (dist < best_dist) {
		    best_dist=dist;
		    best_i=i;
		    best_j=j;
		}
	    }
	}

	// fill trace
	size_type i;
	for (i=1; i<best_i; ++i) {trace[i]=0;}
	for (; i<=best_j; ++i) {trace[i]=1;}
	for (; i<=x.size(); ++i) {trace[i]=0;}
    
	return best_dist;
    }



    pf_t
    FitOnOff::forward(double c0, double c1) {
	v[0][0]=1;
	v[1][0]=1;
    
	for (size_type i=1;i<=x.size(); ++i) {
	    v[0][i] = exp(-beta*sqdiff(x[i-1],c0)) * (v[0][i-1] + v[1][i-1]*exp_delta_10);
	    v[1][i] = exp(-beta*sqdiff(x[i-1],c1)) * (v[1][i-1] + v[0][i-1]*exp_delta_01);
	}

	return v[0][x.size()] + v[1][x.size()];
    }

    std::pair<double,double>
    FitOnOff::optimize(double c0, double c1) {

	double stepwidth=0.0025; // length of the gradient vector added in one iteration
	size_type max_steps=400; // number of iterations


	bool converged=false;
    
	std::vector<pf_t> theta; // factor for staying in numeric range
	theta.resize(x.size()+1);
    
	std::vector<std::vector<std::vector<pf_t> > > dv; // partial derivations dv[x][y][i] of v[x][i] by cy
	dv.resize(2);
	dv[0].resize(2);
	dv[1].resize(2);
	dv[0][0].resize(x.size()+1);
	dv[0][1].resize(x.size()+1);
	dv[1][0].resize(x.size()+1);
	dv[1][1].resize(x.size()+1);
    
	pf_t last2_c0=c1;
	pf_t last2_c1=c0;
	pf_t last_c0=c1;
	pf_t last_c1=c0;
    

	// perform iterations of optimization
	// until convergence or maximal number of iterations
	size_type step=0;
	while (!converged && step<max_steps) {
	    step++;

	
	    // ----------------------------------------
	    // compute thetas for c0,c1
	    viterbi(c0,c1,false);
	
	    pf_t Theta=1;
	    theta[0]=1;
	    for (size_type i=1;i<=x.size(); ++i) {
		//theta[i]=1;
		// estimate factor from difference between sum of viterbi scores at i and i-1
		theta[i]=exp(beta*(v[0][i]+v[1][i]-v[0][i-1]-v[1][i-1])); // exp(-beta*viterbi_score/x.size());
		Theta *= theta[i];
	    }
	
	
	    // determine gradient
	
	    // initialize
	    v[0][0]=1;
	    v[1][0]=1;
	
	    dv[0][0][0] = 0;
	    dv[0][1][0] = 0;
	    dv[1][0][0] = 0;
	    dv[1][1][0] = 0;
	
	    // recurse
	    for (size_type i=1;i<=x.size(); ++i) {
		pf_t exp_dev0=exp(-beta*sqdiff(x[i-1],c0)); // boltzmann weight of square devitation xi to c0
		pf_t exp_dev1=exp(-beta*sqdiff(x[i-1],c1)); // boltzmann weight of square devitation xi to c1
	    
		v[0][i] = theta[i]*exp_dev0 * (v[0][i-1] + v[1][i-1]*exp_delta_10);
		v[1][i] = theta[i]*exp_dev1 * (v[1][i-1] + v[0][i-1]*exp_delta_01);
    	    
		dv[0][0][i] = theta[i]*
		    exp_dev0 * ( dv[0][0][i-1] + dv[1][0][i-1]*exp_delta_10
				 +
				 (v[0][i-1]+v[1][i-1]*exp_delta_10)*2*beta*(x[i-1]-c0));
	    
		dv[0][1][i] = theta[i] * exp_dev0 * ( dv[0][1][i-1] + dv[1][1][i-1]*exp_delta_10);
	    
		dv[1][0][i] = theta[i] * exp_dev1 * ( dv[1][0][i-1] + dv[0][0][i-1]*exp_delta_01);
	    
		dv[1][1][i] = theta[i] *
		    exp_dev1 * ( dv[1][1][i-1] + dv[0][1][i-1]*exp_delta_01
				 +
				 (v[1][i-1]+v[0][i-1]*exp_delta_01)*2*beta*(x[i-1]-c1));
	    }
	
	    //pf_t Z = v[0][x.size()]+v[1][x.size()];
	
	
	    pf_t d0 = dv[0][0][x.size()] + dv[1][0][x.size()];	
	    pf_t d1 = dv[0][1][x.size()] + dv[1][1][x.size()];
	
	
	    // normalize vector (d0,d1) in a numerically nice way
	    pf_t dlen_by_d0 = sqrt(1+d1/d0*d1/d0); // = dlen/d0
	    pf_t dlen_by_d1 = sqrt(d0/d1*d0/d1+1); // = dlen/d1
	
	    d0=(d0>0?1:-1)*stepwidth*1/dlen_by_d0;
	    d1=(d1>0?1:-1)*stepwidth*1/dlen_by_d1;
	
	
	    c0+=d0;
	    c1+=d1;
	
	
	    converged=fabs(last2_c0-c0)<stepwidth/2 && fabs(last2_c1-c1)<stepwidth/2;
	    last2_c0=last_c0;
	    last2_c1=last_c1;
	    last_c0=c0;
	    last_c1=c1;
	}
    
	return std::pair<double,double>(c0,c1);
    }


    void
    FitOnOff::write_viterbi_path(std::ostream &out,double c0, double c1) const {
	out.precision(2);
	for (size_type i=1;i<=x.size(); ++i) {
	    out << (trace[i]==0?c0:c1)<<std::endl;
	}
	//out<<std::endl;
    }

    void
    FitOnOff::write_viterbi_path_compact(std::ostream &out,double c0, double c1) {
    
	const int on=(c0<c1)?1:0;
	const int off=1-on;

	trace[0]=off;

	// for (size_type i=0;i<=x.size(); ++i) {
	// 	out << i<<":"<<trace[i] << " ";
	// }
	// out <<std::endl;
    
    
	for (size_type i=1;i<=x.size(); ++i) {
	    if (trace[i-1]!=trace[i]) {
		if (trace[i]==on) {
		    out << i << " ";
		} else {
		    out << (i-1) << " ";
		}
	    }
	}
    
	if (trace[x.size()]!=trace[0]) {
	    out << x.size();
	}
    
	out<<std::endl;
    }



    void
    FitOnOff::print_table(const std::string &name, const std::vector<bool> &v) const {
	std::cout <<name<<" ";
	for (size_type i=1;i<=x.size(); ++i) {std::cout << v[i] <<"  ";}
	std::cout <<std::endl;
    }
    void
    FitOnOff::print_table(const std::string &name, const std::vector<pf_t> &v) const {
	std::cout.precision(2);
	std::cout <<name<<" ";
	for (size_type i=1;i<=x.size(); ++i) {std::cout << v[i] <<" ";}
	std::cout <<std::endl;
    }


    void
    FitOnOff::print_tables() const 
    {
	print_table("v0",v[0]);
	print_table("v1",v[1]);
	print_table("t0",t[0]);
	print_table("t1",t[1]);
	print_table("tr",trace);
    }

} // END namespace LocARNA
