#ifndef LOCARNA_MATRICES_HH
#define LOCARNA_MATRICES_HH

#include <iostream>
#include <vector>
#include <assert.h>

#include <algorithm>

namespace LocARNA {


/*
  Define classes for the dynamic programming
  matrices.

  The structures support offsets and moving
  of the offset, while maintaining the entries
  in the overlapping sub-matrix
*/

//! simple 2D matrix class, provides access via operator (int,int) 
template <class T>
class Matrix {
public:
    typedef T elem_t;
    typedef typename std::vector<elem_t>::size_type size_type;
    
    typedef std::pair<size_type,size_type> size_pair_type;
    
protected:
    std::vector<elem_t> mat_;
    size_type xdim_;
    size_type ydim_;
    
    size_type addr(size_type i, size_type j) const {
	assert(0<=i && i<this->xdim_);
	assert(0<=j && j<this->ydim_);
	return i*ydim_+j;
    }

public:
    Matrix() 
	: mat_(),xdim_(0),ydim_(0) {
    }
    
    // if from given and !=0 initialize from array <from>
    Matrix(size_type xdim, size_type ydim, const elem_t *from=0L)
	: mat_(xdim*ydim),xdim_(xdim),ydim_(ydim) {
	if (from!=0L) {
	    for (size_type i=0; i<xdim_; i++) {
		for (size_type j=0; j<ydim_; j++) {
		    (*this)(i,j)=from[i*ydim+j];
		}
	    }
	}
    }
    
    size_pair_type sizes() const {
	return size_pair_type(xdim_,ydim_);
    }

    //! resize both dimensions
    void
    resize(size_type xdim, size_type ydim) {
	xdim_=xdim;
	ydim_=ydim;
	
	mat_.resize(xdim_*ydim_);
    }
    
    //! read only access to (i,j)
    const elem_t & operator() (size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }
    
    //! read/write access to (i,j)
    elem_t & operator() (size_type i,size_type j) {
	return mat_[addr(i,j)];
    }

    //! read only access to (i,j)
    const elem_t get(size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }

    //! write access to (i,j)
    void
    set(size_type i,size_type j, const elem_t &x) {
	mat_[addr(i,j)]=x;
    }

    
    //! fill the whole matrix with the given value 
    void 
    fill(const elem_t &val) {
	for (size_type i=0; i<xdim_*ydim_; ++i)
	    mat_[i]=val;
    }

    void
    clear() {
	resize(0,0);
	mat_.clear();
    }
    
    template<class UnaryOperator>
    void transform(UnaryOperator f) {
	std::transform(mat_.begin(),mat_.end(),mat_.begin(),f);
    }
};

template <class T>
std::ostream & operator << (std::ostream &out, Matrix<T> mat) {
    typename Matrix<T>::size_pair_type sizes = mat.sizes();
    
    for (typename Matrix<T>::size_type i=0; i<sizes.first; i++) {
	    for (typename Matrix<T>::size_type j=0; j<sizes.second; j++) {
		out << mat(i,j) << " ";
	    }
	    out << std::endl;
    }
    return out;
}

template <class T>
std::istream & operator >> (std::istream &in, Matrix<T> mat) {
    typename Matrix<T>::size_pair_type sizes = mat.sizes();
    for (typename Matrix<T>::size_type i=0; i<=mat.sizes().first; i++) {
	for (typename Matrix<T>::size_type j=0; j<=mat.sizes().second; j++) {
	    in >> mat(i,j);
	}
    }
    return in;
}


// ----------------------------------------
//! simple matrix class with range.
//! The matrix offers a fix maximal range [0..n]x[0..m].
//! It can be restricted to [xl..xr]x[yl..yr]. After such
//! a restriction, the matrix is invalidated and can
//! only be used with indices (i,j): xl<=i<=xr and yl<=j<=yr  
//
// it was planned to use this for the M matrices in LocARNA,
// I didn't see a performance improvement (maybe for very large instances?)

template <class elem_t>
class RMatrix : public Matrix<elem_t> {
    typedef typename Matrix<elem_t>::size_type size_type;
protected:
    size_type xl_;
    size_type xr_;
    size_type yl_;
    size_type yr_;

    size_type xfactor_;
    size_type offset_;
    
    
    size_type addr(size_type i, size_type j) const {
	assert(xl_<=i && i<=xr_);
	assert(yl_<=j && j<=yr_);
	
	return i*xfactor_ + j - offset_;
    }
    
public:
    RMatrix() : Matrix<elem_t>() {
    }
    
    //! reserves memory in sufficient size
    //! removes any existing restrictions
    void
    resize(size_type xdim, size_type ydim) {
	this->xdim_=xdim;
	this->ydim_=ydim;
	this->mat_.resize(xdim*ydim);
	
	restrict(0,xdim-1,0,ydim-1);
    }

    //! set new restricted range
    void restrict(size_type xl,size_type xr,size_type yl,size_type yr) {
	assert(xl>=0);
	assert(yl>=0);
	assert(xr<this->xdim_);
	assert(yr<this->ydim_);
	assert(xl<=xr);
	assert(yl<=yr);
	
	this->xl_=xl;
	this->xr_=xr;
	this->yl_=yl;
	this->yr_=yr;

	this->xfactor_ = (yr-yl)+1;
	this->offset_  = 0;
	this->offset_  = addr(xl,yl);
    }
    
    const elem_t & operator() (size_type i,size_type j) const {
	return this->mat_[addr(i,j)];
    }

    elem_t & operator() (size_type i,size_type j) {
	return this->mat_[addr(i,j)];
    }

    //! read only access to (i,j)
    const elem_t get(size_type i,size_type j) const {
	return this->mat_[addr(i,j)];
    }

    //! write access to (i,j)
    void
    set(size_type i,size_type j, const elem_t &x) {
	this->mat_[addr(i,j)]=x;
    }

    //! fill the restricted area of the matrix with the given value
    void 
    fill(const elem_t &val) {
	for (size_type i=xl_; i<xr_; ++i)
	    for (size_type j=yl_; j<yr_; ++j)
		mat_(i,j)=val;
    }
};


// ----------------------------------------
//! simple matrix class with offset
template <class elem_t>
class OMatrix : public Matrix<elem_t> {
protected:
    size_t off_;
    size_t xoff_;
    size_t yoff_;

    size_t addr(size_t i, size_t j) const {
	assert(xoff_<=i && i<xoff_+this->xdim_);
	assert(yoff_<=j && j<yoff_+this->ydim_);
	return i*this->xdim_ + j - off_;
    }
    
public:
    OMatrix() : Matrix<elem_t>(0) {
    }
        
    void
    resize(size_t xdim, size_t ydim, size_t xoff=0, size_t yoff=0) {
	xoff_=xoff;
	yoff_=yoff;
	off_=xoff*xdim+yoff;
	this->xdim_=xdim;
	this->ydim_=ydim;
	this->mat_.resize(xdim*ydim);
    }

    const elem_t & operator() (size_t i,size_t j) const {
	return this->mat_[addr(i,j)];
    }

    elem_t & operator() (size_t i,size_t j) {
	return this->mat_[addr(i,j)];
    }

    //! read only access to (i,j)
    const elem_t get(size_t i,size_t j) const {
	return this->mat_[addr(i,j)];
    }

    //! write access to (i,j)
    void
    set(size_t i,size_t j, const elem_t &x) {
	this->mat_[addr(i,j)]=x;
    }

};


// ----------------------------------------
//! matrix class with rotation
template <class elem_t>
class RotMatrix: public Matrix<elem_t> {

protected:
    size_t xrot_;
    size_t yrot_;

    size_t rot(size_t x, size_t r, size_t d) {
	assert(r<d);
	return (x+d-r)%d;
    }
    
    size_t addr(size_t i, size_t j) const {
	assert(xrot_<=i && i<xrot_+this->xdim_);
	assert(yrot_<=j && j<yrot_+this->ydim_);
	return rot(i,xrot_,this->xdim_)*this->xdim_ + rot(j,yrot_,this->ydim_);
    }

public:    
    RotMatrix() : Matrix<elem_t>(0) {
    }
        
    void
    resize(size_t xdim, size_t ydim, size_t xrot=0, size_t yrot=0) {
	xrot_=xrot;
	yrot_=yrot;
	this->xdim_=xdim;
	this->ydim_=ydim;
	this->mat_.resize(xdim*ydim);
    }

    void move(size_t xrot, size_t yrot) {
	xrot=xrot_;
	yrot=yrot_;
    }
    
    const elem_t & operator() (size_t i,size_t j) const {
	return this->mat_[addr(i,j)];
    }

    elem_t & operator() (size_t i,size_t j) {
	return this->mat_[addr(i,j)];
    }

    //! read only access to (i,j)
    const elem_t get(size_t i,size_t j) const {
	return this->mat_[addr(i,j)];
    }

    //! write access to (i,j)
    void
    set(size_t i,size_t j, const elem_t &x) {
	this->mat_[addr(i,j)]=x;
    }
};





} // end namespace LocARNA

#endif // LOCARNA_MATRICES_HH
