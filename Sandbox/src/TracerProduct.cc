//
// A test class that makes printout whenever its methods are called.
//
// $Id: TracerProduct.cc,v 1.1 2011/06/04 20:36:59 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/04 20:36:59 $
//
// Original author Rob Kutschke
//

#include "Sandbox/inc/TracerProduct.hh"

#include <iostream>

using namespace std;

namespace mu2e {

  TracerProduct::TracerProduct():
    val_(-1),
    serial_(count()){
    cerr << "TracerProduct default c'tor: " 
	 << val_    << " " 
	 << serial_ << endl;
  }

  TracerProduct::TracerProduct( int aval):
    val_(aval),
    serial_(count()){
    cerr << "TracerProduct  c'tor from int: " 
	 << val_    << " " 
	 << serial_ << endl;
  }

  TracerProduct::TracerProduct( TracerProduct const& rhs):
    val_(rhs.val_),
    serial_(count()){
    cerr << "TracerProduct copy c'tor: " 
	 << val_    << " " 
	 << serial_ << " from: "
	 << rhs.serial_
	 << endl;
  }

  TracerProduct::~TracerProduct(){
    cerr << "TracerProduct d'tor: " << val_ << " " << serial_ << endl;
  }

  TracerProduct& TracerProduct::operator=( TracerProduct const& rhs ){
    cerr << "TracerProduct: operator= : "
	 << "( " << serial_     << " " << val_    << " ),"
	 << "( " << rhs.serial_ << " " << rhs.val_ << " )"
	 << endl;
    val_ = rhs.val_;
    return *this;
  }

  bool TracerProduct::operator<( TracerProduct const&  rhs ){
    cerr << "TracerProduct: operator< : "
	 << "( " << serial_     << " " << val_    << " ),"
	 << "( " << rhs.serial_ << " " << rhs.val_ << " )"
	 << endl;
    bool answer(val_ < rhs.val_);
    return answer;
  } 

  int TracerProduct::count(){
    static int c(0);
    return c++;
  }

}  // namespace mu2e
