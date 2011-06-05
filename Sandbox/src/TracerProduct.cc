//
// A test class that makes printout whenever its methods are called.
//
// $Id: TracerProduct.cc,v 1.2 2011/06/05 16:13:46 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/05 16:13:46 $
//
// Original author Rob Kutschke
//

#include "Sandbox/inc/TracerProduct.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

namespace mu2e {

  TracerProduct::TracerProduct():
    val_(-1),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct default c'tor: " << *this;
  }

  TracerProduct::TracerProduct( int aval):
    val_(aval),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct  c'tor from int: " << *this;
  }

  TracerProduct::TracerProduct( TracerProduct const& rhs):
    val_(rhs.val_),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct copy c'tor: to: " 
                               << *this  << " from: " << rhs;
  }

  TracerProduct::~TracerProduct(){
    mf::LogVerbatim("Tracing") << "TracerProduct d'tor: " << *this;
  }

  TracerProduct& TracerProduct::operator=( TracerProduct const& rhs ){
    mf::LogVerbatim("Tracing") << "TracerProduct: operator= : to: "
                               << *this << " from: " << rhs;
    val_ = rhs.val_;
    return *this;
  }

  bool TracerProduct::operator<( TracerProduct const&  rhs ){
    mf::LogVerbatim("Tracing")  << "TracerProduct: operator< : lhs"
                                << *this << " rhs: " << rhs;
    bool answer(val_ < rhs.val_);
    return answer;
  } 

  void TracerProduct::print( ostream& ost ) const{
    ost << " ( serial number="
        << serial_ << " value="
        << val_ << " ) ";
  }

  // Maintain a monotonically increasing serial number, starting at 0.
  int TracerProduct::count(){
    static int c(-1);
    return ++c;
  }

}  // namespace mu2e
