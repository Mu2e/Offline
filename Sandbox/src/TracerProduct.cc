//
// A test class that makes printout whenever its methods are called.
//
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

  TracerProduct::~TracerProduct(){
    mf::LogVerbatim("Tracing") << "TracerProduct d'tor: " << *this;
  }

#ifndef __GCCXML__

  TracerProduct::TracerProduct( int aval):
    val_(aval),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct  c'tor from int: " << *this;
  }

  TracerProduct::TracerProduct( TracerProduct const& rhs):
    val_(-1),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct copy c'tor: to: "
                               << *this  << " from: " << rhs;
    val_=rhs.val_;
  }

  TracerProduct::TracerProduct( TracerProduct && rhs):
    val_(-1),
    serial_(count()){
    mf::LogVerbatim("Tracing") << "TracerProduct move c'tor: to: "
                               << *this  << " from: " << rhs;
    val_=rhs.val_;
    rhs.val_=-1;
  }


  TracerProduct& TracerProduct::operator=( TracerProduct const& rhs ){
    mf::LogVerbatim("Tracing") << "TracerProduct: operator= : to: "
                               << *this << " from: " << rhs;
    val_ = rhs.val_;
    return *this;
  }

  TracerProduct& TracerProduct::operator=( TracerProduct && rhs ){
    mf::LogVerbatim("Tracing") << "TracerProduct: move operator= : to: "
                               << *this << " from: " << rhs;
    val_ = rhs.val_;
    rhs.val_ = -1;
    return *this;
  }

  void TracerProduct::swap( TracerProduct& rhs ){
    mf::LogVerbatim("Tracing") << "TracerProduct: swap : lhs: "
                               << *this << " rhs: " << rhs;
    std::swap(this->serial_, rhs.serial_ );
    std::swap(this->val_, rhs.val_ );
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

#endif

}  // namespace mu2e
