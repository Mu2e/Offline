#ifndef GeneralUtilities_RMS_hh
#define GeneralUtilities_RMS_hh

//
// Compute mean, RMS and some other stats about a set of numbers.
//
// $Id: RMS.hh,v 1.1 2011/12/16 23:13:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/12/16 23:13:14 $
//
// Contact person Rob Kutschke
//
// Notes
// 1) Consider throwing in places in which the return value is forced to zero.
//    See comments in code. If we do this, we will have to provide saftey features
//    for the operator <<; maybe a "bool nothrow" argument on each accessor?
//

#include <iostream>

class RMS{

 public:

  RMS ( ): n_(0), sum_(0.), sumsq_(0.) { }
  RMS ( double x ): n_(1), sum_(x), sumsq_(x*x) { }

  // Accept compiler written d'tor, copy c'tor and assignment operator.

  // Number of entries.
  int n() const { return n_; }

  // Mean of the entries.
  double mean() const { return ( n_ == 0 ) ? 0. :sum_/n_; }

  // RMS about the mean; the argument is the number of parameters determined by the data.
  double rms  ( int npar = 1 ) const;

  // RMS about 0; the argument is the number of parameters determined by the data.
  double rms0 ( int npar = 0 ) const;

  // The underlying sums.
  double sum()   const { return sum_; }
  double sumSq() const { return sumsq_; }

  // Errors in the mean and the rms.
  double errorMean( int npar = 1 ) const;
  double errorRMS ( int npar = 1 ) const;

  // Accumulate one entry.
  void accumulate ( double x ){
    ++n_;
    sum_   += x;
    sumsq_ += x*x;
  }

  // Clear all entries.
  void clear(){
    n_     = 0.;
    sum_   = 0.;
    sumsq_ = 0.;
  }

private:
  int    n_;
  double sum_;
  double sumsq_;

};

inline std::ostream& operator<<(std::ostream& ost,
                               const RMS& rms ){
  ost << "( "
      << rms.n()         << ", "
      << rms.mean()      << " +/- "
      << rms.errorMean() << ", "
      << rms.rms()       << " +/- "
      << rms.errorRMS()
      << " )";
  return ost;
}

#endif /* GeneralUtilities_RMS_hh */
