//
// Choose a suitable range for histogramming a random variate distributed
// according to a mu2e::PoissonOrFixed distribution.
//
// $Id: PoissonHistogramBinning.cc,v 1.1 2012/02/14 17:09:54 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/02/14 17:09:54 $
//
// Contact author Rob Kutschke
//
// The mu2e::PoissonOrFixed distribution is constructed:
//
//   PoissonOrFixed::PoissonOrFixed( double mean );
//
// If (mean >0 ) then this distribution behaves as a Poisson distribution
// with the given mean.  If (mean <= 0 ) then the distribution returns
// a fixed value for all calls:
//
//   int n0 = static_cast<int>(floor(std::abs(mean)));
//
//  This class chooses binning as follows:
//
//  1) For mean<=0, return 3 bins from n0-1 to n0+1
//
//  2) For small positive values of mean, return one of several
//     fixed ranges.
//
//  3) For large positive values, (> 60 ) the range is
//     mean +/- width*sqrt(mean), rounded to integers with, at most, 100
//     bins.
//
#include "Mu2eUtilities/inc/PoissonHistogramBinning.hh"

#include <iostream>
#include <cmath>

using namespace std;

namespace mu2e {

  PoissonHistogramBinning::PoissonHistogramBinning( double mean, double width ):
    mean_(mean),
    width_(width),
    nbins_(100),
    xlow_(0),
    xhigh_(100.)
  {

    if ( mean <= 0 ){
      int n0 = static_cast<int>(floor(std::abs(mean)));
      nbins_ = 3;
      xlow_  = n0 - 1;
      xhigh_ = n0 + 1;

    } else if ( mean < 2 ){
      nbins_ = 10;
      xlow_  = 0.;
      xhigh_ = 10.;

    } else if ( mean < 25 ){
      nbins_ = 50;
      xlow_  = 0.;
      xhigh_ = 50.;

    } else if ( mean < 60 ){
      nbins_ = 100;
      xlow_  = 0.;
      xhigh_ = 100.;

    } else {
      int xl = static_cast<int>(floor(mean-width_*sqrt(mean)));
      int xh = static_cast<int>(ceil(mean+width_*sqrt(mean)));
      nbins_ = xh-xl;
      xlow_  = xl;
      xhigh_ = xh;
      if ( nbins_ > 100 ) {
        nbins_ = 100;
      }
    }

  } // end PoissonHistogramBinning c'tor

} // end namespace mu2e
