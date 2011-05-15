#ifndef Binning_HH
#define Binning_HH
//
// A definition of the binning of a histogram.
//
// $Id: Binning.hh,v 1.1 2011/05/15 21:14:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/15 21:14:09 $
//
// Original author Rob Kutschke
//

#include <limits>
#include <cmath>

class Binning{
  
public:

  Binning():
    nbins_(1),
    low_(0.),
    high_(1.)
  {}

  Binning( int nbins, double low, double high):
    nbins_(nbins),
    low_(low),
    high_(high){
  }
  
  // Accessors.
  int    nbins() const {return nbins_;}
  double low()   const {return low_;  }
  double high()  const {return high_; }

private:

  // Number of bins.
  int nbins_;

  // Low edge of low bin.
  double low_;

  // High edge of high bin.
  double high_;

};

#endif
