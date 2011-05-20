#ifndef GeneralUtilities_Binning_hh
#define GeneralUtilities_Binning_hh
//
// A definition of the binning of a histogram.
//
// $Id: Binning.hh,v 1.4 2011/05/20 20:18:23 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:18:23 $
//
// Original author Rob Kutschke
//

#include <cmath>
#include <limits>

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

#endif /* GeneralUtilities_Binning_hh */
