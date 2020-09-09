#ifndef GeneralUtilities_Binning_hh
#define GeneralUtilities_Binning_hh
//
// A definition of the binning of a histogram.
//
//
// Original author Rob Kutschke
//

#include <cmath>

class Binning{

public:
  typedef unsigned long IndexType;
  static const IndexType nobin; // like string::npos

  Binning();
  Binning(IndexType nbins, double low, double high);

  // Accessors.
  IndexType nbins() const {return nbins_;}
  double low()   const {return low_;  }
  double high()  const {return high_; }

  // returns nobin for out of range
  IndexType findBin(double x) const;

  double binCenter(IndexType ibin) const;

  double binWidth() const { return binwidth_; }

private:

  // Number of bins.
  IndexType nbins_;

  // Low edge of low bin.
  double low_;

  // High edge of high bin.
  double high_;

  double binwidth_;

};

#endif /* GeneralUtilities_Binning_hh */
