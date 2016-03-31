#include "GeneralUtilities/inc/Binning.hh"
#include <cmath>

const Binning::IndexType Binning::nobin(-1);

Binning::Binning()
  : Binning{1,0.,1.}
{}

Binning::Binning(IndexType nbins, double low, double high)
  : nbins_(nbins)
  , low_(low)
  , high_(high)
  , binwidth_((high_-low_)/nbins_)
{}

Binning::IndexType Binning::findBin(double x) const {
  IndexType ibin = nobin;

  // we use the [) interval.  NB: comparison with a NaN gives false.
  if( (low_ <= x) && (x < high_) ) {

    ibin = (x - low_)/binwidth_;

    if(ibin >= nbins_) {
      ibin = nobin;
    }
  }

  return ibin;
}

double Binning::binCenter(IndexType i) const {
  if((i != nobin) && (i<nbins_)) {
    return low_ + binwidth_*(0.5+i);
  }
  else {
    return nan("");
  }
}
