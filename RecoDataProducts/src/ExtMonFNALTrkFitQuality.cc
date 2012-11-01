#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"

namespace mu2e {
  std::ostream& operator<<(std::ostream&os, const ExtMonFNALTrkFitQuality& q) {
    return os<<"(chi2="<<q.chi2()<<",ndf="<<q.ndf()<<")";
  }
}
