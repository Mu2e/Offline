#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <limits>

namespace mu2e {
  // interface for updating straw hits
  struct KKStrawHitUpdater {
    KKStrawHitUpdater() : maxdoca_(-1.0), minprob_(-1.0), minddoca_(0), maxddoca_(-1) {}
    KKStrawHitUpdater(double maxdoca, double minprob, double minddoca, double maxddoca) :
      maxdoca_(maxdoca), minprob_(minprob),
      minddoca_(minddoca), maxddoca_(maxddoca) {}
    double maxdoca_; // maximum DOCA to still use a hit
    double minprob_; // minimum chisqquared probability to use a hit
    double minddoca_; // minimum DOCA to use drift information
    double maxddoca_; // maximum DOCA to use drift information
  };
}
#endif
