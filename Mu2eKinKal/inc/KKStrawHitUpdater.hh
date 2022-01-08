#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <limits>

namespace mu2e {
  // interface for updating straw hits
  struct KKStrawHitUpdater {
    KKStrawHitUpdater() : mindoca_(std::numeric_limits<float>::max()), maxdoca_(-1.0), minprob_(-1.0) {}
    KKStrawHitUpdater(double mindoca, double maxdoca, double minprob) : mindoca_(mindoca), maxdoca_(maxdoca), minprob_(minprob){}
    double mindoca_; // minimum DOCA value to use drift information
    double maxdoca_; // maximum DOCA to still use a hit
    double minprob_; // minimum chisqquared probability to use a hit
  };
}
#endif
