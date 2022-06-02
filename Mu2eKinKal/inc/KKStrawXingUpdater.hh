#ifndef Mu2eKinKal_KKStrawXingUpdater_hh
#define Mu2eKinKal_KKStrawXingUpdater_hh
#include "KinKal/Detector/WireHitStructs.hh"
namespace mu2e {
  // simple struct to hold crossing calculation configuration parameters
  using KinKal::WireHitState;
  struct KKStrawXingUpdater {
    WireHitState hitstate_; // associated hit state (inactive if no hit)
    double maxdocasig_; // maximum doca error to use non-averaged value
    double maxdoca_; // maximum DOCA to include this straw's material
    double maxddoca_; // maximum DOCA to use 'exact' calculation, otherwise average over all physical impact parameters
    // default constructor is functional but will always use the impact-parameter averaged material
    KKStrawXingUpdater() : hitstate_(WireHitState::null), maxdocasig_(-1.0), maxdoca_(0.0), maxddoca_(0.0) {}
    KKStrawXingUpdater(double maxdocasig, double maxdoca, double maxddoca) : hitstate_(WireHitState::null),
    maxdocasig_(maxdocasig), maxdoca_(maxdoca), maxddoca_(maxddoca){}
  };
}
#endif
