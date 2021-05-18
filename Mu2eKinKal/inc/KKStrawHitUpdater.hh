#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <limits>
using KinKal::ClosestApproachData;
using KinKal::WireHitState;

namespace mu2e {
  // interface for updating straw hits
  class KKStrawHitUpdater {
    public:
      KKStrawHitUpdater() : mindoca_(std::numeric_limits<float>::max()), maxdoca_(-1.0), nulldim_(WireHitState::both) {}
      KKStrawHitUpdater(double mindoca, double maxdoca, double maxchi, WireHitState::Dimension nulldim) : mindoca_(mindoca), maxdoca_(maxdoca), maxchi_(maxchi), nulldim_(nulldim) {}
      template <class KKSTRAWHIT> void update(KKSTRAWHIT& hit) const;
    private:
      double mindoca_; // minimum DOCA value to use drift information
      double maxdoca_; // maximum DOCA to still use a hit
      double maxchi_; // maximum chi to use a hit
      WireHitState::Dimension nulldim_; // constrain dimension for null hits
  };
// update a particular hit
  template <class KKSTRAWHIT> void KKStrawHitUpdater::update(KKSTRAWHIT& hit) const {
    auto& hitstate = hit.hitState();
    auto& poca = hit.closestApproach();

    double absdoca = fabs(poca.doca()); 
    if( absdoca > maxdoca_){ // hit is too far from the wire: disable it
      hitstate.dimension_ = WireHitState::none; // disable the hit
      } else if(absdoca > mindoca_ && absdoca < maxdoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
	hitstate.lrambig_ = poca.doca() > 0.0 ? WireHitState::right : WireHitState::left;
	hitstate.dimension_ = WireHitState::time;
      } else { // hit very close to the wire: ambiguity information is unusable
	hitstate.lrambig_ = WireHitState::null;
	hitstate.dimension_ = nulldim_;
      }
    }
  }
#endif
