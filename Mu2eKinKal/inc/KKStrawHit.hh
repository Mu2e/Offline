#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
// mu2eKinKal classes
#include "Mu2eKinKal/inc/KKStrawHitUpdater.hh"
//KinKal classes
#include "KinKal/Detector/WireHit.hh"
// Mu2e-specific classes
#include "TrackerGeom/inc/Straw.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
namespace mu2e {
  using KinKal::BFieldMap;
  using KinKal::WireHitState;
  using KinKal::Line;
  using KinKal::MetaIterConfig;
  using KinKal::DriftInfo;
  using KinKal::POL2;

  template <class KTRAJ> class KKStrawHit : public KinKal::WireHit<KTRAJ> {
    public:
      using WIREHIT = KinKal::WireHit<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&, 
	  ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
// WireHit and Hit interface implementations
      void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific to KKStrawHit: this has a constant drift speed
      virtual ~KKStrawHit(){}
      // accessors
      ComboHit const& hit() const { return chit_; }
      Straw const& straw() const { return straw_; }
      StrawHitIndex const& strawHitIndex() const { return shindex_; }
    private:
      ComboHit const& chit_; // reference to hit
      StrawHitIndex shindex_; // index to the StrawHit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) : 
    WIREHIT(bfield,ptca,whstate), chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse)
  {
  // make sure this is a single-straw based ComboHit
    if(chit_.mask().level() != StrawIdMask::uniquestraw)
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: ComboHit doesn't correspond to a unique straw"<< endl;
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // set precision
    WIREHIT::setPrecision(miconfig.tprec_);
    // move to the new trajectory; this updates the closest approach
    this->update(pktraj);
    // look for an updater
    const KKStrawHitUpdater* whupdater(0);
    for(auto const& uparams : miconfig.updaters_){
      const KKStrawHitUpdater* whu = std::any_cast<KKStrawHitUpdater>(&uparams);
      if(whu != 0){
	if(whupdater !=0) throw std::invalid_argument("Multiple KKStrawHitUpdaters found");
	whupdater = whu;
      }
    }
    if(whupdater != 0){
      whupdater->updateState(this->hitState(), this->closestApproach());
      // now update again in case the caches changed
      this->update(pktraj);
    }
  }

  // the purpose of this class is to allow computing the drift using calibrated quantities
  template <class KTRAJ> void KKStrawHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    dinfo.tdrift_ = sresponse_.driftDistanceToTime(chit_.strawId(),drift.R(),drift.Phi());
    dinfo.vdrift_ = sresponse_.driftInstantSpeed(chit_.strawId(),drift.R(),drift.Phi());
    auto derr = sresponse_.driftDistanceError(chit_.strawId(),drift.R(),drift.Phi(), this->closestApproach().doca());
    dinfo.tdriftvar_ = std::pow(derr/dinfo.vdrift_,(int)2);
  }

  template<class KTRAJ> void KKStrawHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " KKStrawHit";
    WIREHIT::print(ost,detail);
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
