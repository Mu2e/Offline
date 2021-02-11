#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
//KinKal classes
#include "KinKal/Detector/WireHit.hh"
#include "KinKal/Detector/StrawXing.hh"
// Mu2e-specific classes
#include "TrackerGeom/inc/Straw.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
namespace mu2e {
// struct for updating straw hits.  This is just for testing, use PanelHit updating for best results
  struct KKStrawHitUpdater {
    double mindoca_; // minimum DOCA value to set an ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    bool nulltime_; // constrain time when hit has null ambiguity
    double rcell_; // straw radius
    KKStrawHitUpdater(double mindoca,double maxdoca, bool nulltime) : mindoca_(mindoca), maxdoca_(maxdoca), nulltime_(nulltime) {}
  };
  using KinKal::BFieldMap;
  using KinKal::WireHitState;
  using KinKal::Line;
  using KinKal::MetaIterConfig;
  using KinKal::DriftInfo;
  using KinKal::POL2;
  template <class KTRAJ> class KKStrawHit : public KinKal::WireHit<KTRAJ> {
    public:
      using WIREHIT = KinKal::WireHit<KTRAJ>;
      using STRAWXING = KinKal::StrawXing<KTRAJ>;
      using STRAWXINGPTR = std::shared_ptr<STRAWXING>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, STRAWXINGPTR const& strawxing, WireHitState const&, 
	  ComboHit const& chit, Straw const& straw, StrawResponse const& sresponse);
// WireHit and Hit interface implementations
      void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific to KKStrawHit: this has a constant drift speed
      virtual ~KKStrawHit(){}
      // accessors
      ComboHit const& hit() const { return chit_; }
      Straw const& straw() const { return straw_; }
    private:
      ComboHit const& chit_; // reference to hit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, STRAWXINGPTR const& strawxing, WireHitState const& whstate,
      ComboHit const& chit, Straw const& straw, StrawResponse const& sresponse) : 
    WIREHIT(bfield,ptca,strawxing,whstate), chit_(chit), straw_(straw), sresponse_(sresponse)
  {
  // make sure this is a single-straw based ComboHit
    if(chit_.mask().level() != StrawIdMask::uniquestraw)
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: ComboHit doesn't correspond to a unique straw"<< endl;
  }

// the purpose of this class is to allow updating using calibrated quantities
  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // set precision
    WIREHIT::setPrecision(miconfig.tprec_);
    // move to the new trajectory; this updates the closest approach
    this->update(pktraj);
    if(miconfig.updatehits_){
      // find the wire hit updater in the update params.  There should be 0 or 1
      const KKStrawHitUpdater* whupdater(0);
      for(auto const& uparams : miconfig.updaters_){
	auto const* whu = std::any_cast<KKStrawHitUpdater>(&uparams);
	if(whu != 0){
	  if(whupdater !=0) throw std::invalid_argument("Multiple KKStrawHitUpdaters found");
	  whupdater = whu;
	}
      }
      // crude updating of ambiguity and activity based on DOCA
      if(whupdater != 0){
	// start with existing state
	WireHitState newstate = WIREHIT::hitState();
	newstate.nullvar_ = whupdater->mindoca_*whupdater->mindoca_/3.0; // RMS of flat distribution beteween +- mindoca
	double doca = fabs(WIREHIT::closestApproach().doca());
	if( fabs(doca) > whupdater->maxdoca_){
	  newstate.dimension_ = WireHitState::none; // disable the hit
	} else if(fabs(doca) > whupdater->mindoca_){
	  newstate.lrambig_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
	  newstate.dimension_ = WireHitState::time;
	} else {
	  newstate.lrambig_ = WireHitState::null;
	  if(whupdater->nulltime_)
	    newstate.dimension_ = WireHitState::both;
	  else
	    newstate.dimension_ = WireHitState::distance;
	}
	this->setHitState(newstate);
	// now update again in case the caches changed
	this->update(pktraj);
      }
      // OK if no updater is found, hits may be frozen this meta-iteration
    }
  }

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
