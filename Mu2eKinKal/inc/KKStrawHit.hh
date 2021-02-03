#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Detector/WireHit.hh"
#include <memory>
namespace mu2e {
// struct for updating straw hits.  This is just for testing, use PanelHit for best results
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
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<PKTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, ComboHit const& chit, Tracker const& tracker, PTCA const& ptca);
// WireHit and Hit interface implementations
      void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      // specific to KKStrawHit: this has a constant drift speed
      virtual ~KKStrawHit(){}
      double timeVariance() const { return tvar_; }
    private:
      double dvel_; // current local drift velocity: this changes during update
      double tvar_; // current variance estimate of time measruement: will change during update
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, ComboHit const& chit, Tracker const& tracker, PTCA const& ptca) :
    WIREHIT(bfield,  tracker.straw(chit.strawId())
    
    
    wire,dxing,whstate), dvel_(driftspeed), tvar_(tvar) {}


  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // set precision
    WIREHIT::setPrecision(miconfig.tprec_);
    // update to move to the new trajectory
    this->update(pktraj);
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
      if(fabs(doca) > whupdater->mindoca_){
	newstate.lrambig_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
	newstate.dimension_ = WireHitState::time;
      } else if( fabs(doca) > whupdater->maxdoca_){
	newstate.dimension_ = WireHitState::none; // disable the hit
      } else {
	newstate.lrambig_ = WireHitState::null;
	if(whupdater->nulltime_)
	  newstate.dimension_ = WireHitState::both;
	else
	  newstate.dimension_ = WireHitState::distance;
      }
      WIREHIT::setHitState(newstate);
      // now update again in case the hit changed
      this->update(pktraj);
    }
    // OK if no updater is found, hits may be frozen this meta-iteration
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // simply translate distance to time using the fixed velocity
    dinfo.tdrift_ = drift.R()/dvel_;
    dinfo.vdrift_ = dvel_;
    dinfo.tdriftvar_ = tvar_;
  }

}
#endif
