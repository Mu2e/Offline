#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
// mu2eKinKal classes
//KinKal classes
#include "KinKal/Detector/ResidualHit.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
// Mu2e-specific classes
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eKinKal/inc/CADSHU.hh"
#include "Offline/Mu2eKinKal/inc/DriftANNSHU.hh"
#include "Offline/Mu2eKinKal/inc/BkgANNSHU.hh"
#include "Offline/Mu2eKinKal/inc/Chi2SHU.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
namespace mu2e {
  using KinKal::BFieldMap;
  using KinKal::Line;
  using KinKal::MetaIterConfig;
  using KinKal::POL2;
  using KinKal::Residual;
  using KinKal::VEC3;
  using KinKal::DVEC;
  using KinKal::CAHint;
  using RESIDCOL = std::array<Residual,2>; // should be a struct FIXME

  template <class KTRAJ> class KKStrawHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using HIT = KinKal::Hit<KTRAJ>;
      using RESIDHIT = KinKal::ResidualHit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PCA const& pca,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // Hit interface implementations
      void updateState(MetaIterConfig const& config,bool first) override;
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      Residual const& refResidual(unsigned ires=Mu2eKinKal::tresid) const override;
      auto const& refResiduals() const { return resids_; }
      auto const& timeResidual() const { return resids_[Mu2eKinKal::tresid];}
      auto const& distResidual() const { return resids_[Mu2eKinKal::dresid];}
      double time() const override { return ca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return ca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      auto const& closestApproach() const { return ca_; }
      auto const& hitState() const { return whstate_; }
      auto const& wire() const { return wire_; }
      auto const& bfield() const { return bfield_; }
      auto precision() const { return ca_.precision(); }
      auto const& hit() const { return chit_; }
      auto const& straw() const { return straw_; }
      auto const& strawId() const { return straw_.id(); }
      auto const& strawHitIndex() const { return shindex_; }
      auto const& strawResponse() const { return sresponse_; }
      // Functions used in updating
      void setResiduals(MetaIterConfig const& miconfig, WireHitState const& whstate, RESIDCOL& resids) const; // compute residuals WRT current reference given the state
      CA unbiasedClosestApproach() const;
      auto updater() const { return whstate_.algo_; }
      void setState(WireHitState const& whstate); // allow cluster updaters to set the state directly
      DriftInfo fillDriftInfo() const;
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      WireHitState whstate_; // current state
      double dVar_; // drift distance variance value
      double dDdT_; // drift distance time derivative, crudely the drift velocity
      Line wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
      // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
      // (when multiplied by the propagation velocity).
      CA ca_; // reference time and position of closest approach to the wire; this is generally biased by the hit
      RESIDCOL resids_; // residuals WRT most recent reference
      ComboHit const& chit_; // reference to hit
      StrawHitIndex shindex_; // index to the StrawHit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
      // utility functions
      void updateWHS(MetaIterConfig const& miconfig);
  };

  // struct to sort hits by time
  template <class KTRAJ> struct StrawHitTimeSort {
    using KKSTRAWHIT = KKStrawHit<KTRAJ>;
    using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
    bool operator ()( const KKSTRAWHITPTR& hit1, const KKSTRAWHITPTR& hit2) {
      return hit1->time() < hit2->time(); }
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PCA const& pca,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    bfield_(bfield), whstate_(WireHitState::null), wire_(pca.sensorTraj()),
    ca_(pca.localTraj(),wire_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()),
    chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse)
  {
    // make sure this is a single-straw based ComboHit
    if(chit_.mask().level() != StrawIdMask::uniquestraw)
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: ComboHit doesn't correspond to a unique straw"<< std::endl;
  }

  template <class KTRAJ> KinKal::ClosestApproach<KTRAJ,Line> KKStrawHit<KTRAJ>::unbiasedClosestApproach() const {
    // compute the unbiased closest approach; this is brute force, but works
    auto const& ca = this->closestApproach();
    auto uparams = HIT::unbiasedParameters();
    KTRAJ utraj(uparams,ca.particleTraj());
    return CA(utraj,wire_,ca.hint(),ca.precision());
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // if we already computed PCA in the previous iteration, use that to set the hint.  This speeds convergence
    // otherwise use the time at the center of the wire, corrected for drift
    CAHint tphint = ca_.usable() ?  ca_.hint() : CAHint(wire_.timeAtMidpoint()-chit_.driftTime(),wire_.timeAtMidpoint());
    ca_ = CA(ktrajptr,wire_,tphint,precision());
    // check that we're on the right branch: we can move off if t0 changes a lot between iterations
    double dz = straw().origin().z() - ca_.particlePoca().Z();
    if((!ca_.usable()) || fabs(dz) >  100) { // need a better absolute scale; should come from KTRAJ FIXME
      tphint = CAHint(Mu2eKinKal::zTime(*ktrajptr,straw().origin().z(),wire_.timeAtMidpoint()), wire_.timeAtMidpoint());
      ca_ = CA(ktrajptr,wire_,tphint,precision());
      dz = straw().origin().z() - ca_.particlePoca().Z();
      if((!ca_.usable()) || fabs(dz) >  100) whstate_.state_ = WireHitState::unusable;// give up on 2nd try
    }
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateWHS(MetaIterConfig const& miconfig) {
    // search for updaters that work directly on StrawHits (not StrawHitClusters)
    auto cashu = miconfig.findUpdater<CADSHU>();
    auto driftshu = miconfig.findUpdater<DriftANNSHU>();
    auto bkgshu = miconfig.findUpdater<BkgANNSHU>();
    CA ca = unbiasedClosestApproach();
    if(ca.usable()){
      auto dinfo = fillDriftInfo();
      // there can be multiple updaters: apply them all
      if(bkgshu)whstate_ = bkgshu->wireHitState(whstate_,ca.tpData(),dinfo,chit_);
      if(cashu)whstate_ = cashu->wireHitState(whstate_,ca.tpData(),dinfo);
      if(driftshu)whstate_ = driftshu->wireHitState(whstate_,ca.tpData(),dinfo,chit_);
      if(whstate_.driftConstraint()){
        dVar_ = dinfo.driftHitVar();
        if(whstate_.constrainDriftDt()){
          dDdT_ = dinfo.driftVelocity_;
        } else{
          dDdT_ = 0.0;
        }
      } else {
        if(whstate_.nullDriftVar()) {
          dVar_ = dinfo.nullHitVar();
        } else {
          dVar_ = DriftInfo::maxdvar_;
        }
      }
    } else {
      whstate_.algo_ = StrawHitUpdaters::unknown;
      whstate_.state_ = WireHitState::unusable;
    }
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // first iteration of a new meta-iteration, update the wire hit state
    if(first)updateWHS(miconfig);
    // update residuals and weights every iteration, regardless of updater algorithm
    setResiduals(miconfig, whstate_, resids_);
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setState(WireHitState const& whstate) {
    whstate_ = whstate;
  }

  template <class KTRAJ> DriftInfo KKStrawHit<KTRAJ>::fillDriftInfo() const {
    double lorentzAngle = Mu2eKinKal::LorentzAngle(ca_.tpData(),ca_.particleTraj().bnom().Unit());
    return sresponse_.driftInfo(strawId(),ca_.deltaT(),lorentzAngle);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setResiduals(MetaIterConfig const& miconfig, WireHitState const& whstate, RESIDCOL& resids) const {
    // reset the residuals, using the fixed state from the last update
    resids[Mu2eKinKal::tresid] = resids[Mu2eKinKal::dresid] = Residual();
    if(whstate.active()){
      auto dinfo = fillDriftInfo();
      // optionally constrain DeltaT using the ComboHit TOT drift time or the absolute drift time
      if(whstate.constrainTOT()){
        double tvar = chit_.timeVar();
        double dt = ca_.deltaT() - chit_.driftTime();
        resids[Mu2eKinKal::tresid] = Residual(dt,tvar,0.0,true,ca_.dTdP());
      }
      // distance residual
      if(whstate.driftConstraint()){
        double dr = whstate.lrSign()*dinfo.rDrift_ - ca_.doca();
        DVEC dRdP = whstate.lrSign()*dDdT_*ca_.dTdP() -ca_.dDdP();
        resids[Mu2eKinKal::dresid] = Residual(dr,dVar_,0.0,true,dRdP);
      } else {
        // Null LR ambiguity. interpret DOCA against the wire directly as the spatial residual
        resids[Mu2eKinKal::dresid] = Residual(ca_.doca(),dVar_,0.0,true,ca_.dDdP());
        // optionally use the null hit time measurement to constrain t0
        if(whstate.constrainAbsDriftDt()){
          double dt = ca_.deltaT() - sresponse_.strawDrift().D2T(fabs(ca_.doca()),dinfo.LorentzAngle_);
          double tvar = dinfo.driftTimeVar();
          // this overwrites the TOT time constraint, in principle both can be used TODO
          resids[Mu2eKinKal::tresid] = Residual(dt,tvar,0.0,true,ca_.dTdP());
        }
      }
    }
  }

  template <class KTRAJ> Residual const& KKStrawHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >Mu2eKinKal::tresid)throw cet::exception("RECO")<<"mu2e::KKStrawHit: Invalid residual" << std::endl;
    return resids_[ires];
  }

  template<class KTRAJ> void KKStrawHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " KKStrawHit time " << this->time();
    switch(whstate_.state_) {
      case WireHitState::inactive:
        ost << "inactive";
        break;
      case WireHitState::left:
        ost << "left";
        break;
      case WireHitState::right:
        ost << "right";
        break;
      case WireHitState::null: default:
        ost << "null";
        break;
    }
    if(detail > 0){
      if(resids_[Mu2eKinKal::tresid].active())
        ost << " Active Time " << resids_[Mu2eKinKal::tresid];
      if(resids_[Mu2eKinKal::dresid].active())
        ost << " Active Dist " << resids_[Mu2eKinKal::dresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << " Ref " << ca_.tpData() << std::endl;
    }
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
