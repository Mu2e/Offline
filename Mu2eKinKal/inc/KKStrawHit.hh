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
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
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
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // Hit interface implementations
      void updateState(MetaIterConfig const& config,bool first) override;
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      Residual const& refResidual(unsigned ires=Mu2eKinKal::tresid) const override;
      auto const& refResiduals() const { return resids_; }
      auto const& timeResidual() const { return resids_[Mu2eKinKal::tresid];}
      auto const& distResidual() const { return resids_[Mu2eKinKal::dresid];}
      double time() const override { return ptca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return ptca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      auto const& closestApproach() const { return ptca_; }
      auto const& hitState() const { return whstate_; }
      auto const& wire() const { return wire_; }
      auto const& bfield() const { return bfield_; }
      auto precision() const { return ptca_.precision(); }
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
      DriftInfo driftInfo() const;
   private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      WireHitState whstate_; // current state
      Line wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
      // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
      // (when multiplied by the propagation velocity).
      CA ptca_; // reference time and position of closest approach to the wire; this is generally biased by the hit
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

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    bfield_(bfield), whstate_(WireHitState::null), wire_(ptca.sensorTraj()),
    ptca_(ptca.localTraj(),wire_,ptca.precision(),ptca.tpData(),ptca.dDdP(),ptca.dTdP()),
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
    return CA(utraj,this->wire(),ca.hint(),ca.precision());
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // if we already computed PCA in the previous iteration, use that to set the hint.  This speeds convergence
    // otherwise use the time at the center of the wire
    CAHint tphint = ptca_.usable() ?  ptca_.hint() : CAHint(wire_.range().mid(),wire_.range().mid());
    ptca_ = CA(ktrajptr,wire_,tphint,precision());
    // check that we're on the right branch: we can move off if t0 changes a lot between iterations
    double dz = straw().origin().z() - ptca_.particlePoca().Z();
    if(fabs(dz) >  100) { // need a better absolute scale; should come from KTRAJ FIXME
      tphint = CAHint(Mu2eKinKal::zTime(*ktrajptr,straw().origin().z(),wire_.range().mid()), wire_.range().mid());
      ptca_ = CA(ktrajptr,wire_,tphint,precision());
      dz = straw().origin().z() - ptca_.particlePoca().Z();
      if(fabs(dz) >  100){
        whstate_.state_ = WireHitState::inactive;
        whstate_.usable_ = false; // failed to find the right branch: disable this hit for now, but we should do an explicit search TODO
      }
    }
    if(!ptca_.usable()){
      whstate_.state_ = WireHitState::inactive;
      whstate_.usable_ = false;
    }
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateWHS(MetaIterConfig const& miconfig) {
    if(!whstate_.frozen_){
      unsigned nupdaters(0);
      // only search for updaters that work directly on StrawHits (not StrawHitClusters)
      auto pshu = miconfig.findUpdater<PTCAStrawHitUpdater>();
      auto nshu = miconfig.findUpdater<NullStrawHitUpdater>();
      StrawHitUpdater const* shu(0);
      if(pshu){ shu = pshu; ++nupdaters; }
      if(nshu){ shu = nshu; ++nupdaters; }
      if(nupdaters > 1)throw cet::exception("RECO")<<"mu2e::KKStrawHit: multiple updaters" << std::endl;
      if(shu != 0){
        CA ca = shu->useUnbiasedClosestApproach() ? unbiasedClosestApproach() : ptca_;
        whstate_.usable_ = ca.usable();
        if(whstate_.usable_){
          auto dinfo = DriftInfo();
          whstate_ = shu->wireHitState(ca.tpData(),dinfo,chit_);
        }else
          whstate_.state_ = WireHitState::inactive;
      } // no updater means don't change the state
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

  template <class KTRAJ> DriftInfo KKStrawHit<KTRAJ>::driftInfo() const {
    DriftInfo dinfo;
    dinfo.LorentzAngle_ = Mu2eKinKal::LorentzAngle(ptca_.tpData(),ptca_.particleTraj().bnom().Unit());
    dinfo.driftDistance_ = sresponse_.driftTimeToDistance(strawId(),ptca_.deltaT(),dinfo.LorentzAngle_);
//    dinfo.driftDistanceError_ = sresponse_.driftDistanceError(strawId(),fabs(ptca_.doca()),dinfo.LorentzAngle_);
//    dinfo.driftVelocity_ = sresponse_.driftInstantSpeed(strawId(),fabs(ptca_.doca()),dinfo.LorentzAngle_,true);
    dinfo.driftDistanceError_ = sresponse_.driftDistanceError(strawId(),fabs(dinfo.driftDistance_),dinfo.LorentzAngle_);
    dinfo.driftVelocity_ = sresponse_.driftInstantSpeed(strawId(),fabs(dinfo.driftDistance_),dinfo.LorentzAngle_,true);
    return dinfo;
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setResiduals(MetaIterConfig const& miconfig, WireHitState const& whstate, RESIDCOL& resids) const {
    // reset the residuals
    resids[Mu2eKinKal::tresid] = resids[Mu2eKinKal::dresid] = Residual();
    auto const& nhinfo = whstate.nhinfo_;
    if(whstate.active()){
      if(whstate.useDrift()){
        DriftInfo dinfo = driftInfo();
        double rvar = dinfo.driftDistanceError_*dinfo.driftDistanceError_;
        double dr = whstate.lrSign()*dinfo.driftDistance_ - ptca_.doca();
        DVEC dRdP = ptca_.lSign()*ptca_.dDdP() - whstate.lrSign()*dinfo.driftVelocity_*ptca_.dTdP();
        resids[Mu2eKinKal::dresid] = Residual(dr,rvar,0.0,true,dRdP);
      } else {
        // Null state. interpret DOCA against the wire directly as a residual.  We have to take the DOCA sign out of the derivatives
        resids[Mu2eKinKal::dresid] = Residual(ptca_.doca(),nhinfo.dvar_,0.0,true,-ptca_.lSign()*ptca_.dDdP());
        // optionally also constrain the time
        if(nhinfo.useTime()){
          // time residual is deltaT
          double dt = ptca_.deltaT() - nhinfo.toff_;
          if(nhinfo.useComboDriftTime())dt -= chit_.driftTime();  // Null updater uses combo hit time
          resids[Mu2eKinKal::tresid] = Residual(dt,nhinfo.tvar_,0.0,true,-ptca_.dTdP());
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
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << ptca_.tpData() << std::endl;
    }
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
