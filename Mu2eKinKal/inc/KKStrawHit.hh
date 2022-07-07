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
#include "Offline/Mu2eKinKal/inc/DOCAStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
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
  //  class CombinatoricStrawHitUpdater;

  template <class KTRAJ> class KKStrawHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using HIT = KinKal::Hit<KTRAJ>;
      using RESIDHIT = KinKal::ResidualHit<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions; should be defined outside this class FIXME
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // Hit interface implementations
      void updateState(MetaIterConfig const& config,bool first) override;
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      Residual const& refResidual(unsigned ires=tresid) const override;
      auto const& refResiduals() const { return resids_; }
      auto const& timeResidual() const { return resids_[tresid];}
      auto const& distResidual() const { return resids_[dresid];}
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
      auto updaterAlgorithm() const { return algo_; }
      // Functions used in updating
      void setResiduals(MetaIterConfig const& miconfig, WireHitState const& whstate, RESIDCOL& resids) const; // compute residuals WRT current reference given the state
      void updateResiduals(MetaIterConfig const& miconfig);
      bool insideStraw(CA const& ca) const; // decide if a CA is inside the straw
      CA unbiasedClosestApproach() const;
      auto updater() const { return algo_; }
      void setState(MetaIterConfig const& miconfig,WireHitState const& whstate);
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
      StrawHitUpdaters::algorithm algo_; // most recent algorithm
      // utility functions
      void updateWHS(MetaIterConfig const& miconfig);
      double minDOCA(MetaIterConfig const& miconfig) const;
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
    chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse), algo_(StrawHitUpdaters::none)
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
    if(!ptca_.usable()) throw cet::exception("RECO")<<"mu2e::KKStrawHit: WireHit TPOCA failure" << std::endl;
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateWHS(MetaIterConfig const& miconfig) {
    // look for updaters; there must be exactly 1
    unsigned nupdaters(0);
    auto cshu = miconfig.findUpdater<CombinatoricStrawHitUpdater>();
    auto dshu = miconfig.findUpdater<DOCAStrawHitUpdater>();
    auto nshu = miconfig.findUpdater<NullStrawHitUpdater>();
    if(cshu){
      algo_ = StrawHitUpdaters::Combinatoric;
      nupdaters++;
    }
    if(dshu){
      algo_ = dshu->algorithm();
      nupdaters++;
    }
    if(nshu){
      algo_ = nshu->algorithm();
      nupdaters++;
    }
    if(nupdaters != 1)throw cet::exception("RECO")<<"mu2e::KKStrawHit: StrawHit updater count error" << std::endl;
    // Combo updater sets WireHitState in StrawHitCluster; leave it unchanged here
    //  For locally-operating updaters, actually update the state
    if(StrawHitUpdaters::updateStrawHits(algo_)){
      CA uca = unbiasedClosestApproach();
      if(uca.usable() && insideStraw (uca)){
        if(dshu) whstate_ = dshu->wireHitState(uca.tpData());
        if(nshu) whstate_ = nshu->wireHitState(uca.tpData());
      } else {
        whstate_ = WireHitState::forcedinactive;
      }
    }
  }

  template <class KTRAJ> double KKStrawHit<KTRAJ>::minDOCA(MetaIterConfig const& miconfig) const {
    double dmin(0.0);
    if(algo_ == StrawHitUpdaters::DOCA){
      auto dshu = miconfig.findUpdater<DOCAStrawHitUpdater>();
      if(!dshu)throw cet::exception("RECO")<<"mu2e::KKStrawHit: missing updater" << std::endl;
      dmin = dshu->minDOCA();
    } else if(algo_ == StrawHitUpdaters::Combinatoric){
      auto cshu = miconfig.findUpdater<CombinatoricStrawHitUpdater>();
      if(!cshu)throw cet::exception("RECO")<<"mu2e::KKStrawHit: missing updater" << std::endl;
      dmin = cshu->minDOCA();
    } else
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: missing updater" << std::endl;
    return dmin;
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // first iteration of a new meta-iteratin, update the wire hit state
    if(first)updateWHS(miconfig);
    // update residuals and weights if the algorithm operates directly on StrawHits
    if(StrawHitUpdaters::updateStrawHits(algo_))updateResiduals(miconfig);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setState(MetaIterConfig const& miconfig, WireHitState const& whstate) {
    whstate_ = whstate;
    updateResiduals(miconfig);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateResiduals(MetaIterConfig const& miconfig) {
    setResiduals(miconfig, whstate_, resids_);
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setResiduals(MetaIterConfig const& miconfig, WireHitState const& whstate, RESIDCOL& resids) const {
    // reset the residuals
    resids[tresid] = resids[dresid] = Residual();
    if(whstate.active()){
      if(whstate.useDrift()){
        // Transldate DOCA to a drift time. ignore phi (Lorentz effects) for now: TODO
        auto dinfo = sresponse_.driftInfoAtDistance(strawId(), fabs(ptca_.doca()), 0.0);
        // residual itself MUST be computed WRT the reference parameters
        double dsign = whstate.lrSign()*ptca_.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
        double dt = ptca_.deltaT()-dinfo.time*dsign;
        // time differnce affects the residual both through the drift distance (DOCA) and the particle arrival time at the wire (TOCA)
        // temporary fix
        DVEC dRdP = ptca_.dDdP()*dsign*dinfo.invSpeed - ptca_.dTdP();
        resids[tresid] = Residual(dt,dinfo.variance,0.0,true,dRdP);
        // distance residual isn't used when drift is
      } else {
        // Null state. interpret DOCA and TOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
        DVEC dRdP = -ptca_.lSign()*ptca_.dDdP();
        // distance residual is always WRT the wire
        double dd = ptca_.doca();
        // variances and time residuals depend on the algorithm
        double ddvar, dt, dtvar;
        if(algo_ == StrawHitUpdaters::null){
          auto nshu = miconfig.findUpdater<NullStrawHitUpdater>();
          if(!nshu)throw cet::exception("RECO")<<"mu2e::KKStrawHit: missing updater" << std::endl;
          // use the combo-hit time to set the time residual; the correction should be calibrated out TODO
          dt = ptca_.deltaT() - chit_.driftTime() - 0.85;
          dtvar = 50.0; // should come from ComboHit TODO
          ddvar = nshu->distVariance(); // this should come from a prodition TODO
        } else {
          // other null updaters are based on an 'effective' DOCA
          double dmin = minDOCA(miconfig);
          // get drift properties at this effective DOCA
          auto dinfo = sresponse_.driftInfoAtDistance(strawId(),0.5*dmin,0.0);
          // distance variance is geometric, based on the effective distance
          static double invthree(1.0/3.0);
          ddvar = invthree*dmin*dmin;
          // time residual is delta-T corrected for the average drift time
          dt = ptca_.deltaT() -dinfo.time;
          // time variance has 2 parts: intrinsic and residual drift
          dtvar = dinfo.variance + 0.25*ddvar*dinfo.invSpeed*dinfo.invSpeed;
        }
        resids[dresid] = Residual(dd,ddvar,0.0,true,dRdP);
        resids[tresid] = Residual(dt,dtvar,0.0,true,-ptca_.dTdP());
        // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
      }
    }
  }

  template <class KTRAJ> Residual const& KKStrawHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >dresid)throw cet::exception("RECO")<<"mu2e::KKStrawHit: Invalid residual" << std::endl;
    return resids_[ires];
  }

  template <class KTRAJ> bool KKStrawHit<KTRAJ>::insideStraw(CA const& ca) const {
    static const double ubuffer(10.0); // should be a parameter FIXME
    // compute the position along the wire and compare to the 1/2 length
    // have to translate from CLHEP, should be native to Straw FIXME
    double upos = VEC3(straw_.wireDirection()).Dot((ca.sensorPoca().Vect() - VEC3(straw_.origin())));
    return fabs(upos) < straw_.halfLength() + ubuffer;
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
      if(resids_[tresid].active())
        ost << " Active Time " << resids_[tresid];
      if(resids_[dresid].active())
        ost << " Active Dist " << resids_[dresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << ptca_.tpData() << std::endl;
    }
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
