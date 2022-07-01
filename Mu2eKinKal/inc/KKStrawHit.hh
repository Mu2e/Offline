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
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
// Mu2e-specific classes
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/StrawProperties.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/DOCAStrawHitUpdater.hh"
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
  class CombinatoricStrawHitUpdater;

  template <class KTRAJ> class KKStrawHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using HIT = KinKal::Hit<KTRAJ>;
      using RESIDHIT = KinKal::ResidualHit<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions; should be defined outside this class FIXME
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&, StrawProperties const& sprops,
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
      double strawRadius() const { return rstraw_; }
      auto updaterAlgorithm() const { return algo_; }
      // Functions used in updating
      void setResiduals(WireHitState const& whstate, RESIDCOL& resids) const; // compute residuals WRT current reference given the state
      void updateResiduals() { setResiduals(whstate_, resids_); }
      bool insideStraw(CA const& ca) const; // decide if a CA is inside the straw
      CA unbiasedClosestApproach() const;
      auto updater() const { return algo_; }
      DriftInfo distanceToTime() const;
      void setState(WireHitState const& whstate,StrawHitUpdaters::algorithm algo);

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
      double rstraw_; // straw radius; this isn't a property of the straw
      ComboHit const& chit_; // reference to hit
      StrawHitIndex shindex_; // index to the StrawHit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
      StrawHitUpdaters::algorithm algo_; // record which updater was last used on this hit
      double mindoca_; // temporary
  };

  // struct to sort hits by time
  template <class KTRAJ> struct StrawHitTimeSort {
    using KKSTRAWHIT = KKStrawHit<KTRAJ>;
    using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
    bool operator ()( const KKSTRAWHITPTR& hit1, const KKSTRAWHITPTR& hit2) {
      return hit1->time() < hit2->time(); }
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate,
      StrawProperties const& sprops,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    bfield_(bfield), whstate_(whstate), wire_(ptca.sensorTraj()),
    ptca_(ptca.localTraj(),wire_,ptca.precision(),ptca.tpData(),ptca.dDdP(),ptca.dTdP()),
    rstraw_(sprops.strawInnerRadius()),
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

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // look for updaters
    unsigned nupdaters(0);
    auto cshu = miconfig.findUpdater<CombinatoricStrawHitUpdater>(); if(cshu != 0)nupdaters++;
    auto dshu = miconfig.findUpdater<DOCAStrawHitUpdater>(); if(dshu != 0)nupdaters++;
    auto nshu = miconfig.findUpdater<NullStrawHitUpdater>(); if(nshu != 0)nupdaters++;
    if(nupdaters != 1)throw cet::exception("RECO")<<"mu2e::KKStrawHit: StrawHit updater count error" << std::endl;
    if(nshu != 0 || dshu != 0){ // single-hit updaters: do the work here
      CA uca = unbiasedClosestApproach();
      if(uca.usable() && insideStraw (uca)){
        if(nshu != 0){
          if(first){
            algo_ = nshu->algorithm();
            whstate_ = nshu->wireHitState(uca.tpData());
            mindoca_ = strawRadius();
          }
        } else if(dshu != 0){
          if(first){
            algo_ = dshu->algorithm();
            whstate_ = dshu->wireHitState(uca.tpData());
            mindoca_ = dshu->minDOCA();
          }

        }
      } else {
        whstate_ = WireHitState::forcedinactive;
      }
    }
    // update residuals
    updateResiduals();
    // update the weight using the new residuals
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setState(WireHitState const& whstate,StrawHitUpdaters::algorithm algo) {
    whstate_ = whstate;
    algo_ = algo;
    updateResiduals();
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::setResiduals(WireHitState const& whstate, RESIDCOL& resids) const {
    if(whstate.active()){
      // compute drift parameters from PTCA
      DriftInfo dinfo = distanceToTime();
      if(whstate.useDrift()){
        // translate PCA to residual. Use ambiguity to convert drift time to a time difference.
        double dsign = whstate.lrSign()*ptca_.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
        double dt = ptca_.deltaT()-dinfo.tdrift_*dsign;
        // time differnce affects the residual both through the drift distance (DOCA) and the particle arrival time at the wire (TOCA)
        DVEC dRdP = ptca_.dDdP()*dsign/dinfo.vdrift_ - ptca_.dTdP();
        resids[tresid] = Residual(dt,dinfo.tdriftvar_,0.0,true,dRdP);
        resids[dresid] = Residual(); // distance residual isn't used when drift is
      } else {
        // interpret DOCA and TOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
        DVEC dRdP = -ptca_.lSign()*ptca_.dDdP();
        double dd,ddvar, dt, dtvar;
        if(algo_ == StrawHitUpdaters::null){
          dt = ptca_.deltaT() - chit_.driftTime() + 2.2; // correct using TOT drift time
          dtvar = 40.0;
          dd = ptca_.doca();
          ddvar = 2.7; // calibrated
        } else {
          dt = -0.5*mindoca_/dinfo.vdrift_; // adjust for the average drift time given the cut.
          dtvar = dinfo.tdriftvar_ + (mindoca_*mindoca_)/dinfo.vdrift_*dinfo.vdrift_*12.0;
          dd = ptca_.doca();
          ddvar =  (mindoca_*mindoca_)/3.0; //
        }
        resids[dresid] = Residual(dd,ddvar,0.0,true,dRdP);
        resids[tresid] = Residual(dt,dtvar,0.0,true,-ptca_.dTdP());
        // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
      }
    } else {
      resids[tresid] = resids[dresid] = Residual();
    }
  }

  template <class KTRAJ> Residual const& KKStrawHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >dresid)throw cet::exception("RECO")<<"mu2e::KKStrawHit: Invalid residual" << std::endl;
    return resids_[ires];
  }

  // compute estimated drift time using calibrated quantities (StrawResponse)
  template <class KTRAJ> DriftInfo KKStrawHit<KTRAJ>::distanceToTime() const {
    DriftInfo dinfo;
//    VEC3 bvec = bfield_.fieldVect(ptca_.particlePoca().Vect());
//    auto pdir = bvec.Cross(wire_.direction()).Unit(); // direction perp to wire and BFieldMap
//    VEC3 dvec = ptca_.delta().Vect();
//    double phi = asin(double(dvec.Unit().Dot(pdir))); // azimuth around the wire WRT the BField
//    POL2 drift(fabs(ptca_.doca()), phi);
    // for now, use simplified response model.
    dinfo.vdrift_ = sresponse_.driftConstantSpeed();
//    dinfo.tdrift_ = drift.R()/dinfo.vdrift_;
    dinfo.tdrift_ = fabs(ptca_.doca())/dinfo.vdrift_;
    dinfo.tdriftvar_ = 20.0; // temporary hack FIXME
    return dinfo;
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
