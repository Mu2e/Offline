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
#include "KinKal/Detector/WireHitStructs.hh"
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
#include "Offline/TrackerGeom/inc/StrawProperties.hh"
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
  using KinKal::Residual;
  using KinKal::VEC3;
  using KinKal::DVEC;
  using KinKal::CAHint;

  template <class KTRAJ> class KKStrawHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using HIT = KinKal::Hit<KTRAJ>;
      using RESIDHIT = KinKal::ResidualHit<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&, StrawProperties const& sprops,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // Hit interface implementations
      void updateState(MetaIterConfig const& config,bool first) override;
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      Residual const& refResidual(unsigned ires=tresid) const override;
      double time() const override { return ptca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return ptca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific to KKStrawHit:
      auto const& closestApproach() const { return ptca_; }
      CA unbiasedClosestApproach() const;
      auto const& hitState() const { return whstate_; }
      auto const& wire() const { return wire_; }
      auto const& bfield() const { return bfield_; }
      auto precision() const { return ptca_.precision(); }
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const;
      double nullVariance(Dimension dim,DriftInfo const& dinfo) const;
      double nullOffset(Dimension dim,DriftInfo const& dinfo) const;
     // accessors
      auto const& hit() const { return chit_; }
      auto const& straw() const { return straw_; }
      auto const& strawId() const { return straw_.id(); }
      auto const& strawHitIndex() const { return shindex_; }
      double minDOCA() const { return mindoca_; }
      double strawRadius() const { return rstraw_; }
      bool nullUpdater() const { return null_; }
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      WireHitState whstate_; // current state
      Line wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
      // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
      // (when multiplied by the propagation velocity).
      CA ptca_; // reference time and position of closest approach to the wire; this is generally biased by the hit
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
      double mindoca_; // minimum doca: used in variance and offset for null hits
      double rstraw_; // straw radius
      ComboHit const& chit_; // reference to hit
      StrawHitIndex shindex_; // index to the StrawHit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
      bool null_; // record if null updater was used on this hit
      // uitility functions
      void updateResiduals(WireHitState const& whstate);
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate,
      StrawProperties const& sprops,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    bfield_(bfield), whstate_(whstate), wire_(ptca.sensorTraj()),
    ptca_(ptca.localTraj(),wire_,ptca.precision(),ptca.tpData(),ptca.dDdP(),ptca.dTdP()),
    mindoca_(sprops.strawInnerRadius()), rstraw_(sprops.strawInnerRadius()),
    chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse), null_(false)
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
    if(!ptca_.usable())throw std::runtime_error("WireHit TPOCA failure");
  }

  template <class KTRAJ> double KKStrawHit<KTRAJ>::nullVariance(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case dresid: default:
        return (mindoca_*mindoca_)/3.0; // doca is signed
      case tresid:
        return (mindoca_*mindoca_)/(dinfo.vdrift_*dinfo.vdrift_*12.0); // TOCA is always larger than the crossing time
    }
  }

  template <class KTRAJ> double KKStrawHit<KTRAJ>::nullOffset(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case dresid: default:
        return 0.0; // not sure if there's a better answer
      case tresid:
        return -0.5*mindoca_/dinfo.vdrift_;
    }
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    WireHitState whstate = this->hitState();
    if(first){
      // look for an updater; if it's there, update the state
      auto dshu = miconfig.findUpdater<DOCAStrawHitUpdater>();
      auto nshu = miconfig.findUpdater<NullStrawHitUpdater>();
      if(nshu != 0 && dshu != 0)throw std::invalid_argument(">1 StrawHit updater specified");
      if(nshu != 0 || dshu != 0){
        // compute the unbiased closest approach; this is brute force, but works
        auto uca = this->unbiasedClosestApproach();
        if(uca.usable()){
          if(nshu != 0){
            null_ = true;
            mindoca_ = strawRadius();
            whstate = nshu->wireHitState(uca.doca());
          } else if(dshu != 0){
            null_ = false;
            // update minDoca (for null ambiguity error estimate)
            mindoca_ = std::min(dshu->minDOCA(),strawRadius());
            whstate = dshu->wireHitState(uca.doca());
          }
        } else {
          whstate = WireHitState();
        }
      }
    }
    // update residuals
    this->updateResiduals(whstate);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::updateResiduals(WireHitState const& whstate) {
    // update the state
    whstate_ = whstate;
    if(whstate_.active()){
      // compute drift parameters.  These are used even for null-ambiguity hits
      VEC3 bvec = bfield_.fieldVect(ptca_.particlePoca().Vect());
      auto pdir = bvec.Cross(wire_.direction()).Unit(); // direction perp to wire and BFieldMap
      VEC3 dvec = ptca_.delta().Vect();
      double phi = asin(double(dvec.Unit().Dot(pdir))); // azimuth around the wire WRT the BField
      //use reference DOCA to set drift info, as that's what the residual expects (relative to reference parameters)
      POL2 drift(fabs(ptca_.doca()), phi);
      DriftInfo dinfo;
      distanceToTime(drift, dinfo);
      if(whstate_.useDrift()){
        // translate PCA to residual. Use ambiguity to convert drift time to a time difference.
        double dsign = whstate_.lrSign()*ptca_.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
        double dt = ptca_.deltaT()-dinfo.tdrift_*dsign;
        // time differnce affects the residual both through the drift distance (DOCA) and the particle arrival time at the wire (TOCA)
        DVEC dRdP = ptca_.dDdP()*dsign/dinfo.vdrift_ - ptca_.dTdP();
        rresid_[tresid] = Residual(dt,dinfo.tdriftvar_,0.0,true,dRdP);
        rresid_[dresid] = Residual();
      } else {
        // interpret DOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
        DVEC dRdP = -ptca_.lSign()*ptca_.dDdP();
        double dd = ptca_.doca() + nullOffset(dresid,dinfo);
        double nulldvar = nullVariance(dresid,dinfo);
        rresid_[dresid] = Residual(dd,nulldvar,0.0,true,dRdP);
        //  interpret TOCA as a residual
        double dt = ptca_.deltaT() + nullOffset(tresid,dinfo);
        // the time constraint variance is the sum of the variance from maxdoca and from the intrinsic measurement variance
        double nulltvar = dinfo.tdriftvar_ + nullVariance(tresid,dinfo);
        rresid_[tresid] = Residual(dt,nulltvar,0.0,true,-ptca_.dTdP());
        // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
      }
    } else {
      // inactive residuals
      rresid_[tresid] = rresid_[dresid] = Residual();
    }
  }

  template <class KTRAJ> Residual const& KKStrawHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >=2)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }


  // the purpose of this class is to allow computing the drift using calibrated quantities (StrawResponse)
  template <class KTRAJ> void KKStrawHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    if(hitState().useDrift()){
      dinfo.vdrift_ = sresponse_.driftInstantSpeed(straw_.id(),drift.R(),drift.Phi());
      dinfo.tdrift_ = sresponse_.driftDistanceToTime(straw_.id(),drift.R(),drift.Phi());
      double tderr = sresponse_.driftTimeError(straw_.id(),drift.R(),drift.Phi(),drift.R());
      dinfo.tdriftvar_ = tderr*tderr;
    } else {
      if(null_){
        dinfo.vdrift_ = sresponse_.driftConstantSpeed();
        dinfo.tdrift_ = chit_.driftTime(); // use TOT to estimate drift time
        dinfo.tdriftvar_ =  50.0; // temporary hack FIXME!
      } else {
        dinfo.vdrift_ = sresponse_.driftInstantSpeed(straw_.id(),0.5*mindoca_,drift.Phi());
        dinfo.tdrift_ = 0.5*mindoca_/dinfo.vdrift_; // average drift time, assuming a flat distribution
        dinfo.tdriftvar_ = (mindoca_*mindoca_)/(dinfo.vdrift_*dinfo.vdrift_*12.0); //
      }
    }
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
      if(rresid_[tresid].active())
        ost << " Active Time Residual " << rresid_[tresid];
      if(rresid_[dresid].active())
        ost << " Active Distance Residual " << rresid_[dresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << ptca_.tpData() << std::endl;
    }
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
