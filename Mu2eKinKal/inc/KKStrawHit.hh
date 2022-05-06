#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
// mu2eKinKal classes
//KinKal classes
#include "Offline/Mu2eKinKal/inc/KKStrawHitUpdater.hh"
#include "KinKal/Detector/WireHit.hh"
// Mu2e-specific classes
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
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
  class KKStrawHitUpdater;

  template <class KTRAJ> class KKStrawHit : public KinKal::WireHit<KTRAJ> {
    public:
      using WIREHIT = KinKal::WireHit<KTRAJ>;
      using Dimension = typename WIREHIT::Dimension;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&, double rstraw,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // WireHit and Hit interface implementations
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      double nullVariance(Dimension dim,DriftInfo const& dinfo) const override;
      double nullOffset(Dimension dim,DriftInfo const& dinfo) const override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific to KKStrawHit: this has a constant drift speed
      virtual ~KKStrawHit(){}
      // accessors
      ComboHit const& hit() const { return chit_; }
      Straw const& straw() const { return straw_; }
      StrawHitIndex const& strawHitIndex() const { return shindex_; }
      double minDOCA() const { return mindoca_; }
      double strawRadius() const { return rstraw_; }
    private:
      double mindoca_; // minimum doca: used in variance and offset for null hits
      double rstraw_; // straw radius
      ComboHit const& chit_; // reference to hit
      StrawHitIndex shindex_; // index to the StrawHit
      Straw const& straw_; // reference to straw of this hit
      StrawResponse const& sresponse_; // straw calibration information
      friend class KKStrawHitUpdater;
  };

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate, double rstraw,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    WIREHIT(bfield,ptca,whstate), mindoca_(rstraw), rstraw_(rstraw), chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse)
  {
    // make sure this is a single-straw based ComboHit
    if(chit_.mask().level() != StrawIdMask::uniquestraw)
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: ComboHit doesn't correspond to a unique straw"<< std::endl;
  }

  template <class KTRAJ> double KKStrawHit<KTRAJ>::nullVariance(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case WIREHIT::dresid: default:
        return (mindoca_*mindoca_)/3.0; // doca is signed
      case WIREHIT::tresid:
        return (mindoca_*mindoca_)/(dinfo.vdrift_*dinfo.vdrift_*12.0); // TOCA is always larger than the crossing time
    }
  }

  template <class KTRAJ> double KKStrawHit<KTRAJ>::nullOffset(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case WIREHIT::dresid: default:
        return 0.0; // not sure if there's a better answer
      case WIREHIT::tresid:
        return -0.5*mindoca_/dinfo.vdrift_;
    }
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    WIREHIT::update(pktraj);
  }



  template <class KTRAJ> void KKStrawHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    WIREHIT::update(pktraj,miconfig);
    // look for an updater; if it's there, update the state
    auto kkshu = miconfig.findUpdater<KKStrawHitUpdater>();
    if(kkshu != 0){
      kkshu->update(*this);
      WIREHIT::update(pktraj,miconfig);
    }
  }

  // the purpose of this class is to allow computing the drift using calibrated quantities (StrawResponse)
  template <class KTRAJ> void KKStrawHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // for now, use simplified response model.
    dinfo.vdrift_ = sresponse_.driftConstantSpeed();
    dinfo.tdrift_ = drift.R()/dinfo.vdrift_;
    dinfo.tdriftvar_ = 16.0; // temporary hack FIXME
    // std::cout << "tdrift " << dinfo.tdrift_ << " VDrift = "<< dinfo.vdrift_ << " derr " << derr << " tvar " << dinfo.tdriftvar_ << std::endl;
  }

  template<class KTRAJ> void KKStrawHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " KKStrawHit time " << this->time();
    WIREHIT::print(ost,detail);
    if(detail > 0)chit_.print(ost,true);
  }
}
#endif
