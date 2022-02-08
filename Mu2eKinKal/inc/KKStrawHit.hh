#ifndef Mu2eKinKal_KKStrawHit_hh
#define Mu2eKinKal_KKStrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
// mu2eKinKal classes
#include "Offline/Mu2eKinKal/inc/KKStrawHitUpdater.hh"
//KinKal classes
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

  template <class KTRAJ> class KKStrawHit : public KinKal::WireHit<KTRAJ> {
    public:
      using WIREHIT = KinKal::WireHit<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&, double mindoca,
          ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse);
      // WireHit and Hit interface implementations
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
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

  template <class KTRAJ> KKStrawHit<KTRAJ>::KKStrawHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate, double mindoca,
      ComboHit const& chit, Straw const& straw, StrawHitIndex const& shindex, StrawResponse const& sresponse) :
    WIREHIT(bfield,ptca,whstate,mindoca), chit_(chit), shindex_(shindex), straw_(straw), sresponse_(sresponse)
  {
    // make sure this is a single-straw based ComboHit
    if(chit_.mask().level() != StrawIdMask::uniquestraw)
      throw cet::exception("RECO")<<"mu2e::KKStrawHit: ComboHit doesn't correspond to a unique straw"<< endl;
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    WIREHIT::update(pktraj);
  }

  template <class KTRAJ> void KKStrawHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    using KinKal::ClosestApproachData;
    using KinKal::WireHitState;
    PTCA tpoca = WIREHIT::wirePTCA(pktraj);
    if(tpoca.usable()){
      this->tpdata_ = tpoca.tpData();
      this->setRefParams(pktraj.nearestPiece(tpoca.particleToca()));
      bool updated(false);
      for(auto const& updater : miconfig.updaters_){
        auto const shupdater = std::any_cast<KKStrawHitUpdater>(&updater);
        if(shupdater != 0){
          if(updated) throw std::invalid_argument("Multiple updaters found");
          // update the null mindoca
          this->mindoca_ = std::min(shupdater->minddoca_,2.5);  // FIXME!
          // update the internal hit state (activity, LR ambiguity, intrinsic error, ...)
          auto const& poca = this->closestApproach();
          auto chisq = this->chisq();
          double absdoca = fabs(poca.doca());
          if( absdoca > shupdater->maxdoca_ || chisq.probability() < shupdater->minprob_){ // hit is too far from the wire or has too small a probability: disable it
            this->wstate_.state_ = WireHitState::inactive; // disable the hit
          } else if(absdoca > shupdater->minddoca_ && absdoca < shupdater->maxddoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
            this->wstate_.state_ = poca.doca() > 0.0 ? WireHitState::right : WireHitState::left;
          } else { // hit too close to the wire to resolve ambiguity, or with a suspiciously large drift: just use the raw wire position and time to constrain the track
            this->wstate_.state_ = WireHitState::null;
          }

          // now update the cache again in case the caches changed
          this->update(pktraj);
          WIREHIT::setResiduals(tpoca);
          updated = true;
       //   print(std::cout,1); // test FIXME
        }
      }
    } else {
      throw std::runtime_error("PTCA failure");
      // OK if no updater is found, hits may be frozen this meta-iteration
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
