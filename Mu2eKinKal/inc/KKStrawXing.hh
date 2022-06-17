#ifndef Mu2eKinKal_KKStrawXing_hh
#define Mu2eKinKal_KKStrawXing_hh
//
//  StrawXing using Mu2e-specific StrawMaterial class.  Otherwise it's the same as KinKal::StrawXing
//
#include "KinKal/Detector/ElementXing.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXingUpdater.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  template <class KTRAJ> class KKStrawXing : public KinKal::ElementXing<KTRAJ> {
    public:
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      // construct without an associated StrawHit
      KKStrawXing(PTCA const& ptca, KKStrawMaterial const& smat, StrawId sid);
      // construct with an associated StrawHit
      KKStrawXing(KKSTRAWHITPTR const& strawhit, KKStrawMaterial const& smat);
      virtual ~KKStrawXing() {}
      // ElementXing interface
      void updateReference(KTRAJPTR const& ktrajptr) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      // offset time WRT TOCA to avoid exact overlapp with the wire hit.  Note: the offset must be POSITIVE to insure
      // Xing is updated after the associated hit
      double time() const override { return tpca_.particleToca() + toff_; }
      double transitTime() const override; // time to cross this element
      KTRAJ const& referenceTrajectory() const override { return tpca_.particleTraj(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      auto const& closestApproach() const { return tpca_; }
      auto const& strawMaterial() const { return smat_; }
      auto const& config() const { return sxconfig_; }
      auto precision() const { return tpca_.precision(); }
      auto const& strawId() const { return sid_; }
    private:
      StrawId sid_; // StrawId
      KKSTRAWHITPTR shptr_; // reference to associated StrawHit
      KinKal::Line axis_; // straw axis, expressed as a timeline
      KKStrawMaterial const& smat_;
      CA tpca_; // result of most recent TPOCA
      double toff_; // small time offset
      KKStrawXingUpdater sxconfig_; // note this must come from an updater during processing
  };

  template <class KTRAJ> KKStrawXing<KTRAJ>::KKStrawXing(PCA const& pca, KKStrawMaterial const& smat, StrawId sid) :
    sid_(sid),
    axis_(pca.sensorTraj()),
    smat_(smat),
    tpca_(pca.localTraj(),axis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()),
    toff_(smat.wireRadius()/pca.particleTraj().speed(pca.particleToca())) // locate the effect to 1 side of the wire to avoid overlap with hits
  {}

  template <class KTRAJ> KKStrawXing<KTRAJ>::KKStrawXing(KKSTRAWHITPTR const& strawhit, KKStrawMaterial const& smat) :
    sid_(strawhit->straw().id()),
    shptr_(strawhit),
    axis_(strawhit->closestApproach().sensorTraj()),
    smat_(smat),
    tpca_(strawhit->closestApproach().particleTraj(),axis_,strawhit->closestApproach().hint(),strawhit->closestApproach().precision()),
    toff_(smat.wireRadius()/strawhit->closestApproach().particleTraj().speed(strawhit->closestApproach().particleToca()))
  {}

  template <class KTRAJ> void KKStrawXing<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    KinKal::CAHint tphint = tpca_.usable() ?  tpca_.hint() : KinKal::CAHint(axis_.range().mid(),axis_.range().mid());
    tpca_ = CA(ktrajptr,axis_,tphint,precision());
    if(!tpca_.usable())throw cet::exception("RECO")<<"mu2e::KKStrawXing: TPOCA failure" << std::endl;
 }

  template <class KTRAJ> void KKStrawXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    if(first) {
      // search for an update to the xing configuration among this meta-iteration payload
      auto sxconfig = miconfig.findUpdater<KKStrawXingUpdater>();
      if(sxconfig != 0){
        sxconfig_ = *sxconfig;
      }
      //  update the associated hit state
      if(shptr_)
        sxconfig_.hitstate_ = shptr_->hitState();
      else
        sxconfig_.hitstate_ = KinKal::WireHitState::inactive;
    }
    // find the material xings from gas, straw wall, and wire
    smat_.findXings(tpca_.tpData(),sxconfig_,EXING::matXings());
  }

  template <class KTRAJ> double KKStrawXing<KTRAJ>::transitTime() const {
    return smat_.transitLength(tpca_.tpData())/tpca_.particleTraj().speed(tpca_.particleToca());
  }

  template <class KTRAJ> void KKStrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"KKStraw Xing time " << this->time();
    if(detail > 0){
      for(auto const& mxing : this->matXings()){
        ost << " " << mxing.dmat_.name() << " pathLen " << mxing.plen_;
      }
    }
    if(detail > 1){
      ost << " Axis ";
      axis_.print(ost,0);
    }
    ost << std::endl;
  }
}
#endif
