#ifndef Mu2eKinKal_KKStrawXing_hh
#define Mu2eKinKal_KKStrawXing_hh
//
//  StrawXing using Mu2e-specific StrawMaterial class
//
#include "KinKal/Detector/ElementXing.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::SVEC3;
  using KinKal::DVEC;
  using KinKal::CAHint;
  using KinKal::DPDV;
  using KinKal::MomBasis;
  using KinKal::NParams;
  using KinKal::SensorLine;
  template <class KTRAJ> class KKStrawXing : public KinKal::ElementXing<KTRAJ> {
    public:
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = KinKal::ClosestApproach<KTRAJ,SensorLine>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      // construct with closest approach
      KKStrawXing(KKSTRAWHITPTR const& strawhit, CA const& ca, KKStrawMaterial const& smat, Straw const& straw,bool active=false);
      // construct with an associated StrawHit
      KKStrawXing(KKSTRAWHITPTR const& strawhit, KKStrawMaterial const& smat);
      virtual ~KKStrawXing() {}
      // clone op for reinstantiation
      KKStrawXing(KKStrawXing const& rhs):
          KKStrawXing(
            rhs.strawHitPtr(),
            rhs.closestApproach(),
            rhs.strawMaterial(),
            rhs.straw()
          ){
        auto shptr = rhs.strawHitPtr();
        if (shptr){
          this->setStrawHitPtr(rhs.strawHitPtr());
        }
      }
      std::shared_ptr< KinKal::ElementXing<KTRAJ> > clone(CloneContext& context) const override{
        auto rv = std::make_shared< KKStrawXing<KTRAJ> >(*this);

        // point to new instance of partner hit, if not null
        KKSTRAWHITPTR shptr;
        if (shptr_){
          shptr = context.get(shptr_);
        }
        rv->setStrawHitPtr(shptr);

        // point to new instance of ClosestApproach
        auto ca = rv->closestApproach();
        auto trajectory = std::make_shared<KTRAJ>(ca.particleTraj());
        ca.setTrajectory(trajectory);
        rv->setClosestApproach(ca);

        return rv;
      };
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      // offset time WRT TOCA to avoid exact overlapp with the wire hit.  Note: the offset must be POSITIVE to insure
      // Xing is updated after the associated hit
      double time() const override { return ca_.particleToca() + toff_; }
      double transitTime() const override; // time to cross this element
      KTRAJ const& referenceTrajectory() const override { return ca_.particleTraj(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      bool active() const override;
      // accessors
      auto const& config() const { return sxconfig_; }
      auto const& closestApproach() const { return ca_; }
      auto const& strawMaterial() const { return smat_; }
      auto precision() const { return ca_.precision(); }
      auto const& straw() const { return straw_; }
      auto const& strawId() const { return straw_.id(); }
      auto const& strawHitPtr() const { return shptr_; }
      // allow overwriting the straw hit
      void setHit(KKSTRAWHITPTR& shptr) { shptr_ = shptr; }
    private:
      StrawXingUpdater sxconfig_; // cache of most recent update
      KKSTRAWHITPTR shptr_; // reference to associated StrawHit
      SensorLine axis_; // straw axis, expressed as a timeline
      KKStrawMaterial const& smat_;
      Straw const& straw_; // reference to straw object, to allow follolwing geometry
      CA ca_; // result of most recent TPOCA
      double toff_; // small time offset
      std::vector<MaterialXing> mxings_;
      Parameters fparams_; // parameter change for forwards time
      bool active_; // inside active region or not
      double varscale_; // variance scale
      // modifiers to support cloning
      void setStrawHitPtr(KKSTRAWHITPTR ptr) { shptr_ = ptr; }
      void setClosestApproach(const CA& ca){ ca_ = ca; }
  };

  template <class KTRAJ> KKStrawXing<KTRAJ>::KKStrawXing(KKSTRAWHITPTR const& strawhit, CA const& ca, KKStrawMaterial const& smat, Straw const& straw,bool active) :
    shptr_(strawhit),
    axis_(ca.sensorTraj()),
    smat_(smat),
    straw_(straw),
    ca_(ca.particleTraj(),axis_,ca.hint(),ca.precision()),
    toff_(smat.wireRadius()/ca.particleTraj().speed(ca.particleToca())), // locate the effect to 1 side of the wire to avoid overlap with hits
    active_(active),
    varscale_(1.0)
  {}

  template <class KTRAJ> KKStrawXing<KTRAJ>::KKStrawXing(KKSTRAWHITPTR const& strawhit, KKStrawMaterial const& smat) :
    shptr_(strawhit),
    axis_(Mu2eKinKal::strawLine(strawhit->straw(),strawhit->closestApproach().particleToca())),
    smat_(smat),
    straw_(strawhit->straw()),
    ca_(strawhit->closestApproach().particleTraj(),axis_,strawhit->closestApproach().hint(),strawhit->closestApproach().precision()),
    toff_(smat.wireRadius()/strawhit->closestApproach().particleTraj().speed(strawhit->closestApproach().particleToca())),
    active_(strawhit->active()),
    varscale_(1.0)
  {}

  template <class KTRAJ> bool KKStrawXing<KTRAJ>::active() const {
    // if the associated hit is active, use it's state. Otherwise use the intrinsic state
    return active_ || (shptr_ && shptr_->hitState().active());
  }

  template <class KTRAJ> void KKStrawXing<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    CAHint tphint(axis_.timeAtMidpoint(),axis_.timeAtMidpoint());
    if(ca_.usable()){
      tphint = ca_.hint();
    }else if(shptr_ && shptr_->closestApproach().usable()){
      tphint = shptr_->closestApproach().hint();
    }
    PCA pca(ptraj,axis_,tphint,precision());
    ca_ = static_cast<CA>(pca);
  }

  template <class KTRAJ> Parameters KKStrawXing<KTRAJ>::params() const {
    return fparams_;
  }

  template <class KTRAJ> void KKStrawXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // reset
    mxings_.clear();
    if(first) {
      // search for an updater among this meta-iteration payload
      auto sxconfig = miconfig.findUpdater<StrawXingUpdater>();
      // cache, as this also sets parameters used in calculating path lengths
      if(sxconfig != 0)sxconfig_ = *sxconfig;
      // require a validat updater
      if(sxconfig_.nsig_ < 0)throw cet::exception("RECO")<<"mu2e::KKStrawXing: invalid updater!" <<  std::endl;
      // update the DOCA range test
      if(sxconfig_.maxdoca_ > 0.0)active_ = fabs(ca_.tpData().doca()) < sxconfig_.maxdoca_;
      // set the variance scale (temperature)
      if(sxconfig_.scalevar_)
        varscale_ = miconfig.varianceScale();
      else
        varscale_ = 1.0;
    }
    // update the material xings from gas, straw wall, and wire
    smat_.findXings(ca_.tpData(),sxconfig_,mxings_);
    // update the effect these have on the parameters
    fparams_ = this->parameterChange(varscale_);
  }

  template <class KTRAJ> double KKStrawXing<KTRAJ>::transitTime() const {
    return smat_.transitLength(ca_.tpData())/ca_.particleTraj().speed(ca_.particleToca());
  }

  template <class KTRAJ> void KKStrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"KKStraw Xing time " << this->time();
    if(detail > 0){
      for(auto const& mxing : mxings_){
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
