#ifndef Mu2e_KKCaloHit_hh
#define Mu2e_KKCaloHit_hh
//
//  hit representing a time measurement from a calorimeter cluster
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
// mu2e includes
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include <stdexcept>
namespace mu2e {
  using KinKal::Line;
  using KinKal::MetaIterConfig;
  using KinKal::Residual;
  using KinKal::CAHint;
  using KinKal::ClosestApproachData;
  using CCPtr = art::Ptr<CaloCluster>;
  template <class KTRAJ> class KKCaloHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = KinKal::ClosestApproach<KTRAJ,Line>;
      using HIT = KinKal::Hit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      // Hit interface overrrides
      unsigned nResid() const override { return 1; } // 1 time residual
      Residual const& refResidual(unsigned ires=0) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      KTRAJPTR const& refTrajPtr() const override { return tpca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // scintHit explicit interface
      KKCaloHit(CCPtr caloCluster,  PCA const& pca, double tvar, double wvar);
      virtual ~KKCaloHit(){}
      Residual const& timeResidual() const { return rresid_; }
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      auto const& sensorAxis() const { return saxis_; }
      auto const& closestApproach() const { return tpca_; }
      auto timeVariance() const { return tvar_; }
      auto widthVariance() const { return wvar_; }
      auto const& caloCluster() const { return caloCluster_; }
      auto precision() const { return tpca_.precision(); }

    private:
      CCPtr caloCluster_;  // associated calorimeter cluster
      Line saxis_; // axis along the crystals, through the COG
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      bool active_; // active or not
      CA tpca_; // reference time and distance of closest approach to the axis.
      // caches
      Residual rresid_; // residual WRT most recent reference parameters
  };

  template <class KTRAJ> KKCaloHit<KTRAJ>::KKCaloHit(CCPtr caloCluster,  PCA const& pca, double tvar, double wvar) :    caloCluster_(caloCluster), saxis_(pca.sensorTraj()), tvar_(tvar), wvar_(wvar), active_(true),
    tpca_(pca.localTraj(),saxis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()) {
    }

  template <class KTRAJ> Residual const& KKCaloHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void KKCaloHit<KTRAJ>::updateReference(std::shared_ptr<KTRAJ> const& ktrajptr) {
    // compute PCA
    CAHint tphint( saxis_.t0(), saxis_.t0());
    // don't update the hint: initial T0 values can be very poor, which can push the CA calculation onto the wrong helix loop,
    // from which it's impossible to ever get back to the correct one.  Active loop checking might be useful eventually too TODO
    //    if(tpca_.usable()) tphint = CAHint(tpca_.particleToca(),tpca_.sensorToca());
    tpca_ = CA(ktrajptr,saxis_,tphint,tpca_.precision());
    if(!tpca_.usable())throw std::runtime_error("TPOCA failure");
  }


  template <class KTRAJ> void KKCaloHit<KTRAJ>::updateState(MetaIterConfig const& config,bool first) {
    // check that TPCA position is consistent with the physical sensor. This can be off if the CA algorithm finds the wrong helix branch
    // early in the fit when t0 has very large errors.
    // If it is unphysical try to adjust it back using a better hint.
    auto ppos = tpca_.particlePoca().Vect();
    auto sstart = saxis_.startPosition();
    auto send = saxis_.endPosition();
    double slen = (send-sstart).R();
    // tolerance should come from the config.  Should also test relative to the error. FIXME
    double tol = slen*1.0;
    if( (ppos-sstart).Dot(saxis_.direction()) < -tol ||
        (ppos-send).Dot(saxis_.direction()) > tol) {
      // adjust hint to the middle and try agian
      double sspeed = tpca_.particleTraj().velocity(tpca_.particleToca()).Dot(saxis_.direction());
      double sdist = (ppos - saxis_.position3(saxis_.range().mid())).Dot(saxis_.direction());
      auto tphint = tpca_.hint();
      tphint.particleToca_ -= sdist/sspeed;
      tpca_ = CA(tpca_.particleTrajPtr(),saxis_,tphint,precision());
      // should check if this is still unphysical and disable the hit if so FIXME
    }
    // residual is just delta-T at CA.
    // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
    double dd2 = tpca_.dirDot()*tpca_.dirDot();
    double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
    rresid_ = Residual(tpca_.deltaT(),totvar,0.0,true,-tpca_.dTdP());
  }

  template<class KTRAJ> void KKCaloHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->active())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " KKCaloHit time " << this->time() << " tvar " << tvar_ << " wvar " << wvar_ << std::endl;
    if(detail > 0){
      ost << "Axis ";
      saxis_.print(ost,detail);
    }
  }

}
#endif
