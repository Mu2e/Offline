#ifndef Mu2eKinKal_KKCaloHit_hh
#define Mu2eKinKal_KKCaloHit_hh
//
//  hit representing a time measurement from a calorimeter cluster.  The geometric information is used to
//  compute propagation time
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
// mu2e includes
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::SensorLine;
  using KinKal::MetaIterConfig;
  using KinKal::Residual;
  using KinKal::CAHint;
  using KinKal::ClosestApproachData;
  using CCPtr = art::Ptr<CaloCluster>;
  template <class KTRAJ> class KKCaloHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = KinKal::ClosestApproach<KTRAJ,SensorLine>;
      using HIT = KinKal::Hit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      // Hit interface overrrides
      unsigned nResid() const override { return 1; } // 1 time residual
      Residual const& refResidual(unsigned ires=0) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      KTRAJPTR const& refTrajPtr() const override { return tpca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // scintHit explicit interface
      KKCaloHit(CCPtr caloCluster,  PCA const& pca, double tvar, double wvar);
      virtual ~KKCaloHit(){}
      Residual const& timeResidual() const { return rresid_; }
      CA unbiasedClosestApproach() const;
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      auto const& sensorAxis() const { return saxis_; }
      auto const& closestApproach() const { return tpca_; }
      auto timeVariance() const { return tvar_; }
      auto widthVariance() const { return wvar_; }
      auto const& caloCluster() const { return caloCluster_; }
      auto precision() const { return tpca_.precision(); }

    private:
      CCPtr caloCluster_;  // associated calorimeter cluster
      SensorLine saxis_; // axis along the crystals, through the COG
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      CA tpca_; // reference time and distance of closest approach to the axis
      Residual rresid_; // residual WRT most recent reference parameters
  };

  template <class KTRAJ> KKCaloHit<KTRAJ>::KKCaloHit(CCPtr caloCluster,  PCA const& pca, double tvar, double wvar) :    caloCluster_(caloCluster), saxis_(pca.sensorTraj()), tvar_(tvar), wvar_(wvar),
    tpca_(pca.localTraj(),saxis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()) {
    }

  template <class KTRAJ> Residual const& KKCaloHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires !=0)throw cet::exception("RECO")<<"mu2e::KKCaloHit: Invalid residual" << std::endl;
    return rresid_;
  }

  template <class KTRAJ> void KKCaloHit<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // compute PCA
    CAHint tphint( saxis_.measurementTime(), saxis_.measurementTime());
    // don't update the hint: initial T0 values can be very poor, which can push the CA calculation onto the wrong helix loop,
    // from which it's impossible to ever get back to the correct one.  Active loop checking might be useful eventually too TODO
    //    if(tpca_.usable()) tphint = CAHint(tpca_.particleToca(),tpca_.sensorToca());
    auto ktrajptr = ptraj.nearestTraj(time());
    tpca_ = CA(ktrajptr,saxis_,tphint,tpca_.precision());
    if(!tpca_.usable())rresid_ = Residual(rresid_.value(),rresid_.variance(),0.0,false,rresid_.dRdP());
  }

  template <class KTRAJ> void KKCaloHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // check that TPCA position is consistent with the physical sensor. This can be off if the CA algorithm finds the wrong helix branch
    // early in the fit when t0 has very large errors.
    // If it is inconsistent with the crystal position try to adjust it back using a better hint.
    auto ppos = tpca_.particlePoca().Vect();
    auto sstart = saxis_.start();
    auto send = saxis_.end();
    // tolerance should come from the miconfig.  Should also test relative to the error. TODO
    double tol = 1.0*saxis_.length();
    double sdist = (ppos-sstart).Dot(saxis_.direction());
    double edist = (ppos-send).Dot(saxis_.direction());
    if( sdist < -tol || edist > tol) {
      // adjust hint to the middle of the crystal and try again
      double sspeed = tpca_.particleTraj().velocity(tpca_.particleToca()).Dot(saxis_.direction());
      auto tphint = tpca_.hint();
      tphint.particleToca_ -= 0.5*(sdist+edist)/sspeed;
      tphint.sensorToca_ = saxis_.timeAtMidpoint();
      tpca_ = CA(tpca_.particleTrajPtr(),saxis_,tphint,precision());
    }
    if(tpca_.usable()){
      // residual is just delta-T at CA.
      double totvar = tvar_;
      // the variance includes the intrinsic measurement variance and the tranvserse position resolution (which couples to the relative direction)
      double speed = saxis_.speed(tpca_.sensorToca()); // signal propagation speed along the crystal
      double s2 = speed*speed;
      double cost = tpca_.dirDot();  // cosine of angle between track and crystal axis
      double cost2 = cost*cost;
      double sint2 = 1.0-cost2;
      // protect against track co-linear with crystal
      if(sint2 > 1e-2) {
        totvar += cost2*wvar_/(s2*sint2);
      } else {
        double ldt = saxis_.length()/speed; // time to go the full length of the crystal
        double invvar2 = std::max(s2*sint2/(cost2*wvar_), 12/(ldt*ldt));
        totvar += 1.0/invvar2;
      }
      rresid_ = Residual(tpca_.deltaT(),totvar,0.0,true,tpca_.dTdP());
    } else {
      rresid_ = Residual(rresid_.value(),rresid_.variance(),0.0,false,rresid_.dRdP());
    }
    // finally update the weight
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> KinKal::ClosestApproach<KTRAJ,SensorLine> KKCaloHit<KTRAJ>::unbiasedClosestApproach() const {
    // compute the unbiased closest approach; this is brute force, but works
    auto const& ca = this->closestApproach();
    auto uparams = HIT::unbiasedParameters();
    KTRAJ utraj(uparams,ca.particleTraj());
    return CA(utraj,saxis_,ca.hint(),ca.precision());
  }

  template<class KTRAJ> void KKCaloHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(rresid_.active())
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
