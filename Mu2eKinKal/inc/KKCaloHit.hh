#ifndef Mu2e_KKCaloHit_hh
#define Mu2e_KKCaloHit_hh
//
//  hit representing a time measurement from a calorimeter cluster
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
// mu2e includes
#include "RecoDataProducts/inc/CaloCluster.hh"
// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include <stdexcept>
namespace mu2e {
  using KinKal::Line;
  using KinKal::MetaIterConfig;
  using KinKal::Residual;
  using KinKal::CAHint;
  using CCPtr = art::Ptr<CaloCluster>;
  template <class KTRAJ> class KKCaloHit : public KinKal::ResidualHit<KTRAJ> {
    public:
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      // Hit interface overrrides
      unsigned nResid() const override { return 1; } // 1 time residual
      bool activeRes(unsigned ires=0) const override;
      Residual const& residual(unsigned ires=0) const override;
      double time() const override { return tpdata_.particleToca(); }
      void update(PKTRAJ const& pktraj) override;
      void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // scintHit explicit interface
      KKCaloHit(CCPtr caloCluster,  Line const& sensorAxis, PKTRAJ const& ptraj, double tvar, double wvar) : 
	caloCluster_(caloCluster), saxis_(sensorAxis), tvar_(tvar), wvar_(wvar), active_(true), precision_(1e-6) {
	  update(ptraj);
	}
      virtual ~KKCaloHit(){}
      Residual const& timeResidual() const { return rresid_; }
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      Line const& sensorAxis() const { return saxis_; }
      ClosestApproachData const& closestApproach() const { return tpdata_; }
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
      CCPtr const& caloCluster() const { return caloCluster_; }

    private:
      CCPtr caloCluster_;  // associated calorimeter cluster
      Line saxis_; // axis along the crystals, through the COG
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      bool active_; // active or not
      ClosestApproachData tpdata_; // reference time and distance of closest approach to the axis.
      // caches
      Residual rresid_; // residual WRT most recent reference parameters
      double precision_; // current precision
  };

  template <class KTRAJ> bool KKCaloHit<KTRAJ>::activeRes(unsigned ires) const {
    if(ires == 0 && active_)
      return true;
    else
      return false;
  }

  template <class KTRAJ> Residual const& KKCaloHit<KTRAJ>::residual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void KKCaloHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute PTCA
    CAHint tphint( saxis_.t0(), saxis_.t0());
    // don't update the hint: initial T0 values can be very poor, which can push the CA calculation onto the wrong helix loop,
    // from which it's impossible to ever get back to the correct one.  Active loop checking might be useful eventually too TODO
    //    if(tpdata_.usable()) tphint = CAHint(tpdata_.particleToca(),tpdata_.sensorToca());
    PTCA tpoca(pktraj,saxis_,tphint,precision_);
    if(tpoca.usable()){
      tpdata_ = tpoca.tpData();
      // residual is just delta-T at CA. 
      // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
      double dd2 = tpoca.dirDot()*tpoca.dirDot();
      double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
      rresid_ = Residual(tpoca.deltaT(),totvar,-tpoca.dTdP());
      this->setRefParams(pktraj.nearestPiece(tpoca.particleToca()));
    } else
      throw std::runtime_error("PTCA failure");
  }

  template <class KTRAJ> void KKCaloHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // for now, no updates are needed.  Eventually could test for consistency, update errors, etc
    precision_ = miconfig.tprec_;
    update(pktraj);
  }

  template<class KTRAJ> void KKCaloHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->active())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " KKCaloHit  tvar " << tvar_ << " wvar " << wvar_ << std::endl;
    if(detail > 0){
      ost << "Line ";
      saxis_.print(ost,detail);
    }
  }

}
#endif
