#ifndef TrkCaloHit_HH
#define TrkCaloHit_HH
// BaBar
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "RecoDataProducts/inc/HitT0.hh"
// Mu2e
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/HitT0.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// c++
#include <vector>
// forward refs
class TrkDifTraj;

namespace mu2e
{
  class TrkCaloHit : public TrkHit {
  public:
    enum TrkStrawHitFlag {weededHit=-5, driftFail=-3, updateFail=-1,addedHit=3,unweededHit=4};
    TrkCaloHit(const CaloCluster& caloCluster, const HitT0& trkt0, double fltlen, double exterr);
    virtual ~TrkCaloHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    virtual int ambig() const { return _iamb; }
//    virtual void invert();

    double hitRMS() const { return _hitRMS;}
// caloCluster specific interface
    const CaloCluster& caloCluster() const { return _caloCluster; }
// correct the hit time for the wire propagation
    double time() const;
    void hitPosition(CLHEP::Hep3Vector& hpos) const;
    HitT0 const& hitT0() const { return _hitt0;}
    void updateHitT0(HitT0 const& t0) { _hitt0 = t0; }
// external hit error (mm); 
    double extErr() const { return _exterr; }
// error ON RDrift and residual coming from hit t0 error
    double t0Err() const { return _hitt0._t0err;}
// total error
    double totalErr() const { return _toterr; }
// intrinsic hit error (mm)
    double hitErr() const { return _hitErr; }
// changing the extneral hit error invalidates the cache, it should invalidate the fit, FIXME!!!!
    void setExtErr(double exterr) { _exterr = exterr; }
// logical operators to allow searching for CaloHits
    // bool operator == (CaloCluster const& cluster) const { return _caloCluster == cluster; }
    // bool operator != (CaloCluster const& cluster) const { return !operator==(cluster); }
    void print(std::ostream& ) const;
  protected:
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);
  //private:
    const CaloCluster& _caloCluster;
    TrkLineTraj*       _hittraj;
    HitT0              _hitt0;
    double             _exterr, _toterr;
    int                _iamb;
    double             _tddist;
    double             _tddist_err;
    double             _hitRMS, _hitErr;
    CLHEP::Hep3Vector  _caloClusterPos;
  };
/// / unary functor to select TrkStrawHit from a given hit
 //  struct FindTrkCaloHit {
 //    FindTrkStrawHit(CaloCluster const& caloCluster) : _caloCluster(caloCluster) {}
 //    bool operator () (TrkCaloHit* const& tcal ) { return tcal->caloCluster() == _caloCluster; }
 //    CaloCluster const& _caloCluster;
 //  };
// define TrkStrawHitVector, to allow explicit conversion and construction
  typedef std::vector<TrkCaloHit*> TrkCaloHitVector;
// utility function to convert vector of TrkHits into TrkCaloHits
  void convert(TrkHitVector const& thv, TrkCaloHitVector& tshv);
}

#endif
