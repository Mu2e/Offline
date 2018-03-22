#ifndef TrkCaloHit_HH
#define TrkCaloHit_HH
// BaBar
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkHit.hh"

#include "RecoDataProducts/inc/HitT0.hh"
// Mu2e
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// c++
#include <vector>
// forward refs
class TrkDifTraj;
class TrkDifPieceTraj;

namespace mu2e
{
  class TrkCaloHit : public TrkHit {
  public:
    TrkCaloHit(const CaloCluster& caloCluster, CLHEP::Hep3Vector &caloClusterPos, 
	       double crystalHalfLength,  CLHEP::Hep3Vector const& clusterAxis,
	       const HitT0& trkt0, double fltlen, double timeWeight, double _dtoffset);
    virtual ~TrkCaloHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }

// correct the hit time
    virtual double time                 () const;
    virtual void   hitPosition          (CLHEP::Hep3Vector& hpos) const;
    virtual bool   signalPropagationTime(double &propTime, double&Doca, 
				 double resid, double &residErr,
				 CLHEP::Hep3Vector trajDirection);//propagation time
    virtual void   trackT0Time          (double& htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt);
     // test the consistincy of this hit with 'physical' limts, with a given # of sigma
    virtual bool isPhysical         (double maxchi) const;

// intrinsic hit error (mm)
    double hitErr               () const { return _hitErr; }

    void print                  (std::ostream& ) const;
// caloCluster specific interface
    const CaloCluster& caloCluster() const { return _caloCluster; }
  protected:
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);

    const CaloCluster& _caloCluster;
    double             _dtoffset;
    TrkLineTraj*       _hittraj;
    double             _hitErr;
  };

// define TrkStrawHitVector, to allow explicit conversion and construction
  typedef std::vector<TrkCaloHit*> TrkCaloHitVector;
// utility function to convert vector of TrkHits into TrkCaloHits
  void convert(TrkHitVector const& thv, TrkCaloHitVector& tshv);
}

#endif
