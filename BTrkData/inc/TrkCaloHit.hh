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
    TrkCaloHit(CaloCluster const& caloCluster, CLHEP::Hep3Vector const& caloClusterPos,
	       double crystalLength,  CLHEP::Hep3Vector const& clusterAxis,
	       const HitT0& trkt0, double fltlen, double timeWeight, 
	       double hiterr, double terr, double _dtoffset);
    virtual ~TrkCaloHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }

// correct the hit time
    virtual double time                 () const;
//    virtual bool   time                 (HitT0&t0);
    virtual void   hitPosition          (CLHEP::Hep3Vector& hpos) const;
    virtual bool signalPropagationTime(	TrkT0& t0);  // this function should be const FIXME!!!
    // the followin function isn't used and should be removed FIXME!
    virtual void   trackT0Time          (double& htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt);
     // test the consistincy of this hit with 'physical' limts, with a given # of sigma
    virtual bool isPhysical         (double maxchi) const;

// intrinsic hit error (mm)
    double hitErr               () const { return _hitErr; }

    void print                  (std::ostream& ) const;
// caloCluster specific interface
    const CaloCluster& caloCluster() const { return _caloCluster; }
    double timeOffset() const { return _dtoffset; }
    double timeErr() const { return _tErr; }
  protected:
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);

    const CaloCluster& _caloCluster;
    double	       _clen; // crystal length
    double             _dtoffset;
    TrkLineTraj*       _hittraj;
    double             _hitErr; // geometric error on the cluster transverse position for POCA calculation
    double		_tErr; // error on the calorimeter time
  };

// define TrkStrawHitVector, to allow explicit conversion and construction
  typedef std::vector<TrkCaloHit*> TrkCaloHitVector;
// utility function to convert vector of TrkHits into TrkCaloHits
  void convert(TrkHitVector const& thv, TrkCaloHitVector& tshv);
}

#endif
