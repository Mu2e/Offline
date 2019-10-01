//
//  Persistent representation of the BTrk Kalman filter fit (KalRep)
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_KalSeed_HH
#define RecoDataProducts_KalSeed_HH
// mu2e
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "RecoDataProducts/inc/TrkCaloHitSeed.hh"
#include "RecoDataProducts/inc/TrkStraw.hh"
#include "RecoDataProducts/inc/KalSegment.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "canvas/Persistency/Common/Ptr.h"
// BTrk
#include "BTrk/TrkBase/TrkParticle.hh"
// Root
#include <Rtypes.h>
// C++
#include <vector>
namespace mu2e {
  class CaloCluster;
  struct KalSeed {
    KalSeed() :  _flt0(0), _chisq(0.0), _fitcon(0.0) {}
    KalSeed(TrkParticle tpart,TrkFitDirection fdir,TrkT0 const& t0, double flt0, TrkFitFlag const& status) :
      _tpart(tpart), _fdir(fdir), _status(status), _t0(t0), _flt0(static_cast<Float_t>(flt0)),
      _chisq(-1.0), _fitcon(-1.0), _nbend(0) {}

    TrkParticle const& particle() const { return _tpart; }
    TrkFitDirection const& fitDirection() const { return _fdir; }
    std::vector<TrkStrawHitSeed> const& hits() const { return _hits;}
    TrkCaloHitSeed const& caloHit() const { return _chit; }
    std::vector<TrkStraw> const& straws() const { return _straws;}
    std::vector<KalSegment> const& segments() const { return _segments; }
    // find the nearest segment to a given GLOBAL flightlength
    std::vector<KalSegment>::const_iterator nearestSegment(float fltlen)  const;
    std::vector<KalSegment>::const_iterator nearestSegment(const XYZVec& pos)  const; // find nearest segment to a GLOBAL position
    TrkFitFlag const& status() const { return _status; }
    Float_t flt0() const { return _flt0; }
    HitT0 const& t0() const { return _t0; }
    Float_t chisquared() const { return _chisq; }
    Float_t fitConsistency() const { return _fitcon; }
    UInt_t nBend() const { return _nbend; }
    bool hasCaloCluster() const { return _chit.caloCluster().isNonnull(); }
    art::Ptr<CaloCluster> const& caloCluster() const { return _chit.caloCluster(); }
    art::Ptr<HelixSeed> const& helix() const { return _helix; }
    art::Ptr<KalSeed> const& kalSeed() const { return _kal; }

    // global information about the track
    TrkParticle			    _tpart; // particle assumed for this fit
    TrkFitDirection	      	    _fdir; // direction in which this particle was fit
    TrkFitFlag			    _status; // status of this fit
    HitT0			    _t0; // track t0; Time particle crosses z=0
    Float_t			    _flt0; // flight distance where the track crosses the tracker midplane (z=0)
    Float_t			    _chisq; // fit chisquared value
    Float_t			    _fitcon; // fit consistency
    UInt_t			    _nbend; // # of field corrections
    //
    // contained content substructure.
    //
    std::vector<KalSegment>	    _segments; // segments of the Kalman filter fit result
    std::vector<TrkStrawHitSeed>    _hits; // hit seeds for all the hits used in this fit
    std::vector<TrkStraw>	    _straws; // straws interesected by this fit
    TrkCaloHitSeed		    _chit;  // CaloCluster-based hit.  If it has no CaloCluster, this is unused
    // add content for BField correction information FIXME!
    // 
    // referenced content.
    //
    art::Ptr<HelixSeed>             _helix; // associated Helix Seed (for seed fits); can be null
    art::Ptr<KalSeed>               _kal; // associated Kalman Seed (for final fits); can be null
  };
  typedef std::vector<mu2e::KalSeed> KalSeedCollection;
}
#endif
