//
//  Persistent representation of the BTrk Kalman filter fit (KalRep)
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_KalSeed_HH
#define RecoDataProducts_KalSeed_HH
// mu2e
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkStraw.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
// Root
#include <Rtypes.h>
// C++
#include <vector>
#include <memory>
namespace mu2e {
  class CaloCluster;
  struct KalSeed {
    using LHPT = KinKal::ParticleTrajectory<KinKal::LoopHelix>;
    using CHPT = KinKal::ParticleTrajectory<KinKal::CentralHelix>;
    using KLPT = KinKal::ParticleTrajectory<KinKal::KinematicLine>;
    using LHPTPtr = std::unique_ptr<LHPT>;
    using CHPTPtr = std::unique_ptr<CHPT>;
    using KLPTPtr = std::unique_ptr<KLPT>;
    KalSeed() {}
    KalSeed(PDGCode::type tpart,TrkFitDirection fdir, TrkFitFlag const& status, double flt0=0.0 ) :
      _tpart(tpart), _fdir(fdir), _status(status), _flt0(static_cast<Float_t>(flt0)){}

    PDGCode::type particle() const { return _tpart; }
    TrkFitDirection const& fitDirection() const { return _fdir; }
    std::vector<TrkStrawHitSeed> const& hits() const { return _hits;}
    TrkCaloHitSeed const& caloHit() const { return _chit; }
    std::vector<TrkStraw> const& straws() const { return _straws;}
    std::vector<KalSegment> const& segments() const { return _segments; }
    TrkFitFlag const& status() const { return _status; }
    double t0Val() const;
    Float_t chisquared() const { return _chisq; }
    Float_t fitConsistency() const { return _fitcon; }
    UInt_t nTrajSegments() const { return _nseg; }
    bool hasCaloCluster() const { return _chit.caloCluster().isNonnull(); }
    art::Ptr<CaloCluster> const& caloCluster() const { return _chit.caloCluster(); }
    std::vector<KalSegment>::const_iterator nearestSeg(double time)  const;
    bool loopHelixFit() const { return _status.hasAllProperties(TrkFitFlag::KKLoopHelix); }
    bool centralHelixFit() const { return _status.hasAllProperties(TrkFitFlag::KKCentralHelix); }
    bool kinematicLineFit() const { return _status.hasAllProperties(TrkFitFlag::KKLine); }
    bool seedBTrkFit() const { return _status.hasAllProperties(TrkFitFlag::KSF); }
    bool finalBTrkFit() const { return _status.hasAllProperties(TrkFitFlag::KFF); }
    // reconstitute (as best as possible) the fit trajectory.  The ptr will be null if the fit wasn't based on the requested trajector type
    // Note these return by value
    // Note that the returned piecetraj may have large gaps, unless the full fit trajectory was stored in the seed.
    LHPTPtr loopHelixFitTrajectory() const;
    CHPTPtr centralHelixFitTrajectory() const;
    KLPTPtr kinematicLineFitTrajectory() const;

    // global information about the track
    PDGCode::type     _tpart = PDGCode::unknown; // particle assumed for this fit
    TrkFitDirection   _fdir = TrkFitDirection::downstream; // direction in which this particle was fit
    TrkFitFlag        _status; // status of this fit: includes alglorithm information
    Float_t           _chisq = -1; // fit chisquared value
    Float_t           _fitcon = -1; // fit consistency
    UInt_t            _nseg = 0; // # of fit trajectory segments
    Float_t           _maxgap = 0;
    Float_t           _avggap = 0; // information about trajectory gaps
    //
    // contained content substructure.
    //
    std::vector<KalSegment>     _segments; // segments of the Kalman filter fit result
    std::vector<TrkStrawHitSeed>    _hits; // hit seeds for all the hits used in this fit
    std::vector<TrkStraw>     _straws; // straws interesected by this fit
    TrkCaloHitSeed        _chit;  // CaloCluster-based hit.  If it has no CaloCluster, this has no content
    std::vector<KalSegment>::const_iterator nearestSegment(float time)  const;
    //
    // deprecated BTrk legacy content, DO NOT write any new code which depends on these functions
    // find the nearest segment to a given the time
    std::vector<KalSegment>::const_iterator nearestSegmentFlt(float fltlen)  const;
    std::vector<KalSegment>::const_iterator nearestSegment(const XYZVectorF& pos)  const; // find nearest segment to a GLOBAL position
    Float_t flt0() const { return _flt0; }
    HitT0 t0() const;
    Float_t         _flt0 = 0.0; // flight distance where the track crosses the tracker midplane (z=0).  Redundant with t0 in KinKal fits, and in the wrong unit
  };
  typedef std::vector<mu2e::KalSeed> KalSeedCollection;
}
#endif
