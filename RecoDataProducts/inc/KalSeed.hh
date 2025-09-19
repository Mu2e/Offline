//
//  Persistent representation of the Kalman filter fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_KalSeed_HH
#define RecoDataProducts_KalSeed_HH
// mu2e
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitCalib.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkStraw.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
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
    using InterIter = std::vector<KalIntersection>::const_iterator;
    using InterIterCol = std::vector<InterIter>;
    KalSeed() {}
    KalSeed(PDGCode::type tpart, TrkFitFlag const& status, double flt0=0.0 ) :
      _tpart(tpart), _status(status), _flt0(static_cast<float>(flt0)){}

    PDGCode::type particle() const { return _tpart; }
    auto const& hits() const { return _hits;}
    auto const& hitCalibInfos() const { return _hitcalibs;}
    auto const& caloHit() const { return _chit; }
    auto const& straws() const { return _straws;}
    auto const& segments() const { return _segments; }
    auto const& intersections() const { return _inters; }
    auto const& status() const { return _status; }
    auto const& domainBounds() const { return _domainbounds; }
    double t0Val() const;
    double t0Var() const;
    float chisquared() const { return _chisq; }
    int nDOF() const { return _ndof; }
    float domainTolerance() const { return _dtol; }
    unsigned nHits(bool active=true) const;
    float fitConsistency() const { return _fitcon; }
    UInt_t nTrajSegments() const { return _segments.size(); }
    KinKal::TimeRange timeRange() const { return KinKal::TimeRange(_segments.front()._tmin,_segments.back()._tmax); }
    bool hasCaloCluster() const { return _chit.caloCluster().isNonnull(); }
    art::Ptr<CaloCluster> const& caloCluster() const { return _chit.caloCluster(); }
    std::vector<KalSegment>::const_iterator nearestSegment(double time)  const;
    std::vector<KalSegment>::const_iterator t0Segment(double& t0) const; // return the segment associated with the t0 value, or the closest to it.  Also return the t0 value
    InterIterCol intersections(SurfaceId const& id)  const; // all intersections with this SurfaceId
    bool loopHelixFit() const { return _status.hasAllProperties(TrkFitFlag::KKLoopHelix); }
    bool centralHelixFit() const { return _status.hasAllProperties(TrkFitFlag::KKCentralHelix); }
    bool kinematicLineFit() const { return _status.hasAllProperties(TrkFitFlag::KKLine); }
    // reconstitute (as best as possible) the fit trajectory.  The ptr will be null if the fit wasn't based on the requested trajector type
    // Note these return by value
    // Note that the returned piecetraj may have large gaps, unless the full fit trajectory was stored in the seed.
    LHPTPtr loopHelixFitTrajectory() const;
    CHPTPtr centralHelixFitTrajectory() const;
    KLPTPtr kinematicLineFitTrajectory() const;

    // global information about the track
    PDGCode::type   _tpart = PDGCode::unknown; // particle assumed for this fit
    TrkFitFlag      _status; // status of this fit: includes alglorithm information
    float           _chisq = -1; // fit chisquared value
    int             _ndof = -1;
    float           _fitcon = -1; // fit consistency
    float           _maxgap = 0;
    float           _avggap = 0; // information about trajectory gaps
    float           _dtol = 0; // tolerance used when creating BField domain
    //
    // contained content substructure.
    //
    std::vector<KalSegment>       _segments; // segments of the Kalman filter fit result
    std::vector<KalIntersection>  _inters; // intersections with materials or reference locations
    std::vector<TrkStrawHitSeed>  _hits; // hit seeds for all the hits used in this fit
    std::vector<TrkStrawHitCalib> _hitcalibs; // extra calibration/alignment info
    std::vector<TrkStraw>         _straws; // straws interesected by this fit
    std::vector<double>           _domainbounds; // domain time boundaries
    TrkCaloHitSeed                _chit;  // CaloCluster-based hit.  If it has no CaloCluster, this has no content
    //
    // deprecated BTrk legacy content, DO NOT write any new code which depends on these functions
    // find the nearest segment to a given the time
    std::vector<KalSegment>::const_iterator nearestSegmentFlt(float fltlen)  const;
    std::vector<KalSegment>::const_iterator nearestSegment(const XYZVectorF& pos)  const; // find nearest segment to a GLOBAL position
    float flt0() const { return _flt0; }
    HitT0 t0() const;
    float _flt0 = 0.0; // flight distance where the track crosses the tracker midplane (z=0).  Redundant with t0 in KinKal fits, and in the wrong unit
  };
  typedef std::vector<mu2e::KalSeed> KalSeedCollection;
  typedef art::Ptr<mu2e::KalSeed> KalSeedPtr;
  typedef std::vector<mu2e::KalSeedPtr> KalSeedPtrCollection;
  typedef std::vector<mu2e::KalSeedPtr> KalSeedCluster;
  typedef std::vector<mu2e::KalSeedCluster> KalSeedClusterCollection;
}
#endif
