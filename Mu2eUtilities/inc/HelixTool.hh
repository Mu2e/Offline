#ifndef Mu2eUtilities_HelixTool_hh
#define Mu2eUtilities_HelixTool_hh
//
// Original author G. Pezzullo
//

#include <functional>

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

namespace mu2e { class Tracker; }
namespace mu2e { class Plane; }
namespace mu2e { struct HelixSeed; }

using namespace ROOT::Math::VectorUtil;

namespace mu2e {

  class HelixTool{

  public:
    // returns true for a straw that can't produce a hit (e.g. TrackerStatus noSignal or suppress).
    // An empty function means all straws are live. A function is used instead of TrackerStatus
    // to avoid a library dependency of Mu2eUtilities on TrackerConditions.
    using StrawNoSignal = std::function<bool(StrawId const&)>;

    HelixTool(const HelixSeed *Helix, const mu2e::Tracker*MyTracker, StrawNoSignal NoSignal = StrawNoSignal());

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    float  nLoops           () const { return _nLoops;            }
    float  meanHitRadialDist() const { return _meanHitRadialDist; }
    float  d0               () const { return _d0;                }
    int    nstrawhits       () const { return _nStrawHits;        }

    //function that evaluates the ratio between the measured straw hits and the
    //expected intersections of the helix with the tracker plane faces inside the tracker
    //active radius. If a StrawNoSignal function is given, each face crossing is weighted by
    //the fraction of live straws (nearest straw in each layer) at the crossing point, so that
    //dead panels and straws don't lower the ratio. A crossing with no straw of the face at that
    //point (outside the panel's straws) keeps the weight 1. Evaluated on the first call.
    float  hitRatio         () const;

    void   dirOfProp(float& slope, float& slopeErr, float& chi2ndof);

  private:
    float  liveFraction(const Plane& Pln, uint16_t Face, float X, float Y) const;

    const HelixSeed* _hel;

    float      _nLoops;
    int        _nStrawHits;
    float      _meanHitRadialDist;
    float      _d0;
    mutable float _hitRatio = -1.; // < 0: not evaluated yet
    float      _trackerRIn, _trackerROut;
    const mu2e::Tracker*_tracker;
    StrawNoSignal _noSignal;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_HelixTool_hh */
