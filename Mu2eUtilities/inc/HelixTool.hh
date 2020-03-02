#ifndef Mu2eUtilities_HelixTool_hh
#define Mu2eUtilities_HelixTool_hh
//
// Original author G. Pezzullo
//

#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "DataProducts/inc/StrawId.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;

namespace mu2e {

  class HelixTool{

  public:
    HelixTool(const HelixSeed *Helix, const mu2e::Tracker*MyTracker);

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    float  nLoops           () const { return _nLoops;            }
    int    nMinHitsLoop     () const { return _nMinHitsLoop;      }
    float  meanHitRadialDist() const { return _meanHitRadialDist; }
    float  d0               () const { return _d0;                }
    float  nstrawhits       () const { return _nStrawHits;        }

    //function that evaluates the ratio between the measured ComboHits and the
    //expected intesections of the helix with the tracker planes. This function
    //models the tracker as a perfect cylinder
    float  hitRatio         () const { return _hitRatio;          }


  private:
    const HelixSeed* _hel;
    int        _nMinHitsLoop;

    float      _nLoops;
    int        _nStrawHits;
    float      _meanHitRadialDist;
    float      _d0;
    float      _hitRatio;
    float      _trackerRIn, _trackerROut, _trackerLength;
    const mu2e::Tracker*_tracker;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_HelixTool_hh */
