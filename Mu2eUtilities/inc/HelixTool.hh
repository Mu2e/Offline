#ifndef Mu2eUtilities_HelixTool_hh
#define Mu2eUtilities_HelixTool_hh
#include "Math/GenVector/VectorUtil.h"  // for VectorUtil

namespace mu2e {
class Tracker;
struct HelixSeed;
}  // namespace mu2e
//
// Original author G. Pezzullo
//


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
