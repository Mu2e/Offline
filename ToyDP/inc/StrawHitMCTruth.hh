#ifndef ToyDP_StrawHitMCTruth_hh
#define ToyDP_StrawHitMCTruth_hh
//
// This is a place to put additional information produced by HitMaker,
// like true drift distance, signal propagation time, etc.
//
// Original author Ivan Logashenko
//

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  struct StrawHitMCTruth{

  public:

    StrawHitMCTruth():
      _t0(0.),
      _driftTime(0.),
      _driftDistance(0.),
      _distanceToMid(0.) {
    }

    // Constructor for a hit that came from an unpacked digi, either
    // from data or from the full MC chain.
    StrawHitMCTruth(float t0,
		    float driftTime,
		    float driftDistance,
		    float distanceToMid) :
      _t0(t0),
      _driftTime(driftTime),
      _driftDistance(driftDistance),
      _distanceToMid(distanceToMid) {
    }

    // Accessors
    float   t0()             const { return _t0;}
    float   driftTime()      const { return _driftTime;}
    float   driftDistance()  const { return _driftDistance;}
    float   distanceToMid()  const { return _distanceToMid; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    float   _t0;            // ns
    float   _driftTime;     // ns
    float   _driftDistance; // mm
    float   _distanceToMid; // mm

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* ToyDP_StrawHitMCTruth_hh */
