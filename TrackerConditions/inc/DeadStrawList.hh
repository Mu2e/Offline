#ifndef HitMakers_DeadStrawList_hh
#define HitMakers_DeadStrawList_hh
//
// Maintain a list of dead straws.  This should be moved to the conditions system.
//
// $Id: DeadStrawList.hh,v 1.1 2014/05/30 19:15:32 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/30 19:15:32 $
//
// Original author Rob Kutschke
// Modified to model dead central section by D. Brown
//

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "TrackerGeom/inc/Straw.hh"

// art and its tool chain
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <set>
#include <iosfwd>

namespace mu2e {

// struct to define a range of deadness in a straw.  This is assumed to be equidistant
// from the straw center, to model radiation damage due to dose (which is concentrated
// at the lowest radius = straw center)
//
  struct DeadStrawRange {
    StrawIndex _strawind; // which straw
    double _range; // how much is dead
// constructor.  Default range is longer than the longest straw, meaning the
// straw is 100% dead
    DeadStrawRange(StrawIndex const& ind,double range=1000.0) : _strawind(ind), _range(range) {}
    DeadStrawRange(Straw const& straw,double range=1000.0) : DeadStrawRange(straw.index(),range) {}
// test if a given position (relative to the straw center) is dead
    bool isDead(double hitpos) const { return fabs(hitpos)<_range; }
// define comparator to ignore range: this allows finding an element based on strawindex quickly
    bool operator < (DeadStrawRange const& other) const { return _strawind < other._strawind; }
  };

  class DeadStrawList{

  public:
    DeadStrawList( fhicl::ParameterSet const& pset );
    // Accept compiler supplied d'tor, copy assignment and copy c'tor.

    // Reset the state when conditions information changes.
    void reset( fhicl::ParameterSet const& pset );

    // Accessors; test if a straw is dead.  The position argument refers to distance from the straw center in mm
    bool isAlive( StrawIndex ind,double hitpos=0.0 ) const { return !(isDead(ind,hitpos)); }
    bool isDead ( StrawIndex ind,double hitpos=0.0 ) const;

    void print( std::ostream& ) const;

  private:

    // Control printout
    int _verbosity;
    std::set<DeadStrawRange> _deadstraws; // sparse list of dead straws.  Ordering allows rapid lookup.

  };

} // namespace mu2e

#endif /* HitMakers_DeadStrawList_hh */
