#ifndef TTrackerGeom_Station_hh
#define TTrackerGeom_Station_hh
//
// Hold information about one station in a tracker.
//
// $Id: Station.hh,v 1.2 2012/01/25 20:29:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/01/25 20:29:23 $
//
// Original author Mark Fischler
//

// C++ includes
#include <vector>
#include <iostream>

// Mu2e includes
#include "TTrackerGeom/inc/StationId.hh"
#include "TrackerGeom/inc/Plane.hh"

// CLHEP includes

namespace mu2e {

  class TTracker;

  class Station{

    friend class TTracker;
    friend class TTrackerMaker;

  public:

    Station():_id(-1){}

    // The next layer down in the hierarchy is Planes
    int nPlanes() const { return _planes.size(); }

    const std::vector<Plane>& getPlanes()  const { return _planes; }
    const Plane& getPlane ( int n)         const { return _planes.at(n); }

    explicit Station( const StationId& id, double z = 0. )
      : _id(id)
      , _z(z)
    {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const StationId& id() const { return _id;}

    // Get geometric abtraction information about this Station
    double midZ() const {return _z;}

    // Formatted string embedding the id of the station.
    std::string name( std::string const& base ) const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& tracker ) const;

  protected:

    StationId          _id;
    double             _z;  // This is Z of the CENTER of the station.
    std::vector<Plane>   _planes;

  };

  inline
  std::ostream & operator<< (std::ostream & os, Station const &x) 
    { return os << x.name("Station "); }

} //namespace mu2e

#endif /* TTrackerGeom_Station_hh */
