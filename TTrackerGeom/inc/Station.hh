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
#include "TTrackerGeom/inc/Plane.hh"
#include "TTrackerGeom/inc/Face.hh"

// CLHEP includes

namespace mu2e {

  class TTracker;

  class Station{

    friend class TTracker;
    friend class TTrackerMaker;

  public:

    Station():_id(-1){}

    explicit Station( const StationId& id, double z = 0. )
      : _id(id)
      , _z(z)
    {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const StationId& id() const { return _id;}

    // TODO MAYBE - If we allow for variances in rotatoinal orientation
    //              and/or XY origin, we ought to have accessors for that
    //              information here.
    
    // The next layer down in the hierarchy is Planes
    int   nPlanes()                        const { return _planes.size(); }
    const std::vector<Plane>& getPlanes () const { return _planes; }
    const Plane& getPlane ( int n)         const { return _planes.at(n); }
    const Plane& getPlane ( const PlaneId& pid ) const{
      // TODO -- throw if pid.getStation() does not match _id
      return _planes.at(pid.getPlane());
    }

    // Often, algorithms skip the Planes layer and work with Faces
    int nFaces()                                 const { return _faces.size(); }
    const std::vector<Face const *>& getFaces () const { return _faces; }
    const Face& getFace ( int n)               const { return *(_faces.at(n)); }
    const Face& getFace ( const FaceId& pid )  const{
      // TODO -- throw if pid.getStation() does not match _id
      return *(_faces.at(pid.getFace()));
    }

    // Unlike the Device/Panel/Layer abstraction, one does not go directly
    // from a Station to a Panel or a Layer

    // Get geometric abtraction information about this Station
    double midZ() const {return _z;}

    // Formatted string embedding the id of the station.
    std::string name( std::string const& base ) const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& tracker ) const;

  protected:

    StationId          _id;
    double             _z;  // This is Z of the CENTER of the station.
 
    // TODO MAYBE - if we need to allow for XY variances, we will have 
    //              to replace _z with a Hep3Vector _origin.
    
    // TODO MAYBE - if we need to allow for rotation variances, we will have 
    //              to add that here.  At this time, we don't know how we will 
    //              handle adjustments off the basic geometry, so we keep this 
    //              simple.  (Similarly for XY origin variances.)
    
    std::vector<Plane>        _planes;
    std::vector<Face const *> _faces;

  };

  inline
  std::ostream & operator<< (std::ostream & os, Station const &x) 
    { return os << x.name("Station "); }

} //namespace mu2e

#endif /* TTrackerGeom_Station_hh */
