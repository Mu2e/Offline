#ifndef TrackerGeom_Plane_hh
#define TrackerGeom_Plane_hh
//
// Holds information about one plane in a tracker.
//

//
// $Id: Plane.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <vector>
#include <iostream>

#include "TTrackerGeom/inc/PlaneId.hh"
#include "TTrackerGeom/inc/Face.hh"

namespace mu2e {

  class TTracker;

  class Plane{

    friend class TTracker;
    friend class TTrackerMaker;

    // Note that _z needs to be set by the TTrackerMaker - thus this friend

  public:

    Plane() : _id(PlaneId(-1,-1)), _z(0.0)               {}
    Plane( const PlaneId& id, double z ): _id(id), _z(z) {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const PlaneId& id() const { return _id;}

    // The next layer down in the hierarchy is Faces
    int   nFaces()                       const { return _faces.size(); }
    const std::vector<Face>& getFaces () const { return _faces; }
    const Face& getFace ( int n)         const { return _faces.at(n); }
    const Face& getFace ( const FaceId& fid ) const{
      return _faces.at(fid.getFace());
    }

    // Unlike the Device/Panel/Layer abstraction, one does not go directly
    // from a Plane to a Panel or a Layer or Straw

    // Get geometric abtraction information about this Plane
    double midZ() const {return _z;}

    // Formatted string embedding the id of the plane.
    std::string name( std::string const& base ) const;


    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& Ttracker ) const;

  protected:

    PlaneId _id;
    double  _z;
    
    // TODO MAYBE - if we need to allow for rotation variances, we may have 
    //              to add that here.  At this time, we don't know how we will 
    //              handle adjustments off the basic geometry, so we keep this 
    //              simple.

    std::vector<Face> _faces;

 };

 inline
 std::ostream & operator<< (std::ostream & os, Plane const &x) 
    { return os << x.name("Plane "); }

}  //namespace mu2e
#endif /* TrackerGeom_Plane_hh */
