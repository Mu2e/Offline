#ifndef TTrackerGeom_PlaneMF_hh
#define TTrackerGeom_PlaneMF_hh
//
// Holds information about one planeMF in a tracker.
//

//
// $Id: PlaneMF.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <vector>
#include <iostream>

#include "TTrackerGeom/inc/PlaneMFId.hh"
#include "TTrackerGeom/inc/Face.hh"

namespace mu2e {

  class TTracker;

  class PlaneMF{

    friend class TTracker;
    friend class TTrackerMaker;

    // Note that _z needs to be set by the TTrackerMaker - thus this friend

  public:

    PlaneMF() : _id(PlaneMFId(-1,-1)), _z(0.0)               {}
    PlaneMF( const PlaneMFId& id, double z ): _id(id), _z(z) {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const PlaneMFId& id() const { return _id;}

    // The next layer down in the hierarchy is Faces
    int   nFaces()                       const { return _faces.size(); }
    const std::vector<Face>& getFaces () const { return _faces; }
    const Face& getFace ( int n)         const { return _faces.at(n); }
    const Face& getFace ( const FaceId& fid ) const{
      return _faces.at(fid.getFace());
    }

    // Unlike the Plane/Panel/Layer abstraction, one does not go directly
    // from a PlaneMF to a Panel or a Layer or Straw

    // Get geometric abtraction information about this PlaneMF
    double midZ() const {return _z;}

    // Formatted string embedding the id of the planeMF.
    std::string name( std::string const& base ) const;


    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& Ttracker ) const;

  protected:

    PlaneMFId _id;
    double  _z;
    
    // TODO MAYBE - if we need to allow for rotation variances, we may have 
    //              to add that here.  At this time, we don't know how we will 
    //              handle adjustments off the basic geometry, so we keep this 
    //              simple.

    std::vector<Face> _faces;

 };

 inline
 std::ostream & operator<< (std::ostream & os, PlaneMF const &x) 
    { return os << x.name("PlaneMF "); }

}  //namespace mu2e
#endif /* TTrackerGeom_PlaneMF_hh */
