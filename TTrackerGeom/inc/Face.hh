#ifndef TrackerGeom_Face_hh
#define TrackerGeom_Face_hh
//
// Holds information about one face in a tracker.
//

//
// $Id: Face.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <vector>
#include <iostream>

#include "TTrackerGeom/inc/FaceId.hh"
#include "TTrackerGeom/inc/Panel.hh"

namespace mu2e {

  class TTracker;

  class Face{

    friend class TTracker;
    friend class TTrackerMaker;

    // Note that _z needs to be set by the TTrackerMaker - thus this friend

  public:

    Face():_id(FaceId()), _z(0)                      {}
    Face( const FaceId& id, double z):_id(id), _z(z) {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const FaceId& id() const { return _id;}

    // The next layer down in the hierarchy is Panels
    int   nPanels()                        const { return _panels.size(); }
    const std::vector<Panel>& getPanels () const { return _panels; }
    const Panel& getPanel ( int n)         const { return _panels.at(n); }
    const Panel& getPanel ( const PanelId& panid ) const{
      return _panels.at(panid.getPanel());
    }

    // Unlike the Device/Sector/Layer abstraction, one does not go directly
    // from a Face to a ZLayer or Straw

    // Formatted string embedding the id of the face.
    std::string name( std::string const& base ) const;

   // Get the Z of the mid-z of this face
   
   double midZ() const { return _z; }
   
   // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& Ttracker ) const;

  protected:

    FaceId _id;
    double _z;

    // TODO MAYBE - if we need to allow for rotation variances, we may have 
    //              to add that here.  At this time, we don't know how we will 
    //              handle adjustments off the basic geometry, so we keep this 
    //              simple.

    std::vector<Panel> _panels;

 };

  inline
  std::ostream & operator<< (std::ostream & os, Face const &x) 
    { return os << x.name("Face "); }

}  //namespace mu2e
#endif /* TrackerGeom_Face_hh */
