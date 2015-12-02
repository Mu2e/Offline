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
#include "TTrackerGeom/inc/PanelMF.hh"

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

    // The next layer down in the hierarchy is PanelMFs
    int   nPanelMFs()                        const { return _panelMFs.size(); }
    const std::vector<PanelMF>& getPanelMFs () const { return _panelMFs; }
    const PanelMF& getPanelMF ( int n)         const { return _panelMFs.at(n); }
    const PanelMF& getPanelMF ( const PanelMFId& panid ) const{
      return _panelMFs.at(panid.getPanelMF());
    }

    // Unlike the Device/Panel/Layer abstraction, one does not go directly
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

    std::vector<PanelMF> _panelMFs;

 };

  inline
  std::ostream & operator<< (std::ostream & os, Face const &x) 
    { return os << x.name("Face "); }

}  //namespace mu2e
#endif /* TrackerGeom_Face_hh */
