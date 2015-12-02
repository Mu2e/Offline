#ifndef TrackerGeom_PanelMF_hh
#define TrackerGeom_PanelMF_hh
//
// Holds information about one face in a tracker.
//

//
// $Id: PanelMF.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <vector>
#include <iostream>

#include "TTrackerGeom/inc/ZLayer.hh"
#include "TTrackerGeom/inc/View.hh"

namespace mu2e {

  class TTracker;

  class PanelMF{

   friend class TTracker;
   friend class TTrackerMaker;

    // Note that _phi and _view need to be set by the 
    // TTrackerMaker - thus this friend

  public:

    PanelMF():_id(PanelMFId()) {}
    PanelMF( const PanelMFId& id, double z ) : _id(id), _z(z) {}

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const PanelMFId& id() const { return _id;}

    // The next layer down in the hierarchy is ZLayers
    int   nZLayers()                         const { return _zlayers.size(); }
    const std::vector<ZLayer>& getZLayers () const { return _zlayers; }
    const ZLayer& getZLayer ( int n)         const { return _zlayers.at(n); }
    const ZLayer& getZLayer ( const ZLayerId& zid ) const{
      return _zlayers.at(zid.getZLayer());
    }

    // Unlike the Device/Panel/Layer abstraction, one does not go directly
    // from a PanelMF to a Straw

   // Get geometric abtraction information about this panelMF
   
   double midZ() const { return _z; } 
   View view()   const { return _view; }
   double phi()  const { return _phi; }
   


    // Formatted string embedding the id of the panelMF.
    std::string name( std::string const& base ) const;

   // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const TTracker& Ttracker ) const;

  protected:

    PanelMFId _id;
    double  _z;
    View    _view;
    double  _phi;
    
    std::vector<ZLayer> _zlayers;

 };

  inline
  std::ostream & operator<< (std::ostream & os, PanelMF const &x) 
    { return os << x.name("PanelMF "); }

}  //namespace mu2e
#endif /* TrackerGeom_PanelMF_hh */
