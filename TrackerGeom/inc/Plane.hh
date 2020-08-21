#ifndef TrackerGeom_Plane_hh
#define TrackerGeom_Plane_hh
//
// Hold information about one plane in a tracker.
//
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>

// Mu2e includes
#include "DataProducts/inc/PlaneId.hh"
#include "TrackerGeom/inc/Panel.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Plane{

    friend class AlignedTrackerMaker;

    friend class TrackerMaker;
    friend class Tracker; // needed for deep copy

  public:

    // A free function, returning void, that takes a const Plane& as an argument.
    typedef void (*PlaneFunction)( const Plane& s);

    Plane():_id(PlaneId()),_origin(),_panels(),_exists(true){}

    explicit Plane( const PlaneId& id,
            CLHEP::Hep3Vector const& origin = CLHEP::Hep3Vector(0.,0.,0.),
            bool exists = true ):
      _id(id),
      _origin(origin),
      _exists(exists) {
    }

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const PlaneId&  id()  const { return _id;}

    const CLHEP::Hep3Vector & origin() const { return _origin; }

    int nPanels() const{
      return _panels.size();
    }

    const std::array<Panel const*,StrawId::_npanels>& getPanels () const{
      return _panels;
    }

    const Panel& getPanel ( int n) const {
      return *_panels.at(n);
    }

    const Panel& getPanel ( const PanelId& pnlid ) const{
      return *_panels.at(pnlid.getPanel());
    }

    const Straw& getStraw ( const StrawId& strid ) const{
      return _panels.at(strid.getPanel())->getStraw(strid);
    }

    // Formatted string embedding the id of the panel.
    std::string name( std::string const& base ) const;

    // On readback from persistency, recursively recompute mutable members.
    void fillPointers ( const Tracker* tracker ) const;

    bool exists() const {
      return _exists;
    }


  protected:

    PlaneId             _id;
    CLHEP::Hep3Vector   _origin;
    std::array<Panel const*,StrawId::_npanels> _panels;
    bool                _exists;
  };

} //namespace mu2e

#endif /* TrackerGeom_Plane_hh */
