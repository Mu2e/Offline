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
#include <array>

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class Tracker;

  class Plane{
    using PanelCollection = std::array<const Panel*,StrawId::_npanels>;
    using TrackerPanelCollection = std::array<Panel,StrawId::_nupanels>;
    using xyzVec = CLHEP::Hep3Vector;

    // Tracker re-points the panel pointers when it is copied
    friend class Tracker;

    public:

    Plane(){} // default object is not functional but needed for storage classes
    // construct from the full set of panels
    explicit Plane( const StrawId& id, TrackerPanelCollection const& panels);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const StrawId&  id()  const { return _id;}

    const xyzVec& origin() const { return _PlanetoDS.displacement(); }
    xyzVec uDirection() const { return _udir; }
    xyzVec vDirection() const { return _vdir; }
    xyzVec wDirection() const { return _wdir; }

    PanelCollection const& panels() const { return _panels; }

    size_t nPanels() const{
      return _panels.size();
    }

    const PanelCollection& getPanels () const{
      return _panels;
    }

    const Panel& getPanel ( int n) const {
      return *_panels.at(n);
    }

    const Panel& getPanel ( const StrawId& pnlid ) const{
      return *_panels.at(pnlid.getPanel());
    }

    const Straw& getStraw ( const StrawId& strid ) const{
      return getPanel(strid).getStraw(strid);
    }

    // Formatted string embedding the id of the panel.
    std::string name( std::string const& base ) const;

    auto const& planeToDS() const { return _PlanetoDS; }
    auto dsToPlane() const { return _PlanetoDS.inverse(); }

    void setPlaneToDS(HepTransform const& planeToDS) {
      _PlanetoDS = planeToDS;
      _udir = _PlanetoDS.rotation()*xyzVec(1.0,0.0,0.0);
      _vdir = _PlanetoDS.rotation()*xyzVec(0.0,1.0,0.0);
      _wdir = _PlanetoDS.rotation()*xyzVec(0.0,0.0,1.0);
    }

    private:
    // Fill the panel pointers from the given collection.  This is used both when
    // building a plane and when a Tracker is copied, since the pointers must then
    // refer to the panels owned by the copy and not to those of the original.
    void rebindPanels( TrackerPanelCollection const& panels );

    StrawId             _id;
    HepTransform        _PlanetoDS; // transform from plane coordinates to DS (just a translation)
    xyzVec _udir, _vdir, _wdir; // direction vectors in DS frame
    PanelCollection     _panels;
    static StrawIdMask  _sidmask; // mask to plane level
  };

} //namespace mu2e

#endif /* TrackerGeom_Plane_hh */
