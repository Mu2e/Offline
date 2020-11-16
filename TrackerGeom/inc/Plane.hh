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
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "TrackerGeom/inc/Panel.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class Plane{
    using PanelCollection = std::array<const Panel*,StrawId::_npanels>;
    using TrackerPanelCollection = std::array<Panel,StrawId::_nupanels>; 
    using xyzVec = CLHEP::Hep3Vector;

    public:

    Plane(){} // default object is not functional but needed for storage classes
    // construct from the full set of panels
    explicit Plane( const StrawId& id, TrackerPanelCollection const& panels);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const StrawId&  id()  const { return _id;}

    const xyzVec& origin() const { return _origin; }
    xyzVec& origin() { return _origin; }

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

    private:
    StrawId             _id;
    xyzVec              _origin; // deprecated
    PanelCollection     _panels;
    static StrawIdMask  _sidmask; // to define equivalence
  };

} //namespace mu2e

#endif /* TrackerGeom_Plane_hh */
