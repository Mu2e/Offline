#ifndef TrackerGeom_Panel_hh
#define TrackerGeom_Panel_hh
//
// Holds information about one panel in a tracker.
//
// Original author Rob Kutschke
//

#include <array>

#include "cetlib_except/exception.h"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"

namespace mu2e {

  class Panel{
    using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

    public:
    using StrawCollection = std::array<const Straw*, StrawId::_nstraws>;
    using TrackerStrawCollection = std::array<Straw,StrawId::_nustraws>;
    Panel():_straws(){} // default object non-function but needed for storage classes
    Panel( const StrawId& id, TrackerStrawCollection const& straws ); // construct from the Id and the full set of straws

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const StrawId& id() const { return _id;}
    const auto& getStrawPointers() const{
      return _straws;
    }

    StrawCollection const& straws() const { return _straws; }

    size_t nStraws() const { return _straws.size(); }

    const Straw& getStraw( uint16_t n ) const {
      return *(_straws.at(n));
    }

    // this getStraw checks if this is the right panel
    const Straw& getStraw ( const StrawId& strid ) const{
      if ( _sidmask.equal(_id,strid) ) {
        return *(_straws.at(strid.straw()));
      } else {
        std::ostringstream msg;
        msg << static_cast<const char*>(__func__) << " Inconsistent straw/panel request "
          << strid << " / " << _id << std::endl;
        throw cet::exception("RANGE") << msg.str();
      }
    }
    xyzVec origin()  const { return _UVWtoDS.displacement(); }
    // local (UVW) coordinate system.
    xyzVec uDirection() const { return _udir; }
    xyzVec vDirection() const { return _vdir; }
    xyzVec wDirection() const { return _wdir; }

    // transform from local to DS coordinates
    auto const& panelToDS() const { return _UVWtoDS; }
    auto dsToPanel() const { return _UVWtoDS.inverse(); }

    // deprecated interface: either use the above local coordinates, or get the straw direction directly
    xyzVec straw0Direction() const { return _straws[0]->wireDirection(); }
    // (The primary straw of each layer is the straw used to establish position.
    //  In the Tracker the primary straw is the innermost straw.)
    // *** In a multi-layer geometry, the straw0MidPoint ***
    // ***        need not lie on any actaul straw       ***
    // this function is misnamed and deprecated, users should use origin instead  FIXME!
    xyzVec straw0MidPoint()  const { return origin(); }
    // deprecated useless function
    size_t nLayers() const { return 2; }

    // Formatted string embedding the id of the panel.
    std::string name( std::string const& base ) const;

    private:
    StrawId _id; // only the plane and panel fields are used to define a panel
    xyzVec _udir, _vdir, _wdir; // direction vectors in DS frame
    HepTransform  _UVWtoDS; // transform from this panel's frame to the DS frame
    // indirection to straws
    StrawCollection _straws;
    static StrawIdMask _sidmask; // mask to panel level
  };

}  //namespace mu2e
#endif /* TrackerGeom_Panel_hh */
