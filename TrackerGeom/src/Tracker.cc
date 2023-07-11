//
// Geometry and identifier info about an Tracker.
//
//
//
// Original author Rob Kutschke
//

#include <utility>
#include "Offline/TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {

  Tracker::Tracker(StrawCollection const& straws, StrawProperties const& sprops,
      const TrackerG4InfoPtr& g4tracker, PEType const& pexists) :
    ProditionsEntity(cxname), _strawprops(sprops), _straws(straws),
    _planeExists(pexists), _g4tracker(g4tracker) {
      // build the panels from the straws
      for(uint16_t plane=0; plane < StrawId::_nplanes; plane++){
        for(uint16_t panel = 0;panel < StrawId::_npanels; panel++){
          StrawId sid(plane,panel,0);
          _panels[sid.uniquePanel()] = Panel(sid, _straws);
        }
      }
      // build the planes from the panels
      for(uint16_t plane=0; plane < StrawId::_nplanes; plane++){
        StrawId sid(plane,0,0);
        _planes[sid.plane()] = Plane(sid, _panels);
      }

    }

} // namespace mu2e
