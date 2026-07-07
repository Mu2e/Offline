//
// Hold information about one Plane in a tracker.
//
//
// Original author Rob Kutschke
//

#include <sstream>

#include "Offline/TrackerGeom/inc/Plane.hh"

using namespace std;
using CLHEP::HepRotation;
using CLHEP::Hep3Vector;
namespace mu2e {
  // define equivalence plane
  StrawIdMask Plane::_sidmask(StrawIdMask::plane);

  string Plane::name( string const& base ) const{
    ostringstream os;

    os << base << _id;

    return os.str();
  }

  Plane::Plane( const StrawId& id, TrackerPanelCollection const& panels) : _id(id) {
    for(auto const& panel : panels ) {
      // pick out all the panels belonging to this plane.  This code relies on the Tracker collection being in order.
      if(_sidmask.equal(_id,panel.id())){
        _panels[panel.id().panel()] = &panel;
      }
    }
    // define the origin (in the tracker nominal frame) as the geometric average of the panel origins
    xyzVec origin;
    for(auto const& panel : _panels) origin += panel->origin();
    origin *= (1.0/6.0);
    // define the plane orientation based on panel 0; some planes are flipped, see doc 888 figure 9, table 4 for details
    auto wdir = _panels[0]->wDirection();
    HepRotation prot; // no default plane rotation
    if(wdir.z() < 0.0) // flip about Y
      prot = HepRotation(Hep3Vector(0.0,1.0,0.0),M_PI);
    _PlanetoDS = HepTransform(origin,prot);
    _udir = _PlanetoDS.rotation()*xyzVec(1.0,0.0,0.0);
    _vdir = _PlanetoDS.rotation()*xyzVec(0.0,1.0,0.0);
    _wdir = _PlanetoDS.rotation()*xyzVec(0.0,0.0,1.0);
  }

}
