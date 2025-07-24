//
// Hold information about one Panel in a tracker.
//
// Original author Rob Kutschke
//

#include "Offline/TrackerGeom/inc/Panel.hh"
#include "cetlib_except/exception.h"
#include <sstream>

using namespace std;
using CLHEP::HepRotation;
using CLHEP::Hep3Vector;

namespace mu2e {
  // define equivalence unique panels
  StrawIdMask Panel::_sidmask(StrawIdMask::uniquepanel);

  // non-functional constructor needed by Proditions system: this should be eliminated FIXME!
  string Panel::name( string const& base ) const{
    ostringstream os;
    os << base
      << _id.getPlane() << "_"
      << _id.getPanel();

    return os.str();
  }

  Panel::Panel( const StrawId& id, TrackerStrawCollection const& straws ) : _id(id) {
    for(auto const& straw : straws ) {
      // pick out all the straws belonging to this panel.
      if(_sidmask.equal(_id,straw.id())){
        _straws[straw.id().straw()] = &straw;
      }
    }
    // compute the panel coordinate axes based on the straw content
    // U points along the straw (Cal to HV), V is radially outward, W is given by right-handedness
    _udir = _straws.front()->wireDirection();
    _vdir = (_straws[StrawId::_nstraws-2]->origin()-_straws[0]->origin()).unit(); // must pick straws in the same layer!
    _wdir = _udir.cross(_vdir);
    // Define the Transform from UVW to the DS frame.  We pick the panel 'middle' to mininimize non-linear rotation effects
    auto origin = 0.5*(_straws[47]->origin() + _straws[48]->origin());
    // nominal rotation is about Z.  We must also account for the flips
    HepRotation prot(0.0,0.0,-_udir.phi());
    HepRotation flip(Hep3Vector(1.0,0.0,0.0),M_PI);
    if(_wdir.z() < 0.0)prot *= flip;
    _UVWtoDS = HepTransform(origin, prot);
    // tests
    auto DStoUVW = _UVWtoDS.inverse();
    auto udir = DStoUVW.rotation()*_udir;
    auto vdir = DStoUVW.rotation()*_vdir;
    auto wdir = DStoUVW.rotation()*_wdir;
//    if( fabs(1.0 - udir.dot(Hep3Vector(1.0,0.0,0.0))) > 1e-6 ||
//        fabs(1.0 - vdir.dot(Hep3Vector(0.0,1.0,0.0))) > 1e-6 ||
//        fabs(1.0 - wdir.dot(Hep3Vector(0.0,0.0,1.0))) > 1e-6 )
//      throw cet::exception("Geom") << "Panel direction error: id " << _id << " udir " << udir << " vdir " << vdir << " wdir " << wdir << std::endl;
//    auto po = DStoUVW*origin;
//    if( fabs(po.r()) > 1e-6 )
//      throw cet::exception("Geom") << "Panel origin error: id " << _id << " po " << po << " transform " << DStoUVW << std::endl;
    //
  }

} // end namespace mu2e
