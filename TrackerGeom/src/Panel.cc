//
// Hold information about one Panel in a tracker.
//
//
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Panel.hh"

using namespace std;

namespace mu2e {

  string Panel::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getPlane() << "_"
       << _id.getPanel();

    return os.str();
  }

  // void Panel::fillPointers ( const Tracker& tracker ) const {
  //   for( size_t i=0; i<_layers.size(); ++i ){
  //     _layers[i].fillPointers(tracker);
  //     _straw0MidPoint += _layers[i].straw0MidPoint();
  //   }
  //   _straw0Direction = _layers[0].straw0Direction();
  //   _straw0MidPoint  /= _layers.size();
  // }

  void Panel::fillPointers ( const Tracker* tracker ) const {

    for(auto& isp : getStrawPointers()) {
      isp->fillPointers(tracker);
    }


    _straw0MidPoint = CLHEP::Hep3Vector();
    for ( uint16_t ilay=0; ilay<StrawId::_nlayers; ++ilay ) {
      const Straw& straw = getStraw( StrawId(id().asUint16() + ilay ) );
      _straw0MidPoint += straw.getMidPoint();
    }
    _straw0Direction = getStraw(0).getDirection();
    _straw0MidPoint  /= StrawId::_nlayers;
  }

} // end namespace mu2e
