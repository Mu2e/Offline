//
// Hold information about one Sector in a tracker.
//
//
// $Id: Sector.cc,v 1.2 2010/04/18 00:31:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:31:56 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Sector.hh"

using namespace std;

namespace mu2e {

  string Sector::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getDevice() << "_"
       << _id.getSector();

    return os.str();
  }

  void Sector::fillPointers ( const Tracker& tracker ) const {
    for( size_t i=0; i<_layers.size(); ++i ){
      _layers[i].fillPointers(tracker);
    }
  }

} // end namespace mu2e
