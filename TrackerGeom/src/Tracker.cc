//
// Geometry and identifier info about an Tracker.
//
//
// $Id: Tracker.cc,v 1.8 2013/01/07 04:01:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/07 04:01:16 $
//
// Original author Rob Kutschke
//

#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {

  void Tracker::fillPointers () const{
    for ( size_t i=0; i<StrawId::_nplanes; ++i){
      _planes[i].fillPointers(this);
    }
  }

} // namespace mu2e
