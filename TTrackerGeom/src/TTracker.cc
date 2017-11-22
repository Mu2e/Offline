//
// Geometry and identifier info about an TTracker.
//
//
// $Id: TTracker.cc,v 1.8 2013/01/07 04:01:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/07 04:01:16 $
//
// Original author Rob Kutschke
//

#include "TTrackerGeom/inc/TTracker.hh"

using namespace std;

namespace mu2e {

  void TTracker::fillPointers () const{
    for ( size_t i=0; i<_planes.size(); ++i){
      _planes[i].fillPointers(*this);
    }
  }

  // void TTracker::fillPointers2 () const{
  //   for ( size_t i=0; i<_planes.size(); ++i){
  //     _planes[i].fillPointers(*this);
  //   }
  // }

} // namespace mu2e
