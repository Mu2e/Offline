//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
//
// $Id: VolumeInfo.cc,v 1.3 2011/11/03 16:28:44 gandr Exp $
// $Author: gandr $
// $Date: 2011/11/03 16:28:44 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

namespace mu2e {

  VolumeInfo::VolumeInfo( const std::string&  pName,
                          const CLHEP::Hep3Vector& inParent,
                          const CLHEP::Hep3Vector& parentInWorld):
    name(pName),
    solid(0),
    logical(0),
    physical(0),
    centerInParent(inParent){
    centerInWorld = centerInParent + parentInWorld;
  }

  const CLHEP::Hep3Vector& VolumeInfo::mu2eOriginInWorld() {
    // AG: Hiding this static in the function ensures that its
    // initialization happes on the first access (after geometry svc
    // is available) and not during the library loading, which would
    // be too early.
    static CLHEP::Hep3Vector _origin(GeomHandle<WorldG4>()->mu2eOriginInWorld());
    return _origin;
  }

} // end namespace mu2e
