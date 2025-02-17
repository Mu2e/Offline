//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
//
//
// Original author Rob Kutschke
//

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"

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

  VolumeInfo::VolumeInfo(const std::string&  pName):
    name(pName),
    solid(0),
    logical(0),
    physical(0),
    centerInParent(),
    centerInWorld() {
  }


  const CLHEP::Hep3Vector& VolumeInfo::mu2eOriginInWorld() {
    // AG: Hiding this static in the function ensures that its
    // initialization happes on the first access (after geometry svc
    // is available) and not during the library loading, which would
    // be too early.

    // WorldG4 may not be available in a "study" geometry, check before using.
    static art::ServiceHandle<GeometryService> sg;
    static CLHEP::Hep3Vector _origin(sg->hasElement<WorldG4>() ?
                                     GeomHandle<WorldG4>()->mu2eOriginInWorld()
                                     : CLHEP::Hep3Vector(0,0,0));
    return _origin;
  }

} // end namespace mu2e
