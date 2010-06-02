//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
// 
// $Id: VolumeInfo.cc,v 1.1 2010/06/02 04:00:50 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/06/02 04:00:50 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/VolumeInfo.hh"

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

} // end namespace mu2e
