//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
//
// $Id: VolumeInfo.cc,v 1.2 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  CLHEP::Hep3Vector VolumeInfo::_Mu2eOriginInWorld;


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
