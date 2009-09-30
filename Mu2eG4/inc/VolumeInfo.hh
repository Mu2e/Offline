#ifndef VOLUMEINFO_HH
#define VOLUMEINFO_HH
//
// Return type for the nestBox and nestTubs free functions.
// Just a struct with a default c'tor.
// 
// $Id: VolumeInfo.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

class G4CSGSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


namespace mu2e {

  class VolumeInfo{

  public:

    VolumeInfo():
      solid(0),
      logical(0),
      physical(0){}
    ~VolumeInfo(){}

    G4CSGSolid*           solid;
    G4LogicalVolume*    logical;
    G4VPhysicalVolume* physical;
  };

}

#endif
