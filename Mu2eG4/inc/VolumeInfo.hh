#ifndef VOLUMEINFO_HH
#define VOLUMEINFO_HH
//
// Return type for the nestBox and nestTubs free functions.
// Just a struct with a default c'tor.
// 
// $Id: VolumeInfo.hh,v 1.2 2010/02/05 11:46:00 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2010/02/05 11:46:00 $
//
// Original author Rob Kutschke
//

class G4VSolid;
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

    G4VSolid*           solid;
    G4LogicalVolume*    logical;
    G4VPhysicalVolume* physical;
  };

}

#endif
