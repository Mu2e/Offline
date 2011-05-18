#ifndef Mu2eG4_VolMapType_hh
#define Mu2eG4_VolMapType_hh
//
// In the run data there is a data product that describes all physical
// volumes the run of G4.
//
// Given a pointer to a physical volume, return the index into the data product
// for that volume.
//
// $Id: VolMapType.hh,v 1.4 2011/05/18 15:47:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:47:40 $
//
// Original author Rob Kutschke
//

#include <map>

class G4VPhysicalVolume;

namespace mu2e{

  typedef std::map<G4VPhysicalVolume*,unsigned> VolMapType;
  typedef std::map<G4VPhysicalVolume*,unsigned>::iterator VolMapType_iterator;
  typedef std::map<G4VPhysicalVolume*,unsigned>::const_iterator VolMapType_const_iterator;

}

#endif /* Mu2eG4_VolMapType_hh */
