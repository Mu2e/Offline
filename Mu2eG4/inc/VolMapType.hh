#ifndef VolMapType_H
#define VolMapType_H
//
// In the run data there is a data product that describes all physical
// volumes the run of G4.  
//
// Given a pointer to a physical volume, return the index into the data product
// for that volume.
//
// $Id: VolMapType.hh,v 1.1 2010/03/23 20:34:30 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/23 20:34:30 $
//
// Original author Rob Kutschke
//

#include <map>

class G4VPhysicalVolume;

namespace mu2e{

  typedef std::map<G4VPhysicalVolume*,uint32_t> VolMapType;
  typedef std::map<G4VPhysicalVolume*,uint32_t>::iterator VolMapType_iterator;
  typedef std::map<G4VPhysicalVolume*,uint32_t>::const_iterator VolMapType_const_iterator;

}

#endif
