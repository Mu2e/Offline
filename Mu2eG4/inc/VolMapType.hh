#ifndef Mu2eG4_VolMapType_hh
#define Mu2eG4_VolMapType_hh
//
// In the run data there is a data product that describes all physical
// volumes the run of G4.
//
// Given a pointer to a physical volume, return the index into the data product
// for that volume.
//
// $Id: VolMapType.hh,v 1.5 2013/09/27 14:56:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 14:56:14 $
//
// Original author Rob Kutschke
//

#include <map>
#include "cetlib/map_vector.h"

class G4VPhysicalVolume;

namespace mu2e{

  typedef std::map<G4VPhysicalVolume*,cet::map_vector_key> VolMapType;
  typedef std::map<G4VPhysicalVolume*,cet::map_vector_key>::iterator VolMapType_iterator;
  typedef std::map<G4VPhysicalVolume*,cet::map_vector_key>::const_iterator VolMapType_const_iterator;

}

#endif /* Mu2eG4_VolMapType_hh */
