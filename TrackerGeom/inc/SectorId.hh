#ifndef SECTORID_HH
#define SECTORID_HH

//
// Identifier for a sector.
//

//
// $Id: SectorId.hh,v 1.2 2010/04/14 14:16:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/14 14:16:41 $
//
// Original author Rob Kutschke
//

#include <ostream>
#include "TrackerGeom/inc/DeviceId.hh"

namespace mu2e { 

struct SectorId{

public:

  SectorId():
    _did(-1),
    _sector(-1){
  }
  
  SectorId( DeviceId device,
	    int sector
	   ):
    _did(device),
    _sector(sector){
  }
  
  ~SectorId  (){
  }
  
  // Compiler generated copy and assignment constructors
  // should be OK.

  const int getDeviceId() const {
    return _did;
  }

  const int getDevice() const {
    return _did;
  }

  const int getSector() const {
    return _sector;
  }

  bool operator==(SectorId const& rhs) const{
    return ( _did == rhs._did && _sector == rhs._sector );
  }

  bool operator!=(SectorId const& rhs) const{
    return !( *this == rhs);
  }
  
  DeviceId _did;
  int32_t _sector;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
				const SectorId& s ){
  ost << s._did << " " << s._sector;
  return ost;
}

}  //namespace mu2e

#endif
