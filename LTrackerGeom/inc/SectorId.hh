#ifndef SECTORID_HH
#define SECTORID_HH

//
// Identifier for a sector.
//

//
//
// $Id: SectorId.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <iostream>
#include "LTrackerGeom/inc/DeviceId.hh"

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

  bool operator==(const SectorId s) const{
    return ( _did == s._did && _sector == s._sector );
  }
  
  DeviceId _did;
  int _sector;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
				const SectorId& s ){
  ost << s._did << " " << s._sector;
  return ost;
}

}  //namespace mu2e

#endif
