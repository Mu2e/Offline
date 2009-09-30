#ifndef LAYERID_HH
#define LAYERID_HH

//
// Identifier of one layer in a tracker.
//

//
// $Id: LayerId.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <iostream>

#include "LTrackerGeom/inc/SectorId.hh"

namespace mu2e {

struct LayerId{

public:

  LayerId():
    _sid(SectorId()),
    _layer(-1){
  }

  LayerId( SectorId sector,
	   int layer
	   ):
    _sid(sector),
    _layer(layer){
  }

  LayerId( DeviceId device,
	   int sector,
	   int layer
	   ):
    _sid(SectorId(device,sector)),
    _layer(layer){
  }

  ~LayerId  (){
  }

  const DeviceId getDeviceId () const{
    return _sid._did;
  }
  const SectorId getSectorId () const{
    return _sid;
  }

  const int getDevice () const{
    return _sid._did;
  }

  const int getSector () const{
    return _sid._sector;
  }

  const int getLayer() const{
    return _layer;
  }

  bool operator==(const LayerId l) const{
    return ( _sid == l._sid && _layer == l._layer );
  }

  
  // Compiler generated copy and assignment constructors
  // should be OK.

  SectorId _sid;
  int _layer;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
				const LayerId& l ){
  ost << l._sid << " " << l._layer;
  return ost;
}

} //namespace mu2e

#endif
