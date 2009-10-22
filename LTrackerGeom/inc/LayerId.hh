#ifndef LAYERID_HH
#define LAYERID_HH

//
// Identifier of one layer in a tracker.
//

//
// $Id: LayerId.hh,v 1.2 2009/10/22 16:43:11 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/22 16:43:11 $
//
// Original author Rob Kutschke
//

#include <ostream>

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

  bool operator==(LayerId const& rhs) const{
    return ( _sid == rhs._sid && _layer == rhs._layer );
  }

  bool operator!=(LayerId const& rhs) const{
    return !( *this == rhs);
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
