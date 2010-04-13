#ifndef LAYERID_HH
#define LAYERID_HH
// $Id: LayerId.hh,v 1.3 2010/04/13 17:15:33 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:15:33 $

// original authors Julie Managan and Robert Bernstein

namespace mu2e{
  namespace calorimeter{

//
// C++ includes
#include <iostream>

//
// Mu2e includes
#include "CalorimeterGeom/inc/DeviceId.hh"

struct LayerId{

public:

  LayerId():
    _did(DeviceId()),
    _layer(-1){
  }

  LayerId( DeviceId device,
	   int layer
	   ):
    _did(device),
    _layer(layer){
  }

  ~LayerId  (){
  }

  const DeviceId getDeviceId () const{
    return _did;
  }

  const int getDevice () const{
    return _did;
  }

  const int getLayer() const{
    return _layer;
  }

  bool operator==(const LayerId l) const{
    return ( _did == l._did && _layer == l._layer );
  }

  
  // Compiler generated copy and assignment constructors
  // should be OK.

  DeviceId _did;
  int _layer;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
				const LayerId& l ){
  ost << l._did << " " << l._layer;
  return ost;
}

  } //namespace calorimeter
} //namespace mu2e

#endif
