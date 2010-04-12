#ifndef LAYERID_HH
#define LAYERID_HH

//
// C++ includes
#include <iostream>

//
// Mu2e includes
#include "Calorimeter/inc/DeviceId.hh"

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
#endif
