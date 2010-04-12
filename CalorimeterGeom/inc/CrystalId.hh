#ifndef CRYSTALID_HH
#define CRYSTALID_HH

#include "Calorimeter/inc/LayerId.hh"

struct CrystalId{

public:

  CrystalId():
    _lid(LayerId()),
    _n(-1){
  }
  
  CrystalId( LayerId layer,
	   int n
	   ):
    _lid(layer),
    _n(n){
  }
  
  CrystalId( DeviceId deviceid,
	   int layer,
	   int n
	   ):
    _lid(deviceid,layer),
    _n(n){
  }

  ~CrystalId  (){
  }

  // Compiler generated copy c'tor and assignment
  // operators should be should be OK.

  const DeviceId& getDeviceId() const {
    return _lid._did;
  }

  const LayerId& getLayerId() const {
    return _lid;
  }
  
  const int getDevice() const{
    return _lid._did;
  }

  const int getLayer() const{
    return _lid._layer;
  }

  const int getCrystal() const{
    return _n;
  }

  bool operator==(const CrystalId c) const{
    return ( _lid == c._lid && _n == c._n );
  }

  LayerId _lid;
  int _n;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
				const CrystalId& c ){
  ost << "Crystal Id: ("
      << c.getLayerId() << " "
      << c._n
      << " )";
  return ost;
}

#endif
