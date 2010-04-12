#ifndef CALORIMETER_HH
#define CALORIMETER_HH

//
// C++ includes
#include <vector>

//
// Mu2e includes
#include "CalorimeterGeom/inc/Device.hh"

class Calorimeter{

  friend class CalorimeterMaker;

public:
  Calorimeter(){}
  ~Calorimeter(){};
 
  // Compiler generated copy and assignment constructors
  // should be OK.

  enum CalorimeterDeviceId { undefined=-1, vane0, vane1, vane2, vane3};

  double r0() const { return _r0;}
  double rInscribed() const { return _rInscribed;}

  int nVanes() const { return _devices.size(); }

  double crystalHalfSide() const{
    return getCrystal(CrystalId(0,0,0)).getDetail().yhalfLength();
  }

  // Check for legal identifiers.
  bool isLegal(DeviceId d) const{
    return ( d>-1 && 
	     std::vector<Device>::size_type(d) <_devices.size() 
	     );
  };

  typedef std::vector<Device>::size_type stypeLayer;
  bool isLegal(const LayerId& lid ) const{
    return ( isLegal(lid._did) &&
	     lid._layer > -1   &&
	     std::vector<Layer>::size_type(lid._layer) < getDevice(lid._did).getLayers().size()
	     );
  }

  bool isLegal(const CrystalId& cid) const{
    return ( isLegal(cid._lid) &&
	     cid._n > -1       &&
	     std::vector<Crystal>::size_type(cid._n) < getLayer(cid._lid).getCrystals().size()
	     );
  }

  // Accessors
  const std::vector<Device>& getDevices() const{ 
    return _devices;
  }

  const Device& getDevice ( DeviceId id) const{ 
    return _devices.at(id);
  }

  const Layer& getLayer ( const LayerId& lid ) const{
    return _devices.at(lid.getDevice()).getLayer(lid);
  }

  const Crystal& getCrystal ( const CrystalId& cid ) const{
    return _devices.at(cid.getDevice()).getCrystal(cid);
  }

  const std::vector<Crystal>& getAllCrystals() const {return _allCrystals;}

#ifndef __CINT__

  // Loop over all crystals and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void Calorimeter::forAllCrystals ( F& f) const{
    for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	  i !=e; ++i){
      i->forAllCrystals(f);
    }
  }

  template <class F>
  inline void Calorimeter::forAllLayers ( F& f) const{
    for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	  i !=e; ++i){
      i->forAllLayers(f);
    }
  }

  template <class F>
  inline void Calorimeter::forAllDevices ( F& f) const{
    for ( std::vector<Device>::const_iterator i=_devices.begin(), e=_devices.end();
	  i !=e; ++i){
      f(*i);
    }
  }

#endif


protected:

  // Nominal values.  
  // _r0 = Nominal radius of the center of the vane.
  double _r0;
  double _rInscribed;

  // Detailed info about each type of crystal.
  std::vector<CrystalDetail> _crystalDetail;

  // A Calorimeter is made of four devices,vanes 0-3.
  std::vector<Device> _devices;

  // There will be pointers to the objects in this container.
  std::vector<Crystal>  _allCrystals;

};

#endif
