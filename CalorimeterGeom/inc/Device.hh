#ifndef DEVICE_HH
#define DEVICE_HH

//
// C++ includes
#include <vector>

//
// Mu2e includes
#include "CalorimeterGeom/inc/DeviceId.hh"
#include "CalorimeterGeom/inc/Layer.hh"

class Device{

  friend class Calorimeter;
  friend class CalorimeterMaker;

public:

  // A free function, returning void, that takes a const Device& as an argument.
  typedef void (*DeviceFunction)( const Device& s);

  Device():_id(-1){}
  Device( const DeviceId& id ):_id(id){}
  ~Device(){}
 
  // Compiler generated copy and assignment constructors
  // should be OK.
  
  const DeviceId Id() const { return _id;}

  const std::vector<Layer>& getLayers () const{ 
    return _layers;
  }

  const Layer& getLayer ( int n) const { 
    return _layers.at(n);
  }

  const Layer& getLayer ( const LayerId& lid ) const{
    return _layers.at(lid._layer);
  }

  int nLayers() const { return _layers.size(); }

  const Crystal& getCrystal ( const CrystalId& lid ) const{
    return _layers.at(lid.getLayer()).getCrystal(lid);
  }

#ifndef __CINT__

  // Loop over all crystals and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void Device::forAllCrystals ( F& f) const{
    for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
	  i !=e; ++i){
      i->forAllCrystals(f);
    }
  }

  // Loop over all layers and call F
  template <class F>
  inline void Device::forAllLayers ( F& f) const{
    for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
	  i !=e; ++i){
      f(*i);
    }
  }

#endif

protected:

  DeviceId _id;
  std::vector<Layer> _layers;

};

#endif
