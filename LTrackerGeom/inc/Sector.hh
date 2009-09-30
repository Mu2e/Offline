#ifndef SECTOR_HH
#define SECTOR_HH
//
// Holds information about one sector in a tracker.
//

//
// $Id: Sector.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "LTrackerGeom/inc/SectorId.hh"
#include "LTrackerGeom/inc/Layer.hh"

#include "CLHEP/Vector/ThreeVector.h"
#ifndef __CINT__
#include "boost/bind.hpp"
#endif


namespace mu2e {

class Sector{

  friend class Device;
  friend class LTracker;
  friend class LTrackerMaker;

public:

  // A free function, returning void, that takes a const Sector& as an argument.
  typedef void (*SectorFunction)( const Sector& s);

  Sector():_id(SectorId(-1,-1)){};
  Sector( const SectorId& id ):_id(id){};
  ~Sector(){}
 
  // Compiler generated copy and assignment constructors
  // should be OK.
  
  const SectorId& Id() const { return _id;}

  const std::vector<Layer>& getLayers() const{ 
    return _layers;
  }

  int nLayers() const{ 
    return _layers.size();
  }

  const Layer& getLayer ( int n ) const { 
    return _layers.at(n);
  }

  const Layer& getLayer ( const LayerId& lid) const { 
    return _layers.at(lid.getLayer());
  }

  const Straw& getStraw ( const StrawId& sid ) const{
    return _layers.at(sid.getLayer()).getStraw(sid);
  }


#ifndef __CINT__

  template <class F>
  void for_each_layer( F f) const{
    std::for_each ( _layers.begin(),
		    _layers.end(),
		    f);
  }

  template <class F>
  void for_each_straw( F f) const {
    for_each_layer( boost::bind( Layer::for_each<F>, _1, f));
  }

  // Loop over all straws and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void Sector::forAllStraws ( F& f) const{
    for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
	  i !=e; ++i){
      i->forAllStraws(f);
    }
  }

  template <class F>
  inline void Sector::forAllLayers ( F& f) const{
    for ( std::vector<Layer>::const_iterator i=_layers.begin(), e=_layers.end();
	  i !=e; ++i){
      f(*i);
    }
  }


#endif
  
protected:
  
  SectorId _id;
  std::vector<Layer> _layers;

  // Vertices of enclosing polygon.
  std::vector<CLHEP::Hep3Vector> corners;

};

}  //namespace mu2e

#endif

