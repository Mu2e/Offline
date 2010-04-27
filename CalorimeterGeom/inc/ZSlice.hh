#ifndef ZSLICE_HH
#define ZSLICE_HH
//
// Holds information about one zslice in a calorimeter.
//

//
// $Id: ZSlice.hh,v 1.1 2010/04/27 18:44:12 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:44:12 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <vector>

#include "CalorimeterGeom/inc/ZSliceId.hh"
#include "CalorimeterGeom/inc/RSlice.hh"

#include "CLHEP/Vector/ThreeVector.h"
#ifndef __CINT__
#include "boost/bind.hpp"
#endif


namespace mu2e {
  namespace calorimeter{
class ZSlice{

  friend class Vane;
  friend class Calorimeter;
  friend class CalorimeterMaker;

public:

  ZSlice():_id(ZSliceId(-1,-1)){};
  ZSlice( const ZSliceId& id ):_id(id){};
  ~ZSlice(){}
 
  // Compiler generated copy and assignment constructors
  // should be OK.
  
  const ZSliceId& Id() const { return _id;}

  const std::vector<RSlice>& getRSlices() const{ 
    return _rslices;
  }

  int nRSlices() const{ 
    return _rslices.size();
  }

  const RSlice& getRSlice ( int n ) const { 
    return _rslices.at(n);
  }

  const RSlice& getRSlice ( const RSliceId& rslid) const { 
    return _rslices.at(rslid.getRSlice());
  }

  const Crystal& getCrystal ( const CrystalId& cid ) const{
    return _rslices.at(cid.getRSlice()).getCrystal(cid);
  }

  // Formatted string embedding the id of the sector.
  std::string name( std::string const& base ) const;

  const std::vector<double>& boxHalfLengths() const { return _boxHalfLengths; }

  const double         boxRxAngle()     const { return _boxRxAngle;     }
  const double         boxRyAngle()     const { return _boxRyAngle;     }
  const double         boxRzAngle()     const { return _boxRzAngle;     }
  const Hep3Vector&    boxOffset()      const { return _boxOffset;      }

  std::vector<CLHEP::Hep3Vector> const& getBasePosition() const{
    return _basePosition;
  }

  CLHEP::Hep3Vector const& getBaseDelta() const{
    return _baseDelta;
  }

#ifndef __CINT__

  template <class F>
  void for_each_rslice( F f) const{
    std::for_each ( _rslices.begin(),
		    _rslices.end(),
		    f);
  }

  template <class F>
  void for_each_crystal( F f) const {
    for_each_rslice( boost::bind( RSlice::for_each<F>, _1, f));
  }

  // Loop over all crystals and call F.
  // F can be a class with an operator() or a free function.
  template <class F>
  inline void ZSlice::forAllCrystals ( F& f) const{
    for ( std::vector<RSlice>::const_iterator i=_rslices.begin(), e=_rslices.end();
	  i !=e; ++i){
      i->forAllCrystals(f);
    }
  }

  template <class F>
  inline void ZSlice::forAllRSlices ( F& f) const{
    for ( std::vector<RSlice>::const_iterator i=_rslices.begin(), e=_rslices.end();
	  i !=e; ++i){
      f(*i);
    }
  }


#endif
  
protected:
  
  ZSliceId _id;
  std::vector<RSlice> _rslices;

  // Vertices of enclosing polygon.
  std::vector<CLHEP::Hep3Vector> corners;

  // Properties of the enclosing logical volume (box).

  // Half lengths of the logical box.
  std::vector<double> _boxHalfLengths;

  std::vector<CLHEP::Hep3Vector> _basePosition;
  CLHEP::Hep3Vector _baseDelta;

  // Rotations and offsets to place the logical box.
  // placedshape = ( offset + RZ*RX*RY*shape );
  //
  double _boxRxAngle;
  double _boxRyAngle;
  double _boxRzAngle;
  CLHEP::Hep3Vector _boxOffset;

};
  } //namespace calorimeter
}  //namespace mu2e

#endif

