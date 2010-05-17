//
// Hold information about one RSlice in the calorimeter.
//
//
// $Id: RSlice.cc,v 1.2 2010/05/17 21:47:33 genser Exp $
// $Author: genser $
// $Date: 2010/05/17 21:47:33 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <sstream>

#include "CalorimeterGeom/inc/RSlice.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#ifndef __CINT__ 

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {
  namespace calorimeter{

    RSlice::RSlice():
      _id(RSliceId()),
      _nCrystals(0),
      _orig(CLHEP::Hep3Vector(0.,0.,0.)),
      _delta(CLHEP::Hep3Vector(0.,0.,0.)){
    }

    RSlice::RSlice(const RSliceId& id,
		   int      nCrystals,
		   const CLHEP::Hep3Vector& origin,
		   const CLHEP::Hep3Vector& delta
		   ):
      _id(id),
      _nCrystals(nCrystals),
      _orig(origin),
      _delta(delta){
    }

    RSlice::RSlice(const RSliceId& id ):
      _id(id){
    }


    string RSlice::name( string const& base ) const{
      ostringstream os;

      os << base
	 << _id.getVane() << "_"
	 << _id.getZSlice() << "_"
	 << _id.getRSlice();
      return os.str();
    }

    /*
      void RSlice::fillPointers ( const Calorimeter& calorimeter ) const{
      _crystals.clear();
      for ( size_t i=0; i<_indices.size(); ++i ){
      CrystalIndex idx = _indices[i];
      const Crystal* crystal =  &calorimeter.getCrystal(idx);
      _crystals.push_back(crystal);
      crystal->fillPointers(calorimeter);
      }
      }
    */
  } //namespace calorimeter
} // namespace mu2e 
#endif
  
