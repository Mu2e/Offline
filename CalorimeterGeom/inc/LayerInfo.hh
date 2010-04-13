#ifndef LAYERINFO_HH
#define LAYERINFO_HH

// $Id: LayerInfo.hh,v 1.2 2010/04/13 17:15:37 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:15:37 $

// original authors Julie Managan and Robert Bernstein

namespace mu2e{
  namespace calorimeter{

    struct LayerInfo{

      // Allows different types for future use
      enum Stype {conductive, nonconductive, undefined};

    public:
      LayerInfo():
	_nCrystals(-1),
	_crystalType(undefined)
      {
      }
      LayerInfo( int nCrystals,
		 Stype crystalType
		 ):
	_nCrystals(nCrystals),
	_crystalType(crystalType){
      }
  
      ~LayerInfo  (){}
  
      // Compiler generated copy and assignment constructors
      // should be OK.
      int _nCrystals;
      Stype _crystalType;
  
    };


  } //namespace calorimeter
} //namespace mu2e

#endif
