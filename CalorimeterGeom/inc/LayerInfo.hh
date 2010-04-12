#ifndef LAYERINFO_HH
#define LAYERINFO_HH


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

#endif
