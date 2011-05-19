#ifndef ITrackerGeom_SuperLayerInfo_hh
#define ITrackerGeom_SuperLayerInfo_hh

namespace mu2e {

class SuperLayerInfo{

  // The cells type
  enum Stype {undefined=-1, hexagonal, square};

public:
  SuperLayerInfo():
    _nLayers(-1),
    _cellType(undefined)
  {
  }
  SuperLayerInfo( int nLayers, Stype cellType
             ):
    _nLayers(nLayers),
    _cellType(cellType){
  }

  // use compiler-generated copy c'tor, copy assignment, and d'tor

private:
  int _nLayers;
  Stype _cellType;

};

}  //namespace mu2e

#endif /* ITrackerGeom_SuperLayerInfo_hh */
