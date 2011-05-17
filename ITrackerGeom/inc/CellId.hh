#ifndef ITrackerGeom_CellId_hh
#define ITrackerGeom_CellId_hh

#include "ITrackerGeom/inc/WireId.hh"

namespace mu2e { 

struct CellId{

public:

  CellId():
    _swid()
  {
  }
  
  CellId( WireId &swid):
    _swid(swid._lid, swid._n)
  {
  }

  CellId( ITLayerId *layer,
           int n
           ):
    _swid(WireId(layer,n))
    {
  }

  ~CellId  (){
  }

  const ITLayerId& getLayerId() const {
    return _swid.getLayerId();
  }
  
  const int getLayer() const{
    return _swid.getLayer();
  }

  const int getCell() const{
    return _swid.getWire();
  }

  bool operator==(const CellId c) const{
    return ( _swid == c._swid );
  }

  WireId _swid;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
                                const CellId& c ){
  ost << "Cell Id: ("
      << c.getLayerId() << " "
      << c._swid._n
      << " )";
  return ost;
}

}
#endif /* ITrackerGeom_CellId_hh */
