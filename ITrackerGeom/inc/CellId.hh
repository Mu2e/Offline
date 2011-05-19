#ifndef ITrackerGeom_CellId_hh
#define ITrackerGeom_CellId_hh

#include "ITrackerGeom/inc/WireId.hh"

namespace mu2e {

class CellId;
inline std::ostream& operator<<(std::ostream& ost,
                                const CellId& c );

class CellId{

public:

  CellId():
    _swid()
  {
  }

  CellId( WireId &swid):
    _swid(swid)
  {
  }

  CellId( ITLayerId *layer,
          int n
        ):
    _swid(layer,n)
  {
  }

  // use compiler-generated copy c'tor, copy assignment, and d'tor

  const ITLayerId& getLayerId() const {
    return _swid.getLayerId();
  }

  int getLayer() const{
    return _swid.getLayer();
  }

  int getCell() const{
    return _swid.getWire();
  }

  bool operator==(const CellId c) const{
    return ( _swid == c._swid );
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const CellId& c ){
    ost << "Cell Id: ("
        << c.getLayerId() << " "
        << c._swid.getWire()
        << " )";
    return ost;
  }

private:

  WireId _swid;

};

}
#endif /* ITrackerGeom_CellId_hh */
