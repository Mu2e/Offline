//
// $Id: WireId.hh,v 1.7 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//
// Original author G. Tassielli
//

#ifndef ITrackerGeom_WireId_hh
#define ITrackerGeom_WireId_hh

#include "ITrackerGeom/inc/ITLayerId.hh"

namespace mu2e {

class WireId;
inline std::ostream& operator<<(std::ostream& ost,
                                const WireId& w );

class WireId{

        friend class ITLayerId;

public:

  WireId():
    _lid(NULL),
    _n(-1){
  }

  WireId( ITLayerId *layer,
           int n
           ):
    _lid(layer),
    _n(n){
  }

  // use compiler-generated copy c'tor, copy assignment, and d'tor

  const ITLayerId& getLayerId() const {
    return *_lid;
  }

  int getLayer() const{
    return _lid->getLayer();
  }

  int getWire() const{
    return _n;
  }

  int getUId() const{
	int wUid = ((_lid->getSuperLayer()+1)*100000)+((_lid->getLayer()+1)*1000);
	wUid*=10;
	wUid+=_n;
    return wUid;
  }

  bool operator==(const WireId w) const{
    return ( *_lid == *(w._lid) && _n == w._n );
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const WireId& w ){
    ost << "Wire Id: ("
        << w.getLayerId() << " "
        << w._n
        << " )";
    return ost;
  }

private:

  ITLayerId *_lid;
  int _n;

};

}
#endif /* ITrackerGeom_WireId_hh */
