#ifndef ITrackerGeom_WireId_hh
#define ITrackerGeom_WireId_hh

#include "ITrackerGeom/inc/ITLayerId.hh"

namespace mu2e { 

struct WireId{

        friend struct ITLayerId;

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
  
  ~WireId  (){
  }

  const ITLayerId& getLayerId() const {
    return *_lid;
  }
  
  const int getLayer() const{
    return _lid->_id;
  }

  const int getWire() const{
    return _n;
  }

  bool operator==(const WireId w) const{
    return ( *_lid == *(w._lid) && _n == w._n );
  }

  ITLayerId *_lid;
  int _n;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
                                const WireId& w ){
  ost << "Wire Id: ("
      << w.getLayerId() << " "
      << w._n
      << " )";
  return ost;
}

}
#endif /* ITrackerGeom_WireId_hh */
