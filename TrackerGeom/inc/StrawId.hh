#ifndef TrackerGeom_StrawId_hh
#define TrackerGeom_StrawId_hh
//
// Identifier of one straw in a tracker.
//

//
// $Id: StrawId.hh,v 1.6 2011/05/18 15:47:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:47:40 $
//
// Original author Rob Kutschke
//
#include <ostream>
#include "TrackerGeom/inc/LayerId.hh"

namespace mu2e {

struct StrawId{

public:

  StrawId():
    _lid(LayerId()),
    _n(-1){
  }

  StrawId( LayerId layer,
           int n
           ):
    _lid(layer),
    _n(n){
  }

  StrawId( SectorId sectorid,
           int layer,
           int n
           ):
    _lid(sectorid,layer),
    _n(n){
  }

  StrawId( DeviceId device,
           int section,
           int layer,
           int n
           ):
    _lid(LayerId(device,section,layer)),
    _n(n){
  }

  ~StrawId  (){
  }

  // Compiler generated copy c'tor and assignment
  // operators should be should be OK.

  const DeviceId& getDeviceId() const {
    return _lid._sid._did;
  }

  const SectorId& getSectorId() const {
    return _lid._sid;
  }

  const LayerId& getLayerId() const {
    return _lid;
  }

  const int getDevice() const{
    return _lid._sid._did;
  }

  const int getSector() const{
    return _lid._sid._sector;
  }

  const int getLayer() const{
    return _lid._layer;
  }

  const int getStraw() const{
    return _n;
  }

  bool operator==( StrawId const& rhs) const{
    return ( _lid == rhs._lid && _n == rhs._n );
  }

  bool operator!=( StrawId const& rhs) const{
    return !( *this == rhs);
  }


  LayerId _lid;
  int     _n;

};

inline std::ostream& operator<<(std::ostream& ost,
                                const StrawId& s ){
  ost << "Straw Id: ("
      << s.getLayerId() << " "
      << s._n
      << " )";
  return ost;
}

}
#endif /* TrackerGeom_StrawId_hh */
