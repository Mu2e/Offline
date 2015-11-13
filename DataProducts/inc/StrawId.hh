#ifndef TrackerGeom_StrawId_hh
#define TrackerGeom_StrawId_hh
//
// Identifier of one straw in a tracker.
//

//
// $Id: StrawId.hh,v 1.7 2011/05/20 19:18:45 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 19:18:45 $
//
// Original author Rob Kutschke
//
#include <ostream>
#include "DataProducts/inc/LayerId.hh"

namespace mu2e {

  class StrawId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const StrawId& s );

class StrawId{

public:

  StrawId():
    _lid(),
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

  // Use compiler-generated copy c'tor, copy assignment, and d'tor.

  const DeviceId& getDeviceId() const {
    return _lid.getDeviceId();
  }

  const SectorId& getSectorId() const {
    return _lid.getSectorId();
  }

  const LayerId& getLayerId() const {
    return _lid;
  }

  int getDevice() const{
    return _lid.getDevice();
  }

  int getSector() const{
    return _lid.getSector();
  }

  int getLayer() const{
    return _lid.getLayer();
  }

  int getStraw() const{
    return _n;
  }

  bool operator==( StrawId const& rhs) const{
    return ( _lid == rhs._lid && _n == rhs._n );
  }

  bool operator!=( StrawId const& rhs) const{
    return !( *this == rhs);
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const StrawId& s ){
    ost << "Straw Id: ("
        << s.getLayerId() << " "
        << s._n
        << " )";
    return ost;
  }

private:

  LayerId _lid;
  int     _n;

};

}
#endif /* TrackerGeom_StrawId_hh */
