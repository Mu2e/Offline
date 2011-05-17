#ifndef TrackerGeom_LayerId_hh
#define TrackerGeom_LayerId_hh

//
// Identifier of one layer in a tracker.
//

//
// $Id: LayerId.hh,v 1.4 2011/05/17 15:41:37 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:37 $
//
// Original author Rob Kutschke
//

#include <ostream>

#include "TrackerGeom/inc/SectorId.hh"

namespace mu2e {

  struct LayerId{

  public:

    LayerId():
      _sid(SectorId()),
      _layer(-1){
    }

    LayerId( SectorId sector,
             int layer
             ):
      _sid(sector),
      _layer(layer){
    }

    LayerId( DeviceId device,
             int sector,
             int layer
             ):
      _sid(SectorId(device,sector)),
      _layer(layer){
    }

    ~LayerId  (){
    }

    const DeviceId getDeviceId () const{
      return _sid._did;
    }
    const SectorId getSectorId () const{
      return _sid;
    }

    const int getDevice () const{
      return _sid._did;
    }

    const int getSector () const{
      return _sid._sector;
    }

    const int getLayer() const{
      return _layer;
    }

    bool operator==(LayerId const& rhs) const{
      return ( _sid == rhs._sid && _layer == rhs._layer );
    }

    bool operator!=(LayerId const& rhs) const{
      return !( *this == rhs);
    }

  
    // Compiler generated copy and assignment constructors
    // should be OK.

    SectorId _sid;
    int32_t _layer;
  
  };

  inline std::ostream& operator<<(std::ostream& ost, 
                                  const LayerId& l ){
    ost << l._sid << " " << l._layer;
    return ost;
  }

} //namespace mu2e

#endif /* TrackerGeom_LayerId_hh */
