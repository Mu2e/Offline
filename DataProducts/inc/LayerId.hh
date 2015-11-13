#ifndef TrackerGeom_LayerId_hh
#define TrackerGeom_LayerId_hh

//
// Identifier of one layer in a tracker.
//

//
// $Id: LayerId.hh,v 1.8 2012/05/14 19:20:45 brownd Exp $
// $Author: brownd $
// $Date: 2012/05/14 19:20:45 $
//
// Original author Rob Kutschke
//

#include <ostream>

#include "DataProducts/inc/SectorId.hh"

namespace mu2e {
  class LayerId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const LayerId& l );

  class LayerId{

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

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const DeviceId& getDeviceId () const{
      return _sid.getDeviceId();
    }
    const SectorId& getSectorId () const{
      return _sid;
    }

    int getDevice () const{
      return _sid.getDevice();
    }

    int getSector () const{
      return _sid.getSector();
    }

    int getLayer() const{
      return _layer;
    }

    bool operator==(LayerId const& rhs) const{
      return ( _sid == rhs._sid && _layer == rhs._layer );
    }

    bool operator!=(LayerId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator < (LayerId const& rhs) const {
      return _sid < rhs._sid || (_sid == rhs._sid && _layer < rhs._layer);
    }
  
    friend std::ostream& operator<<(std::ostream& ost,
                                    const LayerId& l ){
      ost << l._sid << " " << l._layer;
      return ost;
    }

private:

    SectorId _sid;
    int      _layer;

  };

} //namespace mu2e

#endif /* TrackerGeom_LayerId_hh */
