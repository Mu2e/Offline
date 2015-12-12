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

#include "DataProducts/inc/PanelId.hh"

namespace mu2e {
  class LayerId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const LayerId& l );

  class LayerId{

  public:

    LayerId():
      _panelId(PanelId()),
      _layer(-1){
    }

    LayerId( PanelId panel,
             int layer
             ):
      _panelId(panel),
      _layer(layer){
    }

    LayerId( PlaneId plane,
             int panel,
             int layer
             ):
      _panelId(PanelId(plane,panel)),
      _layer(layer){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const PlaneId& getPlaneId () const{
      return _panelId.getPlaneId();
    }
    const PanelId& getPanelId () const{
      return _panelId;
    }

    int getPlane () const{
      return _panelId.getPlane();
    }

    int getStation() const{
      return _panelId.getStation();
    }

    int getPanel () const{
      return _panelId.getPanel();
    }

    int getLayer() const{
      return _layer;
    }

    bool operator==(LayerId const& rhs) const{
      return ( _panelId == rhs._panelId && _layer == rhs._layer );
    }

    bool operator!=(LayerId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator < (LayerId const& rhs) const {
      return _panelId < rhs._panelId || (_panelId == rhs._panelId && _layer < rhs._layer);
    }
  
    friend std::ostream& operator<<(std::ostream& ost,
                                    const LayerId& l ){
      ost << l._panelId << " " << l._layer;
      return ost;
    }

private:

    PanelId _panelId;
    int     _layer;

  };

} //namespace mu2e

#endif /* TrackerGeom_LayerId_hh */
