#ifndef TrackerGeom_ZLayerId_hh
#define TrackerGeom_ZLayerId_hh

//
// Identifier for a layer ordered by Z (ZLayer).
// Although all 3 panels in a face have layers at the identical Z,
// we consider each such layer as a distinct ZLayer.
//

//
// $Id: ZLayerId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>
#include "TTrackerGeom/inc/PanelId.hh"

namespace mu2e {

  struct ZLayerId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const ZLayerId& zlay );

  class ZLayerId{

  public:

    ZLayerId():
      _panid(),
      _zlayer(-1){
    }

    ZLayerId( PanelId const & pan,
              int zlay
              ):
      _panid(pan),
      _zlayer(zlay){
    }

    ZLayerId( int station, int face, int panel, int zlay )
      : _panid(station, face, panel)
      , _zlayer(zlay){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _panid.getStationId(); }
    const PlaneId&   getPlaneId()   const { return _panid.getPlaneId();   }
    const FaceId&    getFaceId()    const { return _panid.getFaceId();    }
    const PanelId&   getPanelId()   const { return _panid;                }
          int   getStation()        const { return _panid.getStation();   }
          int   getPlane()          const { return _panid.getPlane();     }
          int   getFace()           const { return _panid.getFace();      }
          int   getPanel()          const { return _panid.getPanel();     }
          int   getZLayer()         const { return _zlayer;               }


    bool operator==(ZLayerId const& rhs) const{
      return ( _panid == rhs._panid && _zlayer == rhs._zlayer );
    }

    bool operator!=(ZLayerId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (ZLayerId const& rhs) const {
      if (_panid < rhs._panid) return true;
      if (_panid > rhs._panid) return false;
      if (_zlayer < rhs._zlayer) return true;
      return false;
    }
    
    bool operator> (ZLayerId const& rhs) const {
      if (_panid > rhs._panid) return true;
      if (_panid < rhs._panid) return false;
      if (_zlayer > rhs._zlayer) return true;
      return false;
    }
    
    friend std::ostream& operator<<(std::ostream& ost,
                                    const ZLayerId& zlay ){
      ost << zlay._panid << " " << zlay._zlayer;
      return ost;
    }

  private:

    PanelId _panid;
    int     _zlayer;

  };

}  //namespace mu2e

#endif /* TrackerGeom_ZLayerId_hh */
