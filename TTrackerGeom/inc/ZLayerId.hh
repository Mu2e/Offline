#ifndef TrackerGeom_ZLayerId_hh
#define TrackerGeom_ZLayerId_hh

//
// Identifier for a layer ordered by Z (ZLayer).
// Although all 3 panelMFs in a face have layers at the identical Z,
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
#include "TTrackerGeom/inc/PanelMFId.hh"

namespace mu2e {

  struct ZLayerId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const ZLayerId& zlay );

  class ZLayerId{

  public:

    ZLayerId():
      _panMFid(),
      _zlayer(-1){
    }

    ZLayerId( PanelMFId const & pan,
              int zlay
              ):
      _panMFid(pan),
      _zlayer(zlay){
    }

    ZLayerId( int station, int face, int panelMF, int zlay )
      : _panMFid(station, face, panelMF)
      , _zlayer(zlay){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _panMFid.getStationId(); }
    const PlaneMFId&   getPlaneMFId()   const { return _panMFid.getPlaneMFId();   }
    const FaceId&    getFaceId()    const { return _panMFid.getFaceId();    }
    const PanelMFId&   getPanelMFId()   const { return _panMFid;                }
          int   getStation()        const { return _panMFid.getStation();   }
          int   getPlaneMF()          const { return _panMFid.getPlaneMF();     }
          int   getFace()           const { return _panMFid.getFace();      }
          int   getPanelMF()          const { return _panMFid.getPanelMF();     }
          int   getZLayer()         const { return _zlayer;               }


    bool operator==(ZLayerId const& rhs) const{
      return ( _panMFid == rhs._panMFid && _zlayer == rhs._zlayer );
    }

    bool operator!=(ZLayerId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (ZLayerId const& rhs) const {
      if (_panMFid < rhs._panMFid) return true;
      if (_panMFid > rhs._panMFid) return false;
      if (_zlayer < rhs._zlayer) return true;
      return false;
    }
    
    bool operator> (ZLayerId const& rhs) const {
      if (_panMFid > rhs._panMFid) return true;
      if (_panMFid < rhs._panMFid) return false;
      if (_zlayer > rhs._zlayer) return true;
      return false;
    }
    
    friend std::ostream& operator<<(std::ostream& ost,
                                    const ZLayerId& zlay ){
      ost << zlay._panMFid << " " << zlay._zlayer;
      return ost;
    }

  private:

    PanelMFId _panMFid;
    int     _zlayer;

  };

}  //namespace mu2e

#endif /* TrackerGeom_ZLayerId_hh */
