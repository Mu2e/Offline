#ifndef TrackerGeom_PanelId_hh
#define TrackerGeom_PanelId_hh

//
// Identifier for a panel.
//

//
// $Id: PanelId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>
#include "TTrackerGeom/inc/FaceId.hh"

namespace mu2e {

  struct PanelId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PanelId& pan );

  class PanelId{

  public:

    PanelId():
      _fid(),
      _panel(-1){
    }

    PanelId( FaceId const & face,
              int panel
              ):
      _fid(face),
      _panel(panel){
    }

    PanelId( int station, int face, int panel )
      : _fid(station, face)
      , _panel(panel){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _fid.getStationId(); }
    const PlaneId&   getPlaneId()   const { return _fid.getPlaneId();   }
    const FaceId&    getFaceId()    const { return _fid;                }
          int   getStation()        const { return _fid.getStation();   }
          int   getPlane()          const { return _fid.getPlane();     }
          int   getFace()           const { return _fid.getFace();      }
          int   getPanel()          const { return _panel;              }


    bool operator==(PanelId const& rhs) const{
      return ( _fid == rhs._fid && _panel == rhs._panel );
    }

    bool operator!=(PanelId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (PanelId const& rhs) const {
      if (_fid < rhs._fid) return true;
      if (_fid > rhs._fid) return false;
      if (_panel < rhs._panel) return true;
      return false;
    }
    
    bool operator> (PanelId const& rhs) const {
      if (_fid > rhs._fid) return true;
      if (_fid < rhs._fid) return false;
      if (_panel > rhs._panel) return true;
      return false;
    }
    
    friend std::ostream& operator<<(std::ostream& ost,
                                    const PanelId& pan ){
      ost << pan._fid << " " << pan._panel;
      return ost;
    }

  private:

    FaceId _fid;
    int    _panel;

  };

}  //namespace mu2e

#endif /* TrackerGeom_PanelId_hh */
