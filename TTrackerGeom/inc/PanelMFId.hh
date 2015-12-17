#ifndef TrackerGeom_PanelMFId_hh
#define TrackerGeom_PanelMFId_hh

//
// Identifier for a panelMF.
//

//
// $Id: PanelMFId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>
#include "TTrackerGeom/inc/FaceId.hh"

namespace mu2e {

  struct PanelMFId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PanelMFId& pan );

  class PanelMFId{

  public:

    PanelMFId():
      _fid(),
      _panelMF(-1){
    }

    PanelMFId( FaceId const & face,
              int panelMF
              ):
      _fid(face),
      _panelMF(panelMF){
    }

    PanelMFId( int station, int face, int panelMF )
      : _fid(station, face)
      , _panelMF(panelMF){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _fid.getStationId(); }
    const PlaneMFId&   getPlaneMFId()   const { return _fid.getPlaneMFId();   }
    const FaceId&    getFaceId()    const { return _fid;                }
          int   getStation()        const { return _fid.getStation();   }
          int   getPlaneMF()          const { return _fid.getPlaneMF();     }
          int   getFace()           const { return _fid.getFace();      }
          int   getPanelMF()          const { return _panelMF;              }


    bool operator==(PanelMFId const& rhs) const{
      return ( _fid == rhs._fid && _panelMF == rhs._panelMF );
    }

    bool operator!=(PanelMFId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (PanelMFId const& rhs) const {
      if (_fid < rhs._fid) return true;
      if (_fid > rhs._fid) return false;
      if (_panelMF < rhs._panelMF) return true;
      return false;
    }
    
    bool operator> (PanelMFId const& rhs) const {
      if (_fid > rhs._fid) return true;
      if (_fid < rhs._fid) return false;
      if (_panelMF > rhs._panelMF) return true;
      return false;
    }
    
    friend std::ostream& operator<<(std::ostream& ost,
                                    const PanelMFId& pan ){
      ost << pan._fid << " " << pan._panelMF;
      return ost;
    }

  private:

    FaceId _fid;
    int    _panelMF;

  };

}  //namespace mu2e

#endif /* TrackerGeom_PanelMFId_hh */
