#ifndef TrackerGeom_FaceId_hh
#define TrackerGeom_FaceId_hh

//
// Identifier for a face.
//

//
// $Id: FaceId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>

#include "TTrackerGeom/inc/StationId.hh"
#include "TTrackerGeom/inc/PlaneMFId.hh"

namespace mu2e {

  struct FaceId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const FaceId& f );

  class FaceId{

  public:

    FaceId():
      _sid(-1),
      _pid(),
      _face(-1){
    }

    FaceId( StationId station, int face ):
      _sid(station),
      _pid(station,face/2),
      _face(face){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _sid;  }
    const PlaneMFId&   getPlaneMFId()   const { return _pid;  }
          int   getStation()        const { return _sid;  }
          int   getPlaneMF()          const { return _pid.getPlaneMF();  }
          int   getFace()           const { return _face; }

    bool operator==(FaceId const& rhs) const{
      return ( _sid == rhs._sid && _face == rhs._face );
    }

    bool operator!=(FaceId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (FaceId const& rhs) const {
      if (_sid < rhs._sid) return true;
      if (_sid > rhs._sid) return false;
      if (_face < rhs._face) return true;
      return false;
    }

    bool operator> (FaceId const& rhs) const {
      if (_sid > rhs._sid) return true;
      if (_sid < rhs._sid) return false;
      if (_face > rhs._face) return true;
      return false;
    }

    friend std::ostream& operator<<(std::ostream& ost,
                                    const FaceId& f ){
      ost << f._sid << " " << f._face;
      return ost;
    }

  private:

    StationId _sid;
    PlaneMFId   _pid;
    int       _face;

  };

}  //namespace mu2e

#endif /* TrackerGeom_FaceId_hh */
