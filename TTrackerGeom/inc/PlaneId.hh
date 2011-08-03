#ifndef TrackerGeom_PlaneId_hh
#define TrackerGeom_PlaneId_hh

//
// Identifier for a plane.
//

//
// $Id: PlaneId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>
#include "TTrackerGeom/inc/StationId.hh"

namespace mu2e {

  struct PlaneId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PlaneId& s );

  class PlaneId{

  public:

    PlaneId():
      _sid(-1),
      _plane(-1){
    }

    PlaneId( StationId station,
              int plane
              ):
      _sid(station),
      _plane(plane){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _sid; }
          int getStation()          const { return _sid; }
          int getPlane()            const { return _plane; }

    bool operator==(PlaneId const& rhs) const{
      return ( _sid == rhs._sid && _plane == rhs._plane );
    }

    bool operator!=(PlaneId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (PlaneId const& rhs) const {
      if (_sid < rhs._sid) return true;
      if (_sid > rhs._sid) return false;
      if (_plane < rhs._plane) return true;
      return false;
    }

    friend std::ostream& operator<<(std::ostream& ost,
                                    const PlaneId& p ){
      ost << p._sid << " " << p._plane;
      return ost;
    }

  private:

    StationId _sid;
    int      _plane;

  };

}  //namespace mu2e

#endif /* TrackerGeom_PlaneId_hh */
