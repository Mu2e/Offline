#ifndef TTrackerGeom_PlaneMFId_hh
#define TTrackerGeom_PlaneMFId_hh

//
// Identifier for a planeMF.
//

//
// $Id: PlaneMFId.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//

#include <ostream>
#include "TTrackerGeom/inc/StationId.hh"

namespace mu2e {

  struct PlaneMFId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PlaneMFId& s );

  class PlaneMFId{

  public:

    PlaneMFId():
      _sid(-1),
      _planeMF(-1){
    }

    PlaneMFId( StationId station,
              int planeMF
              ):
      _sid(station),
      _planeMF(planeMF){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const StationId& getStationId() const { return _sid; }
          int getStation()          const { return _sid; }
          int getPlaneMF()            const { return _planeMF; }

    bool operator==(PlaneMFId const& rhs) const{
      return ( _sid == rhs._sid && _planeMF == rhs._planeMF );
    }

    bool operator!=(PlaneMFId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator< (PlaneMFId const& rhs) const {
      if (_sid < rhs._sid) return true;
      if (_sid > rhs._sid) return false;
      if (_planeMF < rhs._planeMF) return true;
      return false;
    }

    friend std::ostream& operator<<(std::ostream& ost,
                                    const PlaneMFId& p ){
      ost << p._sid << " " << p._planeMF;
      return ost;
    }

  private:

    StationId _sid;
    int      _planeMF;

  };

}  //namespace mu2e

#endif /* TTrackerGeom_PlaneMFId_hh */
