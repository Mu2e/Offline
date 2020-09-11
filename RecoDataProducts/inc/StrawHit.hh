#ifndef RecoDataProducts_StrawHit_hh
#define RecoDataProducts_StrawHit_hh
//
// First version of a hit as described by Mu2e-doc-900.
//
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <vector>
#include "DataProducts/inc/TrkTypes.hh"

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"

namespace mu2e {

  struct StrawHit{
  private:

    StrawId      _strawId;       
    TrkTypes::TDCTimes	    _time;             // (ns)
    TrkTypes::TOTTimes        _tot;               // (ns)
    float           _energyDep;        // (MeV)

  public:

    StrawHit():
      _strawId(StrawId(-1)),
      _time{0.0,0.0},
      _tot{0.0,0.0},
      _energyDep(0.){
    }

// Constructor for a hit that came from an unpacked digi, either
    // from data or from the full MC chain.
    StrawHit( StrawId       strawId,
              TrkTypes::TDCTimes const& time,
              TrkTypes::TOTTimes const& tot,
              float            energyDep  ):
      _strawId(strawId),_time(time),_tot(tot),
      _energyDep(energyDep) {
    }

    // Accessors
    StrawId strawId() const { return _strawId; }
    float      time(StrawEnd end=StrawEnd::cal)       const { return _time[end];}
    // return the earliest time
    float      dt()         const { return _time[StrawEnd::cal] - _time[StrawEnd::hv]; }
    float      TOT(StrawEnd end=StrawEnd::cal)       const { return _tot[end];}
    float      energyDep()  const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
        bool operator==(StrawHit const& other) const {
	  return (_strawId==other._strawId&&
	      _time[0]==other._time[0]&&
	      _time[1]==other._time[1]&&
	      _tot[0]==other._tot[0]&&
	      _tot[1]==other._tot[1]&&
	      _energyDep==other._energyDep);
	}
    bool operator<( const StrawHit other) const{
      return ( _strawId< other._strawId);
    }
    bool operator>( const StrawHit other) const{
      return ( _strawId>other._strawId);
    }
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
  };
  inline std::ostream& operator<<( std::ostream& ost,
                                   StrawHit const& hit){
    hit.print(ost,false);
    return ost;
  }

   typedef std::vector<mu2e::StrawHit> StrawHitCollection;

} // namespace mu2e

#endif /* RecoDataProducts_StrawHit_hh */
