#include "RecoDataProducts/inc/StrawHitFlag.hh"
#ifndef RecoDataProducts_StrawHit_hh
#define RecoDataProducts_StrawHit_hh
//
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.hh,v 1.2 2012/03/01 19:30:06 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/03/01 19:30:06 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"

namespace mu2e {

  struct StrawHit{
  private:

    StrawIndex       _strawIndex;       // See note 1.
    float            _time;             // (ns)
    float            _dt;               // (ns)
    float            _energyDep;        // (MeV)

  public:

    StrawHit():
      _strawIndex(StrawIndex(-1)),
      _time(0.),
      _dt(0.),
      _energyDep(0.) {
    }

    // Constructor for a hit that came from an unpacked digi, either
    // from data or from the full MC chain.
    StrawHit( StrawIndex       strawIndex,
              float            time,
              float            dt,
              float            energyDep  ):
      _strawIndex(strawIndex),
      _time(time),
      _dt(dt),
      _energyDep(energyDep) {
    }

    // Accessors
    StrawIndex strawIndex() const { return _strawIndex; }
    float      time()       const { return _time;}
    float      dt()         const { return _dt;}
    float      energyDep()  const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
        bool operator==(StrawHit const& other) const {
      return (_strawIndex==other._strawIndex&&
          _time==other._time&&
          _dt==other._dt&&
              _energyDep==other._energyDep);
    }
    bool operator<( const StrawHit other) const{
      return ( _strawIndex< other._strawIndex);
    }
    bool operator>( const StrawHit other) const{
      return ( _strawIndex>other._strawIndex);
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
