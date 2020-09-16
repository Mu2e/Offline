#ifndef DataProducts_CaloHit_hh
#define DataProducts_CaloHit_hh

//
// Original author Ivan Logashenko

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  struct CaloHit{

  public:

    CaloHit():
      _roId(-1),
      _time(0.),
      _energyDep(0.) {
    }

    CaloHit( int   roId,
             float time,
             float energyDep  ):
      _roId(roId),
      _time(time),
      _energyDep(energyDep) {
    }

    // Accessors
    int        id()         const { return _roId; }
    float      time()       const { return _time;}
    float      energyDep()  const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    int       _roId;
    float     _time;             // (ns)
    float     _energyDep;        // (MeV)

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloHit const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* DataProducts_CaloHit_hh */
