#ifndef MCDataProducts_CaloHitMCTruth_hh
#define MCDataProducts_CaloHitMCTruth_hh

// $Id: CaloHitMCTruth.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Ivan Logashenko

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  struct CaloHitMCTruth{

  public:

    CaloHitMCTruth():
      _roId(-1),
      _time(0.),
      _energyDep(0.),
      _charged(0) {
    }

    CaloHitMCTruth( int   roId,
                    float time,
                    float energyDep,
                    int   charged    ):
      _roId(roId),
      _time(time),
      _energyDep(energyDep),
      _charged(charged) {
    }

    // Accessors
    int   id()        const { return _roId; }
    float time()      const { return _time;}
    float energyDep() const { return _energyDep; }
    bool  isCharged() const { return _charged>0; }
    int   charge()    const { return _charged; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    int   _roId;
    float _time;             // (ns)
    float _energyDep;        // (MeV)
    int   _charged;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_CaloHitMCTruth_hh */
