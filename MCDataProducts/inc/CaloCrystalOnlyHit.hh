#ifndef MCDataProducts_CaloCrystalOnlyHit_hh
#define MCDataProducts_CaloCrystalOnlyHit_hh

// $Id: CaloCrystalOnlyHit.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author KLG

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "DataProducts/inc/DPIndex.hh"

namespace art {
  class ProductID;
}

namespace mu2e {

  class CaloCrystalOnlyHit{

  public:

    CaloCrystalOnlyHit():
      _crystalId(-1),
      _time(0.),
      _energyDep(0.)
    {}

    CaloCrystalOnlyHit(int crystalId,double time, double energy):
      _crystalId(crystalId),
      _time(time),
      _energyDep(energy)
    {}

    // Accessors

    int              id() const { return _crystalId; }
    float            time()      const { return _time; }
    float            energyDep() const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    void setEnergyDep(double energy);

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

   private:

    int              _crystalId;
    float            _time;             // (ns)
    float            _energyDep;        // (MeV)

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalOnlyHit const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_CaloCrystalOnlyHit_hh */
