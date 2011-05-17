#ifndef CaloCrystalOnlyHit_H
#define CaloCrystalOnlyHit_H

// $Id: CaloCrystalOnlyHit.hh,v 1.2 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author KLG

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "ToyDP/inc/DPIndex.hh"

class art::ProductID;

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

#endif
