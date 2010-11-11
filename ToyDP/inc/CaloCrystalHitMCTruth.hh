#ifndef CaloCrystalHitMCTruth_H
#define CaloCrystalHitMCTruth_H

// $Id: CaloCrystalHitMCTruth.hh,v 1.3 2010/11/11 21:16:53 genser Exp $
// $Author: genser $
// $Date: 2010/11/11 21:16:53 $
//
// Original author KLG

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "ToyDP/inc/DPIndex.hh"

class edm::ProductID;

namespace mu2e { 

  class CaloHitMCTruth;

  class CaloCrystalHitMCTruth{

  public:

    CaloCrystalHitMCTruth():
      _crystalId(-1),
      _time(0.),
      _energyDep(0.)
    {}

    CaloCrystalHitMCTruth(int crystalId,double time, double energy):
      _crystalId(crystalId),
      _time(time),
      _energyDep(energy)
    {}

    // Accessors

    int              id() const { return _crystalId; }
    float            time()      const { return _time; }
    float            energyDep() const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
    
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

   private:

    int              _crystalId;
    float            _time;             // (ns)
    float            _energyDep;        // (MeV)

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif
