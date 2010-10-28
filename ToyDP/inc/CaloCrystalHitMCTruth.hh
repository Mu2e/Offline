#ifndef CaloCrystalHitMCTruth_H
#define CaloCrystalHitMCTruth_H

// $Id: CaloCrystalHitMCTruth.hh,v 1.1 2010/10/28 20:43:58 genser Exp $
// $Author: genser $
// $Date: 2010/10/28 20:43:58 $
//
// Original author KLG

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "ToyDP/inc/DPIndex.hh"

class edm::ProductID;

namespace mu2e { 

  class CaloHit;

  struct CaloCrystalHitMCTruth{

  public:

    CaloCrystalHitMCTruth():
      _crystalId(-1),
      _time(0.),
      _energyDep(0.),
      _charged(0), 
      _numberOfROIdsUsed(0) 
    {}

    // Accessors

    int              crystalId() const { return _crystalId; }
    float            time()      const { return _time; }
    float            energyDep() const { return _energyDep; }
    bool             isCharged() const { return _charged>0; }
    int              numberOfROIdsUsed() const { return _numberOfROIdsUsed; }
    std::vector<DPIndex> const & roIds() const { return _roIds; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
    
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    int              _crystalId;
    float            _time;             // (ns)
    float            _energyDep;        // (MeV)
    int              _charged;
    int              _numberOfROIdsUsed;  // roIds used to calculate _energyDep
    std::vector<DPIndex> _roIds;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif
