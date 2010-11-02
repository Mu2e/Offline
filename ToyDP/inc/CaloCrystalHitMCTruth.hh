#ifndef CaloCrystalHitMCTruth_H
#define CaloCrystalHitMCTruth_H

// $Id: CaloCrystalHitMCTruth.hh,v 1.2 2010/11/02 03:24:23 genser Exp $
// $Author: genser $
// $Date: 2010/11/02 03:24:23 $
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
      _energyDep(0.),
      _energyDepTotal(0.),
      _numberOfROIdsUsed(0),
      _charge(0)
    {}

    CaloCrystalHitMCTruth(int crystalId, edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit);
    // Accessors

    int              Id() const { return _crystalId; }
    float            time()      const { return _time; }
    float            energyDep() const { return _energyDep; }
    float            energyDepTotal() const { return _energyDepTotal; }
    bool             isCharged() const { return _charge>0; }
    int              charge()    const { return _charge; }
    int              numberOfROIdsUsed() const { return _numberOfROIdsUsed; }
    std::vector<DPIndex> const & roIds() const { return _roIds; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
    
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // almost operator += CaloHitMCTruth
    CaloCrystalHitMCTruth& add(edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit);

    // almost like one of the constructors, plays a role of a two
    // argument assignment operator
    void assign(int crystalId, edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit);

   private:

    int              _crystalId;
    float            _time;             // (ns)
    float            _energyDep;        // (MeV)
    float            _energyDepTotal;     // (MeV)
    int              _numberOfROIdsUsed;  // roIds used to calculate _energyDep
    std::vector<DPIndex> _roIds;
    int              _charge;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif
