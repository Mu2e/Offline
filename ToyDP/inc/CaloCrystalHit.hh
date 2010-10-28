#ifndef CaloCrystalHit_H
#define CaloCrystalHit_H

// $Id: CaloCrystalHit.hh,v 1.1 2010/10/28 20:43:58 genser Exp $
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

  struct CaloCrystalHit{

  public:

    CaloCrystalHit():
      _crystalId(-1),
      _time(0.),
      _energyDep(0.), 
      _numberOfROIdsUsed(0)
    {}

    CaloCrystalHit(int  crystalId, edm::ProductID const & caloHitCollId, CaloHit const & hit);

    // Accessors

    int              crystalId() const { return _crystalId; }
    float            time()      const { return _time;      }
    float            energyDep() const { return _energyDep; }
    int              numberOfROIdsUsed() const { return _numberOfROIdsUsed; }
    std::vector<DPIndex> const & roIds() const { return _roIds; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
    
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // almost operator += CaloHit
    CaloCrystalHit& add(edm::ProductID const & caloHitCollId, CaloHit const & hit);

    // almost like one of the constructors, plays a role of a two
    // argument assignment operator
    void assign(int crystalId, edm::ProductID const & caloHitCollId, CaloHit const & hit);

  private:

    //   iv*ncrys + ic; // Crystal ID

    int              _crystalId;
    float            _time;               // (ns)
    float            _energyDep;          // (MeV)
    int              _numberOfROIdsUsed;  // roIds used to calculate _energyDep
    std::vector<DPIndex> _roIds;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalHit const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif
