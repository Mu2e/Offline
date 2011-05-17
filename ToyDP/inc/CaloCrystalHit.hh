#ifndef ToyDP_CaloCrystalHit_hh
#define ToyDP_CaloCrystalHit_hh

// $Id: CaloCrystalHit.hh,v 1.6 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author KLG

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes
#include "ToyDP/inc/DPIndex.hh"

class art::ProductID;

namespace mu2e { 

  class CaloHit;

  class CaloCrystalHit{

  public:

    CaloCrystalHit():
      _crystalId(-1),
      _time(0.),
      _energyDep(0.), 
      _energyDepTotal(0.), 
      _numberOfROIdsUsed(0)
    {}

    CaloCrystalHit(int  crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit);

    // Accessors

    int              id()                const { return _crystalId; }
    float            time()              const { return _time;      }
    float            energyDep()         const { return _energyDep; }
    float            energyDepTotal()    const { return _energyDepTotal; }
    int              numberOfROIdsUsed() const { return _numberOfROIdsUsed; }
    std::vector<DPIndex> const & roIds() const { return _roIds; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.
    
    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // almost operator += CaloHit
    CaloCrystalHit& add(art::ProductID const & caloHitCollId, CaloHit const & hit);

    CaloCrystalHit& addEnergyToTot(art::ProductID const & caloHitCollId, CaloHit const & hit);

    // almost like one of the constructors, plays a role of a two
    // argument assignment operator
    void assign(int crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit);

    void assignEnergyToTot(int crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit);

    void setEnergyDep(double energy);

  private:

    //   iv*ncrys + ic; // Crystal ID

    int              _crystalId;
    float            _time;               // (ns)
    float            _energyDep;          // (MeV)
    float            _energyDepTotal;     // (MeV)
    int              _numberOfROIdsUsed;  // roIds used to calculate _energyDep
    std::vector<DPIndex> _roIds;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CaloCrystalHit const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* ToyDP_CaloCrystalHit_hh */
