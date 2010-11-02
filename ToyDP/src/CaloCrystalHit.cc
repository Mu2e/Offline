// 
// CaloCrystalHit to be created based on CaloHit's 
//
// $Id: CaloCrystalHit.cc,v 1.2 2010/11/02 03:19:50 genser Exp $
// $Author: genser $
// $Date: 2010/11/02 03:19:50 $
//

// C++ includes
#include <ostream>

// Framework includes.
//#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Provenance/interface/ProductID.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalHit.hh"
#include "ToyDP/inc/CaloHit.hh"
#include "ToyDP/inc/DPIndex.hh"

using namespace std;

namespace mu2e {

  CaloCrystalHit::CaloCrystalHit(int crystalId, edm::ProductID const & caloHitCollId, CaloHit const & hit) :
    _crystalId(crystalId),
    _time(hit.time()),
    _energyDep(hit.energyDep()),
    _energyDepTotal(hit.energyDep()),
    _numberOfROIdsUsed(1)
  {
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
  }

  // operator += CaloHit
  CaloCrystalHit& CaloCrystalHit::add(edm::ProductID const & caloHitCollId, CaloHit const & hit) {
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
    _energyDep += hit.energyDep();
    _energyDepTotal += hit.energyDep();
    ++_numberOfROIdsUsed;
    return *this;
  }

  CaloCrystalHit& CaloCrystalHit::addEnergyToTot(edm::ProductID const & caloHitCollId, CaloHit const & hit) {
    _energyDepTotal += hit.energyDep();
    return *this;
  }

  // almost like one of the constructors, plays a role of a two
  // argument assignment operator
  void CaloCrystalHit::assign(int crystalId, edm::ProductID const & caloHitCollId, CaloHit const & hit) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = hit.energyDep();
    _energyDepTotal = hit.energyDep();
    _roIds.clear();
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
    _numberOfROIdsUsed = 1;
    return;
  }

  // Print the information found in this hit.
  void CaloCrystalHit::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit:   "
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "  << _energyDep
        << " energyDepT: "  << _energyDepTotal
        << " used roids: "  << _numberOfROIdsUsed
        << " crystal roids:";
    for (size_t i=0; i!=_roIds.size(); ++i) {
      ost  << " " << i << ": " <<_roIds[i];
    }


    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
