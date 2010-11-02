// 
// This is a place to put additional information produced by HitMaker,
// like true drift distance, signal propagation time, etc.
//
// $Id: CaloCrystalHitMCTruth.cc,v 1.2 2010/11/02 03:24:39 genser Exp $
// $Author: genser $
// $Date: 2010/11/02 03:24:39 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalHitMCTruth.hh"
#include "ToyDP/inc/CaloHitMCTruth.hh"
#include "ToyDP/inc/DPIndex.hh"

using namespace std;

namespace mu2e {

  CaloCrystalHitMCTruth::CaloCrystalHitMCTruth(int crystalId, edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit) :
    _crystalId(crystalId),
    _time(hit.time()),
    _energyDep(hit.energyDep()),
    _energyDepTotal(hit.energyDep()),
    _numberOfROIdsUsed(1),
    _charge(hit.charge())
  {
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
  }

  // operator += CaloHitMCTruth
  CaloCrystalHitMCTruth& CaloCrystalHitMCTruth::add(edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit) {
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
    _energyDep += hit.energyDep();
    _energyDepTotal += hit.energyDep();
    ++_numberOfROIdsUsed;
    _charge += hit.charge();
    return *this;
  }

  // almost like one of the constructors, plays a role of a two
  // argument assignment operator
  void CaloCrystalHitMCTruth::assign(int crystalId, edm::ProductID const & caloHitCollId, CaloHitMCTruth const & hit) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = hit.energyDep();
    _energyDepTotal = hit.energyDep();
    _roIds.clear();
    _roIds.push_back(DPIndex(caloHitCollId,hit.roId()));
    _numberOfROIdsUsed = 1;
    _charge = hit.charge();
    return;
  }

  // Print the information found in this hit.
  void CaloCrystalHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "  << _energyDep
        << " energyDepT: "  << _energyDepTotal
        << " used roids: "  << _numberOfROIdsUsed
        << " crystal roids:";
    for (size_t i=0; i!=_roIds.size(); ++i) {
      ost  << " " << i << ": " <<_roIds[i];
    }
    ost << " charge: "      << _charge;
 
    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
