//
// CaloCrystalHit to be created based on CaloHit's
//
// $Id: CaloCrystalHit.cc,v 1.6 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//

// C++ includes
#include <ostream>

// Framework includes.
//#include "art/Utilities/Exception.h"
#include "art/Persistency/Provenance/ProductID.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalHit.hh"
#include "ToyDP/inc/CaloHit.hh"
#include "ToyDP/inc/DPIndex.hh"

using namespace std;

namespace mu2e {

  CaloCrystalHit::CaloCrystalHit(int crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit) :
    _crystalId(crystalId),
    _time(hit.time()),
    _energyDep(hit.energyDep()),
    _energyDepTotal(hit.energyDep()),
    _numberOfROIdsUsed(1),
    _roIds(1,DPIndex(caloHitCollId,hit.id()))
  {}

  // operator += CaloHit
  CaloCrystalHit& CaloCrystalHit::add(art::ProductID const & caloHitCollId, CaloHit const & hit) {
    _roIds.push_back(DPIndex(caloHitCollId,hit.id()));
    _energyDep += hit.energyDep();
    _energyDepTotal += hit.energyDep();
    ++_numberOfROIdsUsed;
    return *this;
  }

  CaloCrystalHit& CaloCrystalHit::addEnergyToTot(art::ProductID const & caloHitCollId, CaloHit const & hit) {
    _energyDepTotal += hit.energyDep();
    return *this;
  }

  // almost like one of the constructors, plays a role of a two
  // argument assignment operator
  void CaloCrystalHit::assign(int crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = hit.energyDep();
    _energyDepTotal = hit.energyDep();
    _roIds.clear();
    _roIds.push_back(DPIndex(caloHitCollId,hit.id()));
    _numberOfROIdsUsed = 1;
    return;
  }

  void CaloCrystalHit::assignEnergyToTot(int crystalId, art::ProductID const & caloHitCollId, CaloHit const & hit) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = 0.0;
    _energyDepTotal = hit.energyDep();
    _roIds.clear();
    _numberOfROIdsUsed = 0;
    return;
  }

  void CaloCrystalHit::setEnergyDep(double energy) {

      _energyDep = energy;

    return;

  }



  // Print the information found in this hit.
  void CaloCrystalHit::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit:   "
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "   << _energyDep
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
