//
// CaloCrystalHit to be created based on CaloHit's
//
// $Id: CaloCrystalHit.cc,v 1.3 2012/11/17 00:06:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/17 00:06:25 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "art/Persistency/Provenance/ProductID.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloHit.hh"

using namespace std;

namespace mu2e {

  CaloCrystalHit::CaloCrystalHit(int crystalId, CaloHit const & hit, CaloHitPtr const& chPtr ) :
    _crystalId(crystalId),
    _time(hit.time()),
    _energyDep(hit.energyDep()),
    _energyDepTotal(hit.energyDep()),
    _numberOfROIdsUsed(1),
    _readouts(1,CaloHitPtr(chPtr))
  {}

  // operator += CaloHit
  CaloCrystalHit& CaloCrystalHit::add(CaloHit const & hit, CaloHitPtr const& chPtr ) {
    _readouts.push_back(chPtr);
    _energyDep += hit.energyDep();
    _energyDepTotal += hit.energyDep();
    ++_numberOfROIdsUsed;
    return *this;
  }

  CaloCrystalHit& CaloCrystalHit::addEnergyToTot( CaloHit const & hit) {
    _energyDepTotal += hit.energyDep();
    return *this;
  }

  // almost like one of the constructors, plays a role of a two
  // argument assignment operator
  void CaloCrystalHit::assign(int crystalId, CaloHit const & hit, CaloHitPtr const& chPtr ) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = hit.energyDep();
    _energyDepTotal = hit.energyDep();
    _readouts.clear();
    _readouts.push_back(chPtr);
    _numberOfROIdsUsed = 1;
    return;
  }

  void CaloCrystalHit::assignEnergyToTot(int crystalId, CaloHit const & hit) {
    _crystalId = crystalId;
    _time = hit.time();
    _energyDep = 0.0;
    _energyDepTotal = hit.energyDep();
    _readouts.clear();
    _numberOfROIdsUsed = 0;
    return;
  }

  void CaloCrystalHit::setEnergyDep(double energy) {

      _energyDep = energy;

    return;

  }
  void CaloCrystalHit::setEnergyDepTotal(double energy) {

      _energyDepTotal = energy;

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
    for (size_t i=0; i!=_readouts.size(); ++i) {
      ost  << " " << i << ": " << _readouts[i];
    }


    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
