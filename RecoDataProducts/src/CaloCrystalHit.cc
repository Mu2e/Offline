//
// CaloCrystalHit to be created based on CaloHit's
//
// $Id: CaloCrystalHit.cc,v 1.5 2014/05/14 18:14:21 murat Exp $
// $Author: murat $
// $Date: 2014/05/14 18:14:21 $
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

  //constructor that uses the CaloDigi
  CaloCrystalHit::CaloCrystalHit(int  CrystalId, std::vector<RecoCaloDigiPtr>  RecoCaloDigiCol){
    _crystalId = CrystalId;

    //set the time and teh energy as the mean. It assumes that realtive dealy between 
    //the channels have been already accounted in each RacoCaloDigi
    int     nRecoDigis = RecoCaloDigiCol.size();

    const RecoCaloDigi*      recoDigi;

    //initialize internal parameters
    _time      = 0.;
    _energyDep = 0.;

    for (int i=0; i<nRecoDigis; ++i){
      art::Ptr<RecoCaloDigi>  recoDigiPtr  = RecoCaloDigiCol.at(i);
      recoDigi = recoDigiPtr.operator ->();
      _time      += recoDigi->time();  
      _energyDep += recoDigi->edep();  
    }
    _time      /= nRecoDigis;
    _energyDep /= nRecoDigis;

    _energyDepTotal = _energyDep;

    _numberOfROIdsUsed = nRecoDigis;

    _recoCaloDigis = RecoCaloDigiCol;
  }
  

  // operator += CaloHit
  CaloCrystalHit& CaloCrystalHit::add(CaloHit const & hit, CaloHitPtr const& chPtr ) {
    _readouts.push_back(chPtr);

    double he = hit.energyDep();
    double ht = hit.time();

    _time            = (_time*_energyDep+ht*he)/(_energyDep+he);
    _energyDep      += he;
    _energyDepTotal += he;
    ++_numberOfROIdsUsed;
    return *this;
  }

  //-----------------------------------------------------------------------------
  // hit merging: merge 'Hit' with 'this'
  //-----------------------------------------------------------------------------
  void CaloCrystalHit::add(CaloCrystalHit* Hit) {

    double he = Hit->energyDep();
    double ht = Hit->time();

    _time            = (_time*_energyDep+ht*he)/(_energyDep+he);
    _energyDep      += he;
    _energyDepTotal += he;

    for (std::vector<art::Ptr<mu2e::CaloHit>>::const_iterator k = Hit->readouts().begin(); 
	 k != Hit->readouts().end(); k++) {
      _readouts.push_back(*k);
      ++_numberOfROIdsUsed;
    }
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

  void CaloCrystalHit::setTime     (double Time){
    _time      = Time;
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
