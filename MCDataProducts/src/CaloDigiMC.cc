//
// This is a place to put additional information produced by HitMaker,
//

// Original author Stefano Roberto Soleti
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "MCDataProducts/inc/CaloDigiMC.hh"

using namespace std;

namespace mu2e {
  void          CaloDigiMC::init(){
    _nParticles = 0;
    _totalEDep  = 0.;
    _meanTime   = 0.;
    _timeFirst  = 1e10;
  }

  
  void          CaloDigiMC::addCaloShower(const CaloShowerStep* CaloShower, double HitTimeUnfolded)  {
    
    double    csEdep  = CaloShower->energyMC();
    double    csTime  = HitTimeUnfolded;//CaloShower->time();
    
    int       csSimId = CaloShower->simParticle()->id().asInt();
    
    
    //first: search if the particle is already present
    int      found(0), index(-1);
      
    for (int i=0; i<_nParticles; ++i){
      int   simId = _simParticle.at(i)->id().asInt();
      if (simId == csSimId){
	found = 1;
	index = i;
	break;
      }
    }
      
    if (csTime < _timeFirst){
      _timeFirst = csTime;
    }

    if (found == 1){//update the information
      _time.at(index) = (_time.at(index)*eDep(index) + csEdep*csTime) / (eDep(index) + csEdep);
      _eDep.at(index) += csEdep; 
    }else{   //add particle
	
      _eDep. push_back(csEdep);
      _time. push_back(csTime);
      _pdgId.push_back(CaloShower->simParticle()->pdgId());
	
      _simParticle.push_back(CaloShower->simParticle());
	
      ++_nParticles;
    }
      
    _meanTime     =  (_meanTime * _totalEDep + csEdep*csTime) / (_totalEDep    + csEdep);
    _totalEDep    += csEdep;


  }

} // namespace mu2e
