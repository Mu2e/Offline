//
//
//  Original author MyeongJae Lee
//
// Note : See the included file for particle/direction definitions. 
//
// In TrkParticle :
// eMinus, ePlus, muMinus, muPlus, piMunis, piPlus, KMilus, KPlus
// are defined.
// In TrkFitDirection: 
// Upstream=1, Downstream=0
// are assigned. 
// instance = fdir.name() + tpart.name()
//

// C++ includes.
#include <iostream>
#include <string>
#include <vector>

// Framework includes.

#include "TrkExt/inc/TrkExtInstanceName.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/TrkBase/TrkParticle.hh"


using namespace std;

namespace mu2e {

  TrkExtInstanceNameEntry::TrkExtInstanceNameEntry (TrkParticle::type _hepid, bool _up_down, std::string _fitterName) :
    fitterName(_fitterName) 
  {
    switch (_hepid) {
      case TrkParticle::e_minus:
        charge = -1;
        mass2 = TrkExtParticleMass::MASS_ELECTRON * TrkExtParticleMass::MASS_ELECTRON;
        break;
      case TrkParticle::e_plus:
        charge = 1;
        mass2 = TrkExtParticleMass::MASS_ELECTRON * TrkExtParticleMass::MASS_ELECTRON;
        break;
      case TrkParticle::mu_minus:
        charge = -1;
        mass2 = TrkExtParticleMass::MASS_MUON * TrkExtParticleMass::MASS_MUON;
        break;
      case TrkParticle::mu_plus:
        charge = 1;
        mass2 = TrkExtParticleMass::MASS_MUON * TrkExtParticleMass::MASS_MUON;
        break;
      case TrkParticle::pi_minus:
        charge = -1;
        mass2 = TrkExtParticleMass::MASS_PION * TrkExtParticleMass::MASS_PION;
        break;
      case TrkParticle::pi_plus:
        charge = 1;
        mass2 = TrkExtParticleMass::MASS_PION * TrkExtParticleMass::MASS_PION;
        break;
      case TrkParticle::p_plus:
        charge = 1;
        mass2 = TrkExtParticleMass::MASS_PROTON * TrkExtParticleMass::MASS_PROTON;
        break;
      case TrkParticle::anti_p_minus:
        charge = -1;
        mass2 = TrkExtParticleMass::MASS_PROTON * TrkExtParticleMass::MASS_PROTON;
        break;
      default:
        cerr << "TrkExtInstanceNameEntry Error: not implemented mode" << endl;
        charge = 0;
        mass2 = 0;
     }
     updown = _up_down;
     hepid = (int)_hepid;
     string fitDirName;
     if (updown) fitDirName = TrkFitDirection(TrkFitDirection::upstream).name();
     else         fitDirName = TrkFitDirection(TrkFitDirection::downstream).name();
     string fitPtlName = TrkParticle(_hepid).name();
 
     name = fitDirName + fitPtlName;
     ntrk = 0;
  } 



  TrkExtInstanceName::TrkExtInstanceName () {
    _entries.clear();
  }
  
  void TrkExtInstanceName::construct (TrkParticle::type _hepid, bool _up_down, string _fitterName) 
  {
    TrkExtInstanceNameEntry tmp1 (_hepid, _up_down, _fitterName);
    _entries.push_back(tmp1);
  }



} // end namespace mu2e

