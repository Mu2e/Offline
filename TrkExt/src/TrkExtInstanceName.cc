//
//  $Id: TrkExtInstanceName.cc,v 1.1 2012/10/23 16:26:40 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/10/23 16:26:40 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>
#include <vector>

// Framework includes.

#include "TrkExt/inc/TrkExtInstanceName.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "TrkBase/TrkParticle.hh"

// In TrkParticle :
// eMinus, ePlus, muMinus, muPlus, piMunis, piPlus, KMilus, KPlus
// In TrkFitDirection: 
// Upstream=1, Downstream=0
// instance = fdir.name() + tpart.name()

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
/*
  void TrkExtInstanceName::constructAll() {

    TrkParticle::type fitPtl[8] = {
      TrkParticle::e_minus,
      TrkParticle::e_plus,
      TrkParticle::mu_minus,
      TrkParticle::mu_plus,
      TrkParticle::pi_minus,
      TrkParticle::pi_plus,
      TrkParticle::anti_p_minus,
      TrkParticle::p_plus
    };

    for (int j = 0 ; j <8 ; ++j) {
      TrkExtInstanceNameEntry tmp1 (fitPtl[j], true);
      _entries.push_back(tmp1);
      TrkExtInstanceNameEntry tmp2 (fitPtl[j], false);
      _entries.push_back(tmp2);
    }
  }
*/
  void TrkExtInstanceName::construct (TrkParticle::type _hepid, bool _up_down, string _fitterName) 
  {
    TrkExtInstanceNameEntry tmp1 (_hepid, _up_down, _fitterName);
    _entries.push_back(tmp1);
  }



} // end namespace mu2e

