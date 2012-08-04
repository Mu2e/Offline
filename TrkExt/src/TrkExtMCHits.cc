//
//
//  $Id: TrkExtMCHits.cc,v 1.1 2012/08/04 00:22:09 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:22:09 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtMCHits.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"


using namespace CLHEP;
using namespace std;

namespace mu2e {

 
  TrkExtMCHits::TrkExtMCHits( ) { 
    _g4ModuleLabel = "";
    _instanceName = "";
    _nhits = 0;
    _hitcol.clear();
    setEnergyThreshold(0);
    setClusterHitDistance(0);
  }

  TrkExtMCHits::TrkExtMCHits(std::string g4ModuleLabel, std::string instanceName, double eth, double dist) {
    _g4ModuleLabel = g4ModuleLabel;
    _instanceName = instanceName;
    _nhits = 0;
    _hitcol.clear();
    setEnergyThreshold(eth);
    setClusterHitDistance(dist);
  }

  TrkExtMCHits::~TrkExtMCHits() { }

  unsigned int TrkExtMCHits::readMCHits(const art::Event& event) {

    art::Handle<StepPointMCCollection> pahits;
    event.getByLabel(_g4ModuleLabel, _instanceName, pahits);
    if (pahits->size() <=0) {
      _hitcol.clear();
      _nhits = 0;
      return 0;
    }

    unsigned int i, j = -1;
    _nhits = pahits->size();

    for ( i=0; i<pahits->size(); ++i ){
      const StepPointMC& pahit = (*pahits)[i];
      if (pahit.momentum().mag() > _eth && (pahit.simParticle())->isPrimary()) {
        ++j;
        if(j==0) {
          StepPointMCCollection mccol;
          mccol.push_back(pahit);
          _hitcol.push_back(mccol);
        }
        else {
          bool flag = false;
          for (unsigned int k = 0 ; k < _hitcol.size() ; ++k) {
            StepPointMCCollection & mccol = _hitcol[k];
            for (unsigned int l = 0 ; l < mccol.size() ; ++l) {
              StepPointMC & mc = mccol[l];
              if ((mc.position() - pahit.position()).mag() <_dist) {
                mccol.push_back(pahit);
                flag = true;
                break;
              }
            }
            if (flag) break;
          }
          if (!flag) {
            StepPointMCCollection mccol;
            mccol.push_back(pahit);
            _hitcol.push_back(mccol);
          }
        }
      }
    }

    return _hitcol.size();
  }


} // end namespace mu2e

