//
// $Id: BkgElecUtilities.cc,v 1.2 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Gianni Onorato
//


// Framework includes
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Mu2eUtilities/inc/BkgElecUtilities.hh"
#include <iostream>

using namespace std;

namespace mu2e {

  BkgElecUtilities::BkgElecUtilities(const art::Event & event,
                                       string const &generatorModuleLabel,
                                       string const &g4ModuleLabel,
                                       string const &trackerStepPoints,
                                       string const &caloROModuleLabel):
    _generatorModuleLabel( generatorModuleLabel ),
    _g4ModuleLabel( g4ModuleLabel ),
    _trackerStepPoints( trackerStepPoints ),
    _caloROlabel( caloROModuleLabel ),
    _totEDep(0),
    _stepincalo(false)
  {
    _nbkg = 0;
    checkBkgElec(event);

    if (_nbkg>2) {
      throw cet::exception("RANGE")
        << "More than 2 electrons in the event";
    }
    lookAtHits(event);
    lookAtCalo(event);
  }

  BkgElecUtilities::~BkgElecUtilities()
  {
  }


  //Loop within the Simulated Particles of the event and find the
  //generated bkg electron.
  //This method is called by the constructor

  void BkgElecUtilities::checkBkgElec(const art::Event & event) {

    event.getByLabel(_generatorModuleLabel, _genParticles);
    event.getByLabel(_g4ModuleLabel, _simParticles);
    if (_simParticles.isValid() && !_simParticles->empty()) {

      int n(0);
      for ( SimParticleCollection::const_iterator i=_simParticles->begin();
            i!=_simParticles->end(); ++i ){

        SimParticle const& sim = i->second;
        ++n;

        //        cout << "Id = " << sim.id() << "\t" << "PID = "
        //     << sim.pdgId() << " from "
        //     << (sim.fromGenerator()?"generator":"")
        //     << (sim.madeInG4()?"G4":"") << endl;
        if ( sim.fromGenerator() ) {
          //   cout << "index " << sim.generatorIndex() << " and "
          //      << _genParticles->at(sim.generatorIndex()).generatorId().name() << endl;
          if ( _genParticles->at(sim.generatorIndex()).generatorId().id() != GenId::conversionGun){
            _bkgTrackId = sim.id();
            _simParticle = unique_ptr<SimParticle>( new SimParticle(sim) );
            //  cout << n << " and " << _bkgTrackId << endl;
            _nbkg++;
            //  cout << "bkg track id is " << _bkgTrackId << endl;
          }
        }
      }
    } else { cout << "No SimParticles in the event" << endl; }
  }

  const SimParticle& BkgElecUtilities::simBkgElec() const{
    return *_simParticle;
  }

  const GenParticle& BkgElecUtilities::genBkgElec() {
    return _genParticles->at(_simParticles->at(_bkgTrackId).generatorIndex());
  }

  //Loop through the hits. Identify hits made by bkg electron
  //and earliest hit
  //This method is called by the constructor

  void BkgElecUtilities::lookAtHits(const art::Event & event) {
    event.getByLabel(_g4ModuleLabel, _trackerStepPoints, hits);
    double time = 10e15; // dummy value
    _earliestidx = 10000; // dummy value
    for (size_t i=0; i<hits->size(); ++i) {
      StepPointMC const& hit = (*hits)[i];
      if ( hit.trackId() == _bkgTrackId ) {
        if (hit.time()<time) {
          time = hit.time();
          _earliestidx = i;
        }
        _bkgElecHits.push_back(i);
        //        cout << "Hit straw index " << hit.strawIndex()
        //     << " and time " << hit.time()
        //     << "earliest is " << time << endl;
        vector<StrawIndex>::iterator it = _bkgElecStrawIdx.begin();
        bool itexists = false;
        while (it != _bkgElecStrawIdx.end()) {
          if (*it == hit.strawIndex()) {
            itexists = true;
          } 
          ++it;
        }
        if (!itexists) {
          _bkgElecStrawIdx.push_back(hit.strawIndex());
        }
        _totEDep += hit.totalEDep();
      }
    }
  }

  void BkgElecUtilities::lookAtCalo(const art::Event & event) {

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    //    art::Handle<StepPointMCCollection> steps;
    
    event.getByLabel(_caloROlabel,"CaloHitMCCrystalPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const * collMCPTR = mcptrHandle.product();    
    double time = 10e15; // dummy value
    _earliestcry = 10000; // dummy value
    _earliestvector = 10000;
    for (size_t i1 = 0; i1 < collMCPTR->size(); ++i1) {
      PtrStepPointMCVector::const_iterator i2 = collMCPTR->at(i1).begin();
      while (i2 != collMCPTR->at(i1).end()) {
	//    event.getByLabel(_g4ModuleLabel,"calorimeterRO",calohits);
	//	cout << "Got hits in the calo and they are " << calohits->size() << endl;
	//	cout << "Hit n. " << i+1 << ":" << endl;
	StepPointMC const& hit = **i2;
	//	cout << "Track Id of this hit is " << hit.trackId() 
	//	     << " while the bkg electron has track id " << _bkgTrackId << endl;
	if ( hit.trackId() == _bkgTrackId ) {
	  //	  cout << "Time of this hit is " << hit.time() << endl;
	  if (hit.time()<time) {
	    time = hit.time();
	    _earliestcry = i2->key();
	    _earliestvector = i1;
	    _earliestSPMC = const_cast<StepPointMC&>(hit);
	  }
	}
      	++i2;
      }
    }
    if (_earliestcry != 10000) {
      _stepincalo = true;
    }

  

  }

  bool BkgElecUtilities::gotCaloHit() {
    return _stepincalo;
  }

  StepPointMC const& BkgElecUtilities::firstCaloHit() {
    if (_stepincalo) {
      return (const StepPointMC&)_earliestSPMC; 
    }  else {
    //      return (*calohits)[_earliestcry];
      throw cet::exception("RANGE")
        << "No hit associated to Bkg Electron track in calorimeter.";
    }
  }
  
  //Return a reference to the easliest hit in the calorimeter
  //Return a reference to the earliest BkgElectron Hit
  StepPointMC const& BkgElecUtilities::firstHit() {
    if (hasStepPointMC()) {
      return (*hits)[_earliestidx];
    } else {
      throw cet::exception("RANGE")
        << "No hit associated to Bkg Electron track.";
    }
  }

  //return the strawindex of the earliest hit of the bkg electron
  StrawIndex BkgElecUtilities::earliestStrawIndex() const {
    if (hasStepPointMC()) {
    return (*hits)[_earliestidx].strawIndex();
    } else {
      throw cet::exception("RANGE")
        << "No hit associated to Bkg Electron track.";
    }
  }




} // end namespace mu2e
