//
// $Id: ConvElecUtilities.cc,v 1.2 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Gianni Onorato
//


// Framework includes
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Mu2eUtilities/inc/ConvElecUtilities.hh"
#include <iostream>

using namespace std;

namespace mu2e {

  ConvElecUtilities::ConvElecUtilities(const art::Event & event,
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
    _nconv = 0;
    checkConvElec(event);

    if (_nconv>1) {
      throw cet::exception("RANGE")
        << "More than 1 conversion electron in the event";
    }
    lookAtHits(event);
    lookAtCalo(event);
  }

  ConvElecUtilities::~ConvElecUtilities()
  {
  }


  //Loop within the Simulated Particles of the event and find the
  //generated conversion electron.
  //This method is called by the constructor

  void ConvElecUtilities::checkConvElec(const art::Event & event) {

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
          if (_genParticles->at(sim.generatorIndex()).generatorId() == GenId::conversionGun) {
            _convTrackId = sim.id();
            _simParticle = unique_ptr<SimParticle>( new SimParticle(sim) );
            //  cout << n << " and " << _convTrackId << endl;
            _nconv++;
            //  cout << "conv track id is " << _convTrackId << endl;
          }
        }
      }
    } else { cout << "No SimParticles in the event" << endl; }
  }

  const SimParticle& ConvElecUtilities::simConvElec() const{
    return *_simParticle;
  }

  const GenParticle& ConvElecUtilities::genConvElec() {
    return _genParticles->at(_simParticles->at(_convTrackId).generatorIndex());
  }

  //Loop through the hits. Identify hits made by conv electron
  //and earliest hit
  //This method is called by the constructor

  void ConvElecUtilities::lookAtHits(const art::Event & event) {
    event.getByLabel(_g4ModuleLabel, _trackerStepPoints, hits);
    double time = 10e15; // dummy value
    _earliestidx = 10000; // dummy value
    for (size_t i=0; i<hits->size(); ++i) {
      StepPointMC const& hit = (*hits)[i];
      if ( hit.trackId() == _convTrackId ) {
        if (hit.time()<time) {
          time = hit.time();
          _earliestidx = i;
        }
        _convElecHits.push_back(i);
        //        cout << "Hit straw index " << hit.strawIndex()
        //     << " and time " << hit.time()
        //     << "earliest is " << time << endl;
        vector<StrawIndex>::iterator it = _convElecStrawIdx.begin();
        bool itexists = false;
        while (it != _convElecStrawIdx.end()) {
          if (*it == hit.strawIndex()) {
            itexists = true;
          } 
          ++it;
        }
        if (!itexists) {
          _convElecStrawIdx.push_back(hit.strawIndex());
        }
        _totEDep += hit.totalEDep();
      }
    }
  }

  void ConvElecUtilities::lookAtCalo(const art::Event & event) {

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
	//	     << " while the conversion electron has track id " << _convTrackId << endl;
	if ( hit.trackId() == _convTrackId ) {
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

  bool ConvElecUtilities::gotCaloHit() {
    return _stepincalo;
  }

  StepPointMC const& ConvElecUtilities::firstCaloHit() {
    if (_stepincalo) {
      return (const StepPointMC&)_earliestSPMC; 
    }  else {
    //      return (*calohits)[_earliestcry];
      throw cet::exception("RANGE")
        << "No hit associated to Conversion Electron track in calorimeter.";
    }
  }
  
  //Return a reference to the easliest hit in the calorimeter
  //Return a reference to the earliest ConvElectron Hit
  StepPointMC const& ConvElecUtilities::firstHit() {
    if (hasStepPointMC()) {
      return (*hits)[_earliestidx];
    } else {
      throw cet::exception("RANGE")
        << "No hit associated to Conversion Electron track.";
    }
  }

  //return the strawindex of the earliest hit of the conversion electron
  StrawIndex ConvElecUtilities::earliestStrawIndex() const {
    if (hasStepPointMC()) {
    return (*hits)[_earliestidx].strawIndex();
    } else {
      throw cet::exception("RANGE")
        << "No hit associated to Conversion Electron track.";
    }
  }




} // end namespace mu2e
