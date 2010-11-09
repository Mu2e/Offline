//
// $Id: ConvElecUtilities.cc,v 1.4 2010/11/09 21:06:04 onoratog Exp $
// $Author: onoratog $
// $Date: 2010/11/09 21:06:04 $
//
// Original author Gianni Onorato
//


// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"

// Mu2e includes
#include "Mu2eG4/inc/ConvElecUtilities.hh"



using namespace std;

namespace mu2e {

  ConvElecUtilities::ConvElecUtilities(const edm::Event & event,
                                       string g4ModuleLabel, 
                                       string trackerStepPoints):
    _g4ModuleLabel( g4ModuleLabel ),
    _trackerStepPoints( trackerStepPoints )
  {
    _nconv = 0;
    checkConvElec(event);
    if (_nconv>1) {
      throw cms::Exception("RANGE")
        << "More than 1 conversion electron in the event";
    }
    lookAtHits(event);
  }

  ConvElecUtilities::~ConvElecUtilities()
  {
  }


  //Loop within the Simulated Particles of the event and find the 
  //generated conversion electron.
  //This method is called by the constructor

  void ConvElecUtilities::checkConvElec(const edm::Event & event) {

    event.getByType(_genParticles);
    event.getByType(_simParticles);
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
          //      << _genParticles->at(sim.generatorIndex())._generatorId.name() << endl;
          if (_genParticles->at(sim.generatorIndex())._generatorId == GenId::conversionGun) {
            _convTrackId = sim.id();
            _simParticle = auto_ptr<SimParticle>( new SimParticle(sim) );
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

  const ToyGenParticle& ConvElecUtilities::genConvElec() {
    return _genParticles->at(_simParticles->at(_convTrackId).generatorIndex());
  }

  //Loop through the hits. Identify hits made by conv electron
  //and earliest hit
  //This method is called by the constructor

  void ConvElecUtilities::lookAtHits(const edm::Event & event) {
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
        _convElecStrawIdx.push_back(hit.strawIndex());
      }
    }
  }

  //Return a reference to the earliest ConvElectron Hit
  StepPointMC const& ConvElecUtilities::firstHit() {
    if (hasStepPointMC()) {
      return (*hits)[_earliestidx]; 
    } else { 
      throw cms::Exception("RANGE")
        << "No hit associated to Conversion Electron track.";
    }
  }

  //return the strawindex of the earliest hit of the conversion electron
  StrawIndex ConvElecUtilities::earliestStrawIndex() const {
    if (hasStepPointMC()) {
    return (*hits)[_earliestidx].strawIndex();
    } else { 
      throw cms::Exception("RANGE")
        << "No hit associated to Conversion Electron track.";
    }
  }
  
  
} // end namespace mu2e
