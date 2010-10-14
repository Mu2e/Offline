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
      for (size_t i=0; i<_simParticles->size(); ++i) {
        SimParticle const& sim = _simParticles->at(i);
        //        cout << "Id = " << sim.id() << "\t" << "PID = " 
        //     << sim.pdgId() << " from " 
        //     << (sim.fromGenerator()?"generator":"")
        //     << (sim.madeInG4()?"G4":"") << endl;
        if ( sim.fromGenerator() ) {
          //   cout << "index " << sim.generatorIndex() << " and " 
          //      << _genParticles->at(sim.generatorIndex())._generatorId.name() << endl;
          if (_genParticles->at(sim.generatorIndex())._generatorId.name() == "conversionGun") {
            _convTrackId = sim.id();
            _nconv++;
            //  cout << "conv track id is " << _convTrackId << endl;
          }
        }
      }
    } else { cout << "No SimParticles in the event" << endl; }
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
      //      cout << "Hit " << i+1 << " : " << hit.trackId() 
      //     << "   " << hit.volumeId() << endl;
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

  //Returns how many hits come from the convElectron
  size_t ConvElecUtilities::hasStepPointMC() const {
    return _convElecHits.size();
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
    return (*hits)[_earliestidx].strawIndex();
  }


  //Number of Conversion Electrons. Used as a check
  int ConvElecUtilities::nConvElec() const {
    return _nconv;
  }


  //return a vector of index related to the stepPointMCCollection
  //identifying hits of the conversion electron
  vector<size_t> ConvElecUtilities::convElecHitsIdx() {
    return _convElecHits;
  } //maybe it could be transformed in a vector of references to event hits 



  //return a vector of StrawIndex related to hits of the conversion electron
  vector<StrawIndex> ConvElecUtilities::convElecStrawIdx() {
    return _convElecStrawIdx;
  }


} // end namespace mu2e
