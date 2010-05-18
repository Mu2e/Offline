//
// Maintain multiple independent random number engines,
// including save and restore of seeds and state. 
//
// $Id: RandomNumberService.cc,v 1.4 2010/05/18 21:16:43 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 21:16:43 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The CMS code on which is this modelled is available at
//    http://cmslxr.fnal.gov/lxr/source/IOMC/RandomEngine/src/RandomNumberGeneratorService.cc
// 
// 2) CLHEP specifies that state will be returned as vector<unsigned long>.  
//    The size of a long is machine dependent. If unsigned long is an 8 byte 
//    variable, only the least significant 4 bytes are filled and the most 
//    significant 4 bytes are zero.  We need to store the state as
//    with a machine indpendent size, which we choose to be uint32_t.
//    This conversion really belongs in the RandomEngineState class but
//    we are constrained by the HepRandomGenerator interface from the framework.
//
// 3) We are currently linking to clhep v1.9.3.2. This is out of date and does not
//    contain the correct code for saving state on 64 bit machines.  So I call
//    getState() instead of get() to restore the state.  Change back when we switch
//    to the correct version of CLHEP.  The issue is masking out the high order word
//    returned from crc32ul on a 64 bit machine.
// 
// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e include files
#include "RandomNumberService/inc/RandomNumberService.hh"
#include "ToyDP/inc/RandomEngineState.hh"
#include "GeneralUtilities/inc/vectorTransform.hh"

// CLHEP include files.
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/engineIDulong.h"

using namespace std;

namespace mu2e {

  RandomNumberService::RandomNumberService(edm::ParameterSet const& iPS, 
                                           edm::ActivityRegistry&iRegistry) :
    _globalSeed(iPS.getUntrackedParameter<int>("globalSeed",-1)),
    _restoreStateLabel(iPS.getUntrackedParameter<string>("restoreStateLabel","")),
    _labels(),
    _states(),
    _seeds(),
    _nPrint(iPS.getUntrackedParameter<int>("nPrint",1)),
    _debug(iPS.getUntrackedParameter<bool>("debug",false))
  {

    // Register for callbacks.
    iRegistry.watchPostBeginJob(this, &RandomNumberService::postBeginJob);
    iRegistry.watchPostEndJob(this, &RandomNumberService::postEndJob);

    iRegistry.watchPreSource(this, &RandomNumberService::preSource);
    iRegistry.watchPostSource(this, &RandomNumberService::postSource);

    iRegistry.watchPreProcessEvent(this, &RandomNumberService::preProcessEvent);
    iRegistry.watchPostProcessEvent(this, &RandomNumberService::postProcessEvent);
  }
  
  RandomNumberService::~RandomNumberService(){
  }
  
  // This is called after all modules have had their beginJob methods called
  // but before the event loop starts.  I would prefer a preBeginJob
  // method but the corresponding watch method does not exist.
  void RandomNumberService::postBeginJob(){

    // The HepRandom default random number generator is HepJamesRandom.
    // A comment in the sources says that "seed" should be within 
    // the range [0,900000000]
    if ( _globalSeed > -1 && _globalSeed < 90000000 ){
      CLHEP::HepRandom::setTheSeed(_globalSeed);
      if ( _debug ) {
        edm::LogInfo("RANDOM") 
          << "Setting the seed of the HepRandom at job start: " << _globalSeed;
      }
    } else if ( _globalSeed != -1 ) {

      // Or do we want to just give a warning message?
      throw cms::Exception("RANGE")
        << "Seed for the HepRandom generator is out of bounds: "
        << _globalSeed;
    }

    // Add the HepRandom singleton engine to the list of managed engines.
    _labels.push_back("@HepRandom");

    vector<uint32_t> seed;
    seed.push_back(_globalSeed);
    _seeds.push_back(seed);

    // Empty state vector as a place holder.
    vector<uint32_t> state;
    _states.push_back(state);
    
  }

  void RandomNumberService::preSource(){
    //cerr << "xxxxpreSourceEvent: " << endl;
  }
  void RandomNumberService::postSource(){
    //cerr << "xxxxpostSourceEvent: " << endl;
  }
  void RandomNumberService::preProcessEvent(edm::EventID const& id, edm::Timestamp const& iTime){
    //cerr << "xxxxpreProcesEvent: " << id << " " << iTime.value() << endl;
  }
  void RandomNumberService::postProcessEvent(edm::Event const& event, edm::EventSetup const&){
    //cerr << "xxxxpostProcessEvent: " << event.id() << endl;
  }

  // A candidate for saving state to a file?
  void RandomNumberService::postEndJob(){
  }

  CLHEP::HepRandomEngine& RandomNumberService::getEngine() const{

    // First version just does this singleton engine.
    return *CLHEP::HepRandom::getTheEngine();
  }

  uint32_t RandomNumberService::mySeed() const {
    throw cms::Exception("OBSOLETE")
      << "The method RandomService::mySeed() is obsolete and should not be used.";
    return 1;
  }; 

  void RandomNumberService::snapShot(){
    if ( _states.size() == 0 ) return;

    // See note 2.
    // For now just do the HepRandom singleton.
    vector<uint32_t>& compressed = _states[0];
    vectorTransform< unsigned long, uint32_t>( CLHEP::HepRandom::getTheEngine()->put(), compressed);

  }
  
  void RandomNumberService::restoreState(const edm::Event& event){

    if ( _restoreStateLabel.size() == 0 ) return;

    edm::Handle<vector<RandomEngineState> > rns;
    event.getByLabel(_restoreStateLabel,rns);

    // For now just do the HepRandom singleton, which is stored in slot 0.
    const RandomEngineState&  state0 = (*rns)[0];
      
    // Restore the state. See notes 2 and 3.
    vector<unsigned long> v;
    vectorTransform( state0.getState(), v);
    bool status = CLHEP::HepRandom::getTheEngine()->getState(v);
    if ( !status ){
      throw cms::Exception("RANGE")
        << "Failed during restore of state of engine for: "
        << state0.getLabel() << " "
        << state0.getState().size();
    }

    // Private data must match the restored state of the engines.
    _states[0] = state0.getState();

    // Printout if requested.
    static int ncalls(-1);
    if ( _debug && ++ncalls < _nPrint ) {
      print();
    }

  }

  
  void RandomNumberService::print(){
    
    edm::LogInfo log("RANDOM");
    if ( _restoreStateLabel.size() != 0 ){
      log << "Name of module that created stored state: "
          << _restoreStateLabel
          << "\n";
    }else{
      log << "Will not restore state from the event. \n";
    }
    for ( size_t i=0; i<_labels.size(); ++i ){
      log << "Module label for engine: " 
          << _labels[i]
          << "  Size of state: " 
          << _states[i].size()
          << "  Size of seeds: " 
          << _seeds[i].size()
          << "\n";
    }
  }

  void RandomNumberService::saveEngineState(const std::string& fileName){
  }

  void RandomNumberService::restoreEngineState(const std::string& fileName){
  }
  
} // end namespace mu2e
