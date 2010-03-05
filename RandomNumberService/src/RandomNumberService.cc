//
// Maintain multiple independent random number engines,
// including save and restore of seeds and state. 
//
// $Id: RandomNumberService.cc,v 1.2 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/05 16:07:38 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The CMS code is available at
//    http://cmslxr.fnal.gov/lxr/source/IOMC/RandomEngine/src/RandomNumberGeneratorService.cc

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

// CLHEP include files.
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/JamesRandom.h"

using namespace std;

namespace mu2e {

  RandomNumberService::RandomNumberService(edm::ParameterSet const& iPS, 
					   edm::ActivityRegistry&iRegistry) :
    _globalSeed(iPS.getUntrackedParameter<int>("globalSeed",-1)),
    _labels(),
    _states(),
    _seeds()
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

    // The default random number generator is HepJamesRandom.
    // A comment in the sources says that "seed" should be within 
    // the range [0,900000000]
    if ( _globalSeed > -1 && _globalSeed < 90000000 ){
      CLHEP::HepRandom::setTheSeed(_globalSeed);
    } else if ( _globalSeed != -1 ) {

      // Or do we want to just give a warning message?
      throw cms::Exception("RANGE")
	<< "Seed for the HepRandom generator is out of bounds: "
	<< _globalSeed;
    }
    
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
    // cerr << "xxxxCalling getEngine: " << endl;

    // A hack to allow this to compile.
    static CLHEP::HepJamesRandom hack;
    return hack;
  }

  uint32_t RandomNumberService::mySeed() const {
    throw cms::Exception("Obsolete")
      << "The method RandomService::mySeed() is obsolete and should not be used.";
    return 1;
  };

  void RandomNumberService::snapShot(){
    // cerr << "xxxxCalling snapshot: " << endl;
  }
  
  void RandomNumberService::restoreState(const edm::Event& event){
    // cerr << "xxxxCalling restoreState: " << endl;
  }

  void RandomNumberService::print(){
    edm::LogInfo log("Random");
    log << "Information about the RandomNumberService"
	<< endl;
    log << "Number of saved engines: " << _labels.size() << endl;
  }

  void RandomNumberService::saveEngineState(const std::string& fileName){
  }

  void RandomNumberService::restoreEngineState(const std::string& fileName){
  }
  
} // end namespace mu2e
