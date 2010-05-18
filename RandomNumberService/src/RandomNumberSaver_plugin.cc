//
// Store state of RandomNumberService into the event. 
//
// $Id: RandomNumberSaver_plugin.cc,v 1.3 2010/05/18 21:16:42 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 21:16:42 $
//
// Original author Rob Kutschke
//
// Notes

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "RandomNumberService/inc/RandomNumberService.hh"
#include "ToyDP/inc/RandomEngineState.hh"

using namespace std;

namespace mu2e {

  class RandomNumberSaver : public edm::EDProducer {

  public:
    explicit RandomNumberSaver(edm::ParameterSet const& pset);
    virtual ~RandomNumberSaver() { }
    
    void produce( edm::Event& e, edm::EventSetup const&);
    
    static void fillDescription(edm::ParameterSetDescription& iDesc,
                                string const& moduleLabel) {
      iDesc.setAllowAnything();
    }
    
  private:

    bool _debug;
    
  };
  
  RandomNumberSaver::RandomNumberSaver(edm::ParameterSet const& pset):
    _debug(pset.getUntrackedParameter<bool>("debug",false)){
    produces<std::vector<RandomEngineState> >();
  }

  void RandomNumberSaver::produce( edm::Event& event, edm::EventSetup const&) {

    auto_ptr<vector<RandomEngineState> >  rnState(new vector<RandomEngineState>() );

    // References to access the data inside the RandomNumberService.
    edm::Service<edm::RandomNumberGenerator> rng;
    const RandomNumberService::LabelInfo& labels = rng->getCachedLabels();
    const RandomNumberService::StateInfo& states = rng->getCachedStates();
    const RandomNumberService::SeedInfo&  seeds  = rng->getCachedSeeds();

    // Check self consistency.
    if ( labels.size() != states.size() ||
         labels.size() != seeds.size()     ){
      throw cms::Exception("RANGE")
        << "Inconsistent sizes of objects in the state of the RandomNumberService: "
        << labels.size() << " "
        << states.size() << " "
        << seeds.size()  << "\n";
    }

    // For each saved engine, copy the state information to the data product.
    for ( size_t i=0; i<labels.size(); ++i){
      
      // Two phase construction to avoid the double copy.  Use emplace_back when available.
      rnState->push_back(RandomEngineState());
      RandomEngineState& s = rnState->back();
      s.setLabel( labels[i] );
      s.setState( states[i] );
      s.setSeed (  seeds[i] );

    }

    event.put(rnState);

    if ( _debug ){
      rng->print();
    }
  }

} // end namespace mu2e

using mu2e::RandomNumberSaver;
DEFINE_FWK_MODULE(RandomNumberSaver);
