//
// Store state of RandomNumberService into the event. 
//
// $Id: RandomNumberSaver_plugin.cc,v 1.1 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/05 16:07:38 $
//
// Original author Rob Kutschke
//
// Notes
// 1) In the present implementation both StateInfo and SeedInfo are typedefs
//    to the same underlying type.  The edm sees the underlying type.
//    When we put both StateInfo and SeedInfo into the event, the event
//    sees two data products of identical type; so we are required to 
//    distinguish them using the string argument to produces().
//

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

  };

  RandomNumberSaver::RandomNumberSaver(edm::ParameterSet const& pset){

    throw cms::Exception("NotReady")
      << "The RandomNumberSaver module does not yet work. \n"
      << "Please remove it from your configuration.";

    // This produces three data products.  See note 1.
    /*
    produces<RandomNumberService::LabelInfo>();
    produces<RandomNumberService::StateInfo>("StateInfo");
    produces<RandomNumberService::SeedInfo>("SeedInfo");
    */

  }


  void
  RandomNumberSaver::produce( edm::Event& event, edm::EventSetup const&) {

    // Get the data from the RandomNumberService and put it into a data product.
    // The edm requires that we make copies and put the copies into the event.
    edm::Service<edm::RandomNumberGenerator> rng;
    auto_ptr<RandomNumberService::LabelInfo> labels(new RandomNumberService::LabelInfo(rng->getCachedLabels()));
    auto_ptr<RandomNumberService::StateInfo> states(new RandomNumberService::StateInfo(rng->getCachedStates()));
    auto_ptr<RandomNumberService::SeedInfo>  seeds (new RandomNumberService::SeedInfo (rng->getCachedSeeds()));

    event.put(labels);
    event.put(states);
    event.put(seeds);

    cerr << "Putting stuff into the event." << endl;

  }
  
}


using mu2e::RandomNumberSaver;
DEFINE_FWK_MODULE(RandomNumberSaver);
