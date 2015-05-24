//
//  Module to summarize multiple ProtonBunchIntensity objects in the event.
//  Event Mixing can result in multiple such objects in the event, which can
//  be compared to test for consistency.  This module also creates an event
//  weight that can be used by downstream analysis modules.
//  Original author: David Brown (LBNL) 22 May 2015
//
// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
// general Mu2e includes
// data products produced by this module
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include <iostream>
#include <vector>

namespace mu2e {
  class ProtonBunchIntensitySummarizer : public art::EDProducer {
    public:
      explicit ProtonBunchIntensitySummarizer(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
    private:
      int _printLevel; // level of diagnostic printout
  };


  ProtonBunchIntensitySummarizer::ProtonBunchIntensitySummarizer(const fhicl::ParameterSet& pset) :
    _printLevel(pset.get<int>("PrintLevel",0))
  {
    produces<mu2e::EventWeight>();
  }

  void ProtonBunchIntensitySummarizer::produce(art::Event& event) {
     // Get all of the tracker ProtonBunchIntensity objects from the event:
    typedef std::vector< art::Handle<ProtonBunchIntensity> > HandleVector;
    HandleVector pbiHandles;
    event.getManyByType(pbiHandles);
    if(pbiHandles.empty()){
      throw cet::exception("SIM")<<"mu2e::ProtonBunchIntensitySummarizer: No ProtonBunchIntensity objects found" << std::endl;
    }
    // Loop over ProtonBunchIntensity objects
    bool first(true);
    ProtonBunchIntensity pbi;
    for ( HandleVector::const_iterator ipbi=pbiHandles.begin(), epbi=pbiHandles.end();ipbi != epbi; ++ipbi ){
      art::Handle<ProtonBunchIntensity> const& pbihandle(*ipbi);
      if(_printLevel > 1)
	std::cout << "ProtonBunchIntensity object found, Nprotons = " << pbihandle->intensity() << ", Mean Intensity = " << pbihandle->meanIntensity() << std::endl;
// if this is the 1st object, set the cached value
      if(first){
	pbi = *pbihandle;
	first = false;
      } else {
// check for consistency in the event
	if(pbi != *pbihandle){
	  throw cet::exception("SIM")<<"mu2e::ProtonBunchIntensitySummarizer: Inconsistent ProtonBunchIntensity objects found" << std::endl;
	}
      }
    }
  // downstream modules need to weight any process that depends on the proton bunch intensity, even
    // if they are only generating 1/event (like conversion electrons), as the probability of producing
    // that one event scales with proton intensity
    std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(pbi.weight()) );
    if(_printLevel > 0){
      std::cout << "Found " << pbi.intensity() << " protons in this microbunch, for event weight = "
      << evtwt->weight() << std::endl;
    }
    event.put(std::move(evtwt));
  }
}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensitySummarizer);

