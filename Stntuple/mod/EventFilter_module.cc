//
//  Filter events passing only th ones from teh list (A module to look at the provenance of a product.
//
//  $Id: EventFilter_module.cc,v 1.4 2013/10/21 21:15:46 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/10/21 21:15:46 $
//
//  Original author Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class EventFilter : public art::EDFilter {
  protected: 
    vector<int>       _eventList;  // 3 numbers per event: [run, subrun, event]
    int               _nevents;

  public:
    explicit EventFilter(fhicl::ParameterSet const& pset);

    virtual bool filter(art::Event& anEvent) override;

  private:

    void printProvenance( art::Provenance const& );

  };

  EventFilter::EventFilter(fhicl::ParameterSet const& pset)
    : art::EDFilter(),
      _eventList(pset.get<vector<int>>("eventList"))
  {
    _nevents = _eventList.size()/3;
  }


//-----------------------------------------------------------------------------
  bool EventFilter::filter(art::Event&  anEvent) {

    uint   eventNumber, runNumber, subRunNumber;
    bool   found(true);

    if (_nevents > 0) {
      found = false;

      for (int i=0; i<_nevents; i++) {
	runNumber    = _eventList[3*i+0];
	eventNumber  = _eventList[3*i+1];
	subRunNumber = _eventList[3*i+2];
	
	if ((anEvent.event () == eventNumber ) && 
	    (anEvent.run   () == runNumber   ) && 
	    (anEvent.subRun() == subRunNumber)    ) {
	  found = true;
	}
      }
    }

    return found;
  } 

} // end namespace mu2e

using mu2e::EventFilter;
DEFINE_ART_MODULE(EventFilter);
