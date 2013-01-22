//
// Write the event ids of all events.
//
// $Id: EventLister_module.cc,v 1.1 2013/01/22 00:05:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/22 00:05:10 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

using namespace std;

namespace mu2e {

  class EventLister : public art::EDAnalyzer {
  public:

    explicit EventLister(fhicl::ParameterSet const& pset);

    virtual void analyze ( const art::Event& event);

  private:

  };

  EventLister::EventLister(fhicl::ParameterSet const& pset) {}

  void EventLister::analyze(const art::Event& event) {
    cout << "Event: " << event.id() << endl;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::EventLister);
