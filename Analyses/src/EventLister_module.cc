//
// Write the event ids of all events.
//
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
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

using namespace std;

namespace mu2e {

  class EventLister : public art::EDAnalyzer {
  public:

    explicit EventLister(fhicl::ParameterSet const& pset);

    virtual void analyze ( const art::Event& event);

  private:
    int _ievent;
  };

  EventLister::EventLister(fhicl::ParameterSet const& pset):
        art::EDAnalyzer(pset)
  {
    _ievent = 0;
  }

  void EventLister::analyze(const art::Event& event) {
    _ievent++;
    printf("Event: %8i run: %10i subRun: %5i event: %10i\n",
	   _ievent,event.run(),event.subRun(),event.event());
    //    cout << "Event: " << event.id() << endl;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::EventLister);
