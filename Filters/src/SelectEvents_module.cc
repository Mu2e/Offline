//
// Filter events whit killed tracks.
//

// Mu2e includes.

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Root includes
//#include "TNtuple.h"
#include "TTree.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <vector>

using namespace std;

namespace mu2e {

  class SelectEvents : public art::EDFilter {
  public:
    explicit SelectEvents(fhicl::ParameterSet const& pset);
    virtual ~SelectEvents() { }

    bool filter( art::Event& event);

  private:

    std::vector<unsigned> _events;
    bool _verbose;
  };

  SelectEvents::SelectEvents(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _events(pset.get<std::vector<unsigned> >("events")),
    _verbose(pset.get<bool>("verbose",false)) {
  }

  bool SelectEvents::filter(art::Event& event) {

    unsigned ievent = event.event();
    if(std::find(_events.begin(),_events.end(),ievent) != _events.end()){
      if(_verbose) std::cout << "selected event " << ievent << std::endl;
      return true;
    } else
      return false;

  }

}

using mu2e::SelectEvents;
DEFINE_ART_MODULE(SelectEvents);
