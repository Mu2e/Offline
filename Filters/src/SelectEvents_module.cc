//
// Filter events whit killed tracks.
// $Id: SelectEvents_module.cc,v 1.1 2012/12/13 19:43:12 brownd Exp $
// $Author: brownd $
// $Date: 2012/12/13 19:43:12 $
//
// Contact person G. Tassielli.
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

  };

  SelectEvents::SelectEvents(fhicl::ParameterSet const& pset):
    _events(pset.get<std::vector<unsigned> >("events")) {
  }

  bool SelectEvents::filter(art::Event& event) {

    unsigned ievent = event.event();
    if(std::find(_events.begin(),_events.end(),ievent) != _events.end())
      return true;
    else
      return false;

  } // end of ::analyze.

}

using mu2e::SelectEvents;
DEFINE_ART_MODULE(SelectEvents);
