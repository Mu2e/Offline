//
//  Select events based on filters run in previous processing
// 
// Original author: David Brown (LBNL) 20 Jun 2019
//
// C++
#include <string>
#include <vector>
#include <iostream>
// root
#include "TString.h"
// art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "cetlib_except/exception.h"
// mu2e includes
#include "Mu2eUtilities/inc/TriggerResultsNavigator.hh"

using namespace std;

namespace mu2e {

  //================================================================
  class TriggerResultsFilter : public art::EDFilter {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<int> diagLevel{ Name("DiagLevel"), Comment("Diagonstic Level"), 0};
	fhicl::Atom<bool> printFirst{ Name("PrintFirst"),
	  Comment("Print the TriggerResults on the first event"), false};
	fhicl::Atom<string> processName{Name("ProcessName"), Comment("Process which generated TriggerResults")};
	fhicl::Sequence<string> triggerNames{ Name("TriggerNames"),
	  Comment("Trigger line names to test; if any of these are set the event will pass the filter")};
      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit TriggerResultsFilter(Parameters const& config);
      virtual bool filter(art::Event& event) override;
      virtual void endJob() override;

    private:
      int _diag;
      string _pname; // process name for the TriggerResults object to test
      bool _pfirst; // print lines on first event
      std::vector<string> _tnames; // trigger line names; if any of these lines are set, accept the event
      std::vector<unsigned> _nset; // number of events passing each line
      unsigned _nevts, _npassed;
  };

  //================================================================
  TriggerResultsFilter::TriggerResultsFilter(Parameters const& config)
    : art::EDFilter{config},
    _diag(config().diagLevel()),
    _pname(config().processName()),
    _pfirst(config().printFirst()),
    _tnames(config().triggerNames()),
    _nset(_tnames.size(),0),
    _nevts(0),_npassed(0)
    {}

  //================================================================
  bool TriggerResultsFilter::filter(art::Event& event) {
    _nevts++;
    // find the TriggerResults object for the requested process name
    art::InputTag const tag{Form("TriggerResults::%s", _pname.c_str())};
    auto trigResultsH = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults* trigResults = trigResultsH.product();
    TriggerResultsNavigator tnav(trigResults);
    if(_pfirst){
      _pfirst = false;
      tnav.print(); 
    }
    // loop over all the lines in this TriggerResults and see if any of the requested are set.
    // Count each line separately for diagnostics
    bool passed(false);
    for(size_t iname=0;iname < _tnames.size(); iname++){
      auto const& tname = _tnames[iname];
      size_t itrig = tnav.findTrigPath(tname);
      // require that the line exist
      if(itrig == trigResults->size())throw cet::exception("Filter")<<"mu2e::TriggerResultsFilter: cannot find TriggerResults value for trigger " <<  tname << endl;
    if(_diag>0) cout << "trigger line " << itrig << " found for name " << tname << " with value " << tnav.getTrigPath(itrig) << " status "
      << trigResults->accept(itrig) << endl;

      if(trigResults->accept(itrig)){
	passed = true;
	_nset[iname]++;
      }
    }
    if(passed)_npassed++;
    return passed;
  }

  void TriggerResultsFilter::endJob() {
    cout << "filter passed " << _npassed << " of " << _nevts << " events" << endl;
    for(size_t iname=0;iname<_tnames.size(); iname++){
      cout << "Trigger Line " << _tnames[iname] << " set for " << _nset[iname] << " events" << endl;
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TriggerResultsFilter)
