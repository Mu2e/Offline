//
//  Select events based on filters run in previous processing
//
// Original author: David Brown (LBNL) 20 Jun 2019
//
// C++
#include <string>
#include <vector>
#include <iostream>
#include <format>

// art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e includes
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"

using namespace std;

namespace mu2e {

  //================================================================
  class TriggerResultsFilter : public art::EDFilter {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int>  diagLevel{ Name("DiagLevel"), Comment("Diagonstic Level"), 0};
        fhicl::Atom<bool> printFirst{ Name("PrintFirst"), Comment("Print the TriggerResults on the first event"), false};
        fhicl::Atom<bool> noFilter{ Name("NoFilter"), Comment("If true, do not filter any events"), false};
        fhicl::Atom<string> processName{Name("ProcessName"), Comment("Process which generated TriggerResults")};
        fhicl::Sequence<string> triggerNames{ Name("TriggerNames"),
          Comment("Trigger line names to test; if any of these are set the event will pass the filter"), vector<string>()};
        fhicl::Sequence<unsigned> triggerBits{ Name("TriggerBits"),
          Comment("Trigger line bits to test; if any of these are set the event will pass the filter"), vector<unsigned>()};
      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit TriggerResultsFilter(Parameters const& config);
      virtual bool filter(art::Event& event) override;
      virtual void endJob() override;

    private:
      int _diag;
      string _pname; // process name for the TriggerResults object to test
      bool _pfirst; // print lines on first event
      bool _noFilter; // if true, do not filter any events
      std::vector<string> _tnames; // trigger line names; if any of these lines are set, accept the event
      std::vector<unsigned> _tbits; // trigger line bits; if any of these lines are set, accept the event
      std::vector<unsigned> _nset; // number of events passing each line
      unsigned _nevts, _npassed;
  };

  //================================================================
  TriggerResultsFilter::TriggerResultsFilter(Parameters const& config)
    : art::EDFilter{config},
    _diag(config().diagLevel()),
    _pname(config().processName()),
    _pfirst(config().printFirst()),
    _noFilter(config().noFilter()),
    _tnames(config().triggerNames()),
    _tbits(config().triggerBits()),
    _nset(max(_tnames.size(), _tbits.size()),0),
    _nevts(0),_npassed(0)
    {
      // Only filter on trigger bits or trigger names, not both
      if(_tnames.size()>0 && _tbits.size()>0){
        throw cet::exception("CONFIG")<<"mu2e::TriggerResultsFilter: cannot specify both trigger names and trigger bits to filter on";
      }
      if(_tnames.empty() && _tbits.empty() && !_noFilter){
        throw cet::exception("CONFIG")<<"mu2e::TriggerResultsFilter: must specify at least one trigger name or bit to filter on, or set NoFilter to true";
      }
    }

  //================================================================
  bool TriggerResultsFilter::filter(art::Event& event) {
    _nevts++;
    // find the TriggerResults object for the requested process name
    art::InputTag const tag{format("TriggerResults::{}", _pname.c_str())};
    auto trigResultsH = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults* trigResults = trigResultsH.product();
    TriggerResultsNavigator tnav(trigResults);
    if(_pfirst || _diag > 2){
      _pfirst = false;
      tnav.print();
    }
    // loop over all the lines in this TriggerResults and see if any of the requested are set.
    // Count each line separately for diagnostics
    bool passed(false);
    const bool use_bits = _tnames.empty();
    const size_t nbits = (use_bits) ? _tbits.size() : _tnames.size();
    for(size_t itrig = 0; itrig < nbits; itrig++) {
      try {
        auto const& tname = (use_bits) ? tnav.getTrigNameByBit(_tbits[itrig]) : _tnames[itrig];
        const bool accepted = tnav.accepted(tname);
        if(_diag>1) printf("[TriggerResultsFilter::%s] Trigger %s (bit %zu) accepted = %o\n",
                           __func__, tname.c_str(), tnav.getTrigBit(tname), accepted);
        if(accepted) {
          passed = true;
          _nset[itrig]++;
        }
      } catch (...) { // don't require the trigger to exist, as it may only appear in some runs
        if(_diag>0) printf("[TriggerResultsFilter::%s] Trigger list index %zu not found\n",
                            __func__, itrig);
      }
    }
    if(passed) _npassed++;
    return passed || _noFilter;
  }

  void TriggerResultsFilter::endJob() {
    printf("[TriggerResultFilter::%s] Processed %u events, accepted %u events (%8.4f%%)\n",
           __func__, _nevts, _npassed, _nevts > 0 ? _npassed*100./_nevts : 0.);
    if(_diag > 0) {
      if(_tnames.size()>0) {
        for(size_t iname=0;iname<_tnames.size(); iname++){
          cout << "  Trigger Line " << _tnames[iname] << " set for " << _nset[iname] << " events" << endl;
        }
      } else if(_tbits.size()>0) {
        for(size_t ibit=0;ibit<_tbits.size(); ibit++){
          cout << "  Trigger Bit " << _tbits[ibit] << " set for " << _nset[ibit] << " events" << endl;
        }
      }
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TriggerResultsFilter)
