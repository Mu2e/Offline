//
// This module scans all subruns and events and at the end prints
// metadata suitable for inserting into a SAM metadata record.
// It finds the min and max run/subrun for subruns and 
// min and max run/event for events. It also counts events and keeps 
// a list of runs and subruns.
//
// Ray Culbertson
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"

#include <iostream>
#include <vector>
#include <algorithm>

namespace mu2e {

  class RunEventSubRun : public art::EDAnalyzer {

  public:

    explicit RunEventSubRun(fhicl::ParameterSet const& );
    void beginSubRun( art::SubRun const& subRun ) override;
    void analyze    ( art::Event const&  event  ) override;
    void endJob () override;

  private:

    // keep a list of subruns
    typedef std::vector<art::SubRunID> subvec;
    subvec _subruns;

    art::RunNumber_t _max_run_s; // min and max run number for subruns
    art::RunNumber_t _min_run_s;
    art::SubRunNumber_t _max_sub;
    art::SubRunNumber_t _min_sub;
    art::RunNumber_t _max_run_e; // min and max run number for events
    art::RunNumber_t _min_run_e;
    art::EventNumber_t _max_evt;
    art::EventNumber_t _min_evt;
    int _count;   // event count

  };

}

mu2e::RunEventSubRun::RunEventSubRun(fhicl::ParameterSet const& pset ):
  art::EDAnalyzer(pset),
  _max_run_s(0),_min_run_s(-1), // -1 is maxint for a unsigned int
  _max_sub(0),_min_sub(-1),
  _max_run_e(0),_min_run_e(-1),
  _max_evt(0),_min_evt(-1),
  _count(0){
}

void
mu2e::RunEventSubRun::beginSubRun( art::SubRun const& subRun ){

  // add this subrun to the list, if it is not already there
  art::SubRunID id = subRun.id();
  subvec::iterator beg = _subruns.begin();
  subvec::iterator end = _subruns.end();
  if( std::find(beg,end,id) == end ) _subruns.emplace_back(id);

  //min subRun
  if(subRun.run() < _min_run_s) {
    _min_run_s = subRun.run();
    _min_sub = subRun.subRun();
  } else if(subRun.run() == _min_run_s) {
    if(subRun.subRun() < _min_sub) {
      _min_sub = subRun.subRun();
    }
  }

  //max subRun
  if(subRun.run() > _max_run_s) {
    _max_run_s = subRun.run();
    _max_sub = subRun.subRun();
  } else if(subRun.run() == _max_run_s) {
    if(subRun.subRun() > _max_sub) {
      _max_sub = subRun.subRun();
    }
  }

}

void
mu2e::RunEventSubRun::analyze(art::Event const& event){

  _count++;

  // min event
  if(event.run() < _min_run_e) {
    _min_run_e = event.run();
    _min_evt = event.event();
  } else if(event.run() == _min_run_e) {
    if(event.event() < _min_evt) {
      _min_evt = event.event();
    }
  }

  // max event
  if(event.run() > _max_run_e) {
    _max_run_e = event.run();
    _max_evt = event.event();
  } else if(event.run() == _max_run_e) {
    if(event.event() > _max_evt) {
      _max_evt = event.event();
    }
  }
  
}

void mu2e::RunEventSubRun::endJob () {

  //order the subruns
  std::sort(_subruns.begin(),_subruns.end());

  std::cout << "start RunEventSubRun::endJob summary" << std::endl;
  std::cout << "{" << std::endl;
  std::cout << "  \"event_count\"      : " << _count     << "," << std::endl;
  std::cout << "  \"dh.first_run_subrun\" : " << _min_run_s << "," << std::endl;
  std::cout << "  \"dh.first_subrun\"     : " << _min_sub   << "," << std::endl;
  std::cout << "  \"dh.first_run_event\"  : " << _min_run_e << "," << std::endl;
  std::cout << "  \"dh.first_event\"      : " << _min_evt   << "," << std::endl;
  std::cout << "  \"dh.last_run_subrun\"  : " << _max_run_s << "," << std::endl;
  std::cout << "  \"dh.last_subrun\"      : " << _max_sub   << "," << std::endl;
  std::cout << "  \"dh.last_run_event\"   : " << _max_run_e << "," << std::endl;
  std::cout << "  \"dh.last_event\"       : " << _max_evt   << "," << std::endl;
  std::cout << "  \"runs\"             : [\n";
  for (auto sr : _subruns ){
    std::cout << "\t[";
    std::cout << sr.run() << "," << sr.subRun() << ",\"unknown\"]";
    if(sr!=_subruns.back()) std::cout << ",\n" ;
  }
  std::cout << "\n\t\t]" << std::endl;
  std::cout << "}" << std::endl;
  std::cout << "end RunEventSubRun::endJob summary" << std::endl;

}


DEFINE_ART_MODULE(mu2e::RunEventSubRun)
