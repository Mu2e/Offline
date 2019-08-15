//
// Filter events whit killed tracks.
// $Id: KilledEventFilter_module.cc,v 1.1 2012/04/20 07:04:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/04/20 07:04:28 $
//
// Contact person G. Tassielli.
//

// Mu2e includes.
#include "MCDataProducts/inc/StatusG4.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

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

  struct buffEventID {
          std::uint32_t ev;
          std::uint32_t subRun;
          std::uint32_t run;

          buffEventID():ev(0),subRun(0),run(0){}
          buffEventID(std::uint32_t tEv, std::uint32_t tSubRun, std::uint32_t tRun)
             :ev(tEv),subRun(tSubRun),run(tRun){}
  };

  class KilledEventFilter : public art::EDFilter {
  public:
    explicit KilledEventFilter(fhicl::ParameterSet const& pset);
    virtual ~KilledEventFilter() { }

    bool filter( art::Event& event);

    virtual bool beginRun(art::Run &run);
    virtual void endJob();

  private:

    // Module label of the module that made the StepPointMCCollections.
    std::string g4ModuleLabel_;
    //TNtuple* ntup_;
    TTree* ntup_;

    // Number of events that don't pass the filter.
    std::uint64_t _nBadG4Status, _nOverflow, _nKilled;
    std::uint32_t *_evListBadG4Status;
    std::uint32_t *_subRunListBadG4Status;
    std::uint32_t *_runListBadG4Status;
    std::uint32_t *_evListOverflow;
    std::uint32_t *_subRunListOverflow;
    std::uint32_t *_runListOverflow;
    std::uint32_t *_evListKilled;
    std::uint32_t *_subRunListKilled;
    std::uint32_t *_runListKilled;
    std::vector<mu2e::buffEventID> _buffEvListBadG4Status;
    std::vector<mu2e::buffEventID> _buffEvListOverflow;
    std::vector<mu2e::buffEventID> _buffEvListKilled;

    // Number of events that pass the filter.
    std::uint64_t _nPassed;

  };

  KilledEventFilter::KilledEventFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel")),
    ntup_(0),
    _nBadG4Status(0),
    _nOverflow(0),
    _nKilled(0),
    _evListBadG4Status(0x0),
    _subRunListBadG4Status(0x0),
    _runListBadG4Status(0x0),
    _evListOverflow(0x0),
    _subRunListOverflow(0x0),
    _runListOverflow(0x0),
    _evListKilled(0x0),
    _subRunListKilled(0x0),
    _runListKilled(0x0),
    _nPassed(0){
          _buffEvListBadG4Status.clear();
          _buffEvListOverflow.clear();
          _buffEvListKilled.clear();
  }

  bool KilledEventFilter::beginRun(art::Run& ){

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "EventSummary" );

    //ntup_ = tfdir.make<TNtuple>( "evntNtup", "Event Summary", "nPassed:nBadG4Status:nOverflow:nKilled");
    ntup_ = tfdir.make<TTree>( "evntNtup", "Event Summary");
    ntup_->Branch("nPassed",&_nPassed,"_nPassed/l");
    ntup_->Branch("nBadG4Status",&_nBadG4Status,"_nBadG4Status/l");
    ntup_->Branch("nOverflow",&_nOverflow,"_nOverflow/l");
    ntup_->Branch("nKilled",&_nKilled,"_nKilled/l");
//    ntup_->Branch("evListBadG4Status",_evListBadG4Status,"_evListBadG4Status[_nBadG4Status]/i");
//    ntup_->Branch("subRunListBadG4Status",_subRunListBadG4Status,"_subRunListBadG4Status[_nBadG4Status]/i");
//    ntup_->Branch("runListBadG4Status",_runListBadG4Status,"_runListBadG4Status[_nBadG4Status]/i");
//    ntup_->Branch("evListOverflow",_evListOverflow,"_evListOverflow[_nOverflow]/i");
//    ntup_->Branch("subRunListOverflow",_subRunListOverflow,"_subRunListOverflow[_nOverflow]/i");
//    ntup_->Branch("runListOverflow",_runListOverflow,"_runListOverflow[_nOverflow]/i");
//    ntup_->Branch("evListKilled",_evListKilled,"_evListKilled[_nKilled]/i");
//    ntup_->Branch("subRunListKilled",_subRunListKilled,"_subRunListKilled[_nKilled]/i");
//    ntup_->Branch("runListKilled",_runListKilled,"_runListKilled[_nKilled]/i");

    return true;
  }

  bool KilledEventFilter::filter(art::Event& event) {

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Accept only events with good status from G4.
    if ( g4Status.status() > 1 ) {
            ++_nBadG4Status;
            _buffEvListBadG4Status.push_back(buffEventID(event.event(),event.getSubRun().subRun(),event.getRun().run()));
            mf::LogError("G4")
            <<"Aborting in Run "<<event.getRun().run()<<" SubRun "<<event.getSubRun().subRun()<<" event "<<event.event()
            <<" KilledEventFilter::filter due to G4 status\n"
            << g4Status;
            return false;
    }

    if (g4Status.overflowSimParticles()) {
            ++_nOverflow;
            _buffEvListOverflow.push_back(buffEventID(event.event(),event.getSubRun().subRun(),event.getRun().run()));
            mf::LogError("G4")
            <<"Aborting in Run "<<event.getRun().run()<<" SubRun "<<event.getSubRun().subRun()<<" event "<<event.event()
            <<" KilledEventFilter::filter due to overflow of particles\n"
            << g4Status;
            return false;
    }

    if (g4Status.nKilledStepLimit() > 0) {
            ++_nKilled;
            _buffEvListKilled.push_back(buffEventID(event.event(),event.getSubRun().subRun(),event.getRun().run()));
            mf::LogError("G4")
            <<"Aborting in Run "<<event.getRun().run()<<" SubRun "<<event.getSubRun().subRun()<<" event "<<event.event()
            <<" KilledEventFilter::filter due to nkilledStepLimit reached\n"
            << g4Status;
            return false;
    }

    _nPassed++;
    return true;

  } // end of ::analyze.

  void KilledEventFilter::endJob() {

    /*float nt[ntup_->GetNvar()];
    nt[0] = _nPassed;
    nt[1] = _nBadG4Status;
    nt[2] = _nOverflow;
    nt[3] = _nKilled;

    ntup_->Fill(nt);*/

    _evListBadG4Status      = new std::uint32_t[_nBadG4Status];
    _subRunListBadG4Status  = new std::uint32_t[_nBadG4Status];
    _runListBadG4Status     = new std::uint32_t[_nBadG4Status];
    _evListOverflow         = new std::uint32_t[_nOverflow];
    _subRunListOverflow     = new std::uint32_t[_nOverflow];
    _runListOverflow        = new std::uint32_t[_nOverflow];
    _evListKilled           = new std::uint32_t[_nKilled];
    _subRunListKilled       = new std::uint32_t[_nKilled];
    _runListKilled          = new std::uint32_t[_nKilled];
    for (std::uint64_t iBS=0; iBS<_nBadG4Status; iBS++) {
            _evListBadG4Status[iBS]     = _buffEvListBadG4Status.at(iBS).ev;
            _subRunListBadG4Status[iBS] = _buffEvListBadG4Status.at(iBS).subRun;
            _runListBadG4Status[iBS]    = _buffEvListBadG4Status.at(iBS).run;
    }
    for (std::uint64_t iOF=0; iOF<_nOverflow; iOF++) {
            _evListOverflow[iOF]     = _buffEvListOverflow.at(iOF).ev;
            _subRunListOverflow[iOF] = _buffEvListOverflow.at(iOF).subRun;
            _runListOverflow[iOF]    = _buffEvListOverflow.at(iOF).run;
    }
    for (std::uint64_t iKL=0; iKL<_nKilled; iKL++) {
            _evListKilled[iKL]     = _buffEvListKilled.at(iKL).ev;
            _subRunListKilled[iKL] = _buffEvListKilled.at(iKL).subRun;
            _runListKilled[iKL]    = _buffEvListKilled.at(iKL).run;
    }

    ntup_->Branch("evListBadG4Status",_evListBadG4Status,"_evListBadG4Status[_nBadG4Status]/i");
    ntup_->Branch("subRunListBadG4Status",_subRunListBadG4Status,"_subRunListBadG4Status[_nBadG4Status]/i");
    ntup_->Branch("runListBadG4Status",_runListBadG4Status,"_runListBadG4Status[_nBadG4Status]/i");
    ntup_->Branch("evListOverflow",_evListOverflow,"_evListOverflow[_nOverflow]/i");
    ntup_->Branch("subRunListOverflow",_subRunListOverflow,"_subRunListOverflow[_nOverflow]/i");
    ntup_->Branch("runListOverflow",_runListOverflow,"_runListOverflow[_nOverflow]/i");
    ntup_->Branch("evListKilled",_evListKilled,"_evListKilled[_nKilled]/i");
    ntup_->Branch("subRunListKilled",_subRunListKilled,"_subRunListKilled[_nKilled]/i");
    ntup_->Branch("runListKilled",_runListKilled,"_runListKilled[_nKilled]/i");

    ntup_->Fill();

    mf::LogInfo("Summary") 
      << "KilledEventFilter_module: Number of events passing the filter: "
      << _nPassed
      << "\nNumber of events skipped: "
      << "due to G4 completion status: "
      << _nBadG4Status
      << "\nKilledEventFilter::endJob Number of overflow events "
      << "due to too many particles in G4: "
      << _nOverflow
      << "\nKilledEventFilter::endJob Number of events with killed particles "
      << "due to too many steps in G4: "
      << _nKilled
      << endl;
  }

}

using mu2e::KilledEventFilter;
DEFINE_ART_MODULE(KilledEventFilter);
