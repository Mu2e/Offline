//
// Read the pid
//
// Original author Vadim Rusu
//

// Mu2e includes.

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// ROOT incldues
#include "TH1F.h"

// C++ includes.
#include <iostream>
#include <string>

#include "Offline/RecoDataProducts/inc/PIDProduct.hh"

//ROOT
#include "TTree.h"

using namespace std;

namespace mu2e {

  class ParticleIDRead : public art::EDAnalyzer {
  public:

    explicit ParticleIDRead(fhicl::ParameterSet const& pset);
    virtual ~ParticleIDRead() { }

    void beginJob();
    void analyze(const art::Event& e);

  private:

    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;

    // Control level of printout.
    int _verbosity;
    int _maxPrint;
    // whether or not to include MC info for empty events
    bool _processEmpty;

    // Histograms

    TTree* _diag;

    Int_t _trkid,_eventid;

  };

  ParticleIDRead::ParticleIDRead(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _verbosity(pset.get<int>("verbosity",0)),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _processEmpty(pset.get<bool>("processEmpty",true)),
    _diag(0){
// construct the data product instance name
  }

  void ParticleIDRead::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
    _diag = tfs->make<TTree>("PID", "PID info");
    _diag->Branch("eventid",&_eventid,"eventid/I");
    _diag->Branch("trkid",&_trkid,"trkid/I");

    _eventid = 0;
  }

  // For each event, look at tracker hits and calorimeter hits.
  void ParticleIDRead::analyze(const art::Event& event) {
//    cout << "Enter ParticleIDRead:: analyze: " << _verbosity << endl;

    _eventid++;

    art::Handle<PIDProductCollection> pidsHandle;
    event.getByLabel(_fitterModuleLabel,pidsHandle);
    PIDProductCollection const& pids = *pidsHandle;

    if ( _verbosity > 0 && _eventid <= _maxPrint ){
      cout << "ParticleIDRead  for event: " << event.id() << "  Number of pid objects: " << pids.size() << endl;
    }

    _trkid = -1;

    for ( size_t i=0; i< pids.size(); ++i ){
      _trkid = i;
      PIDProduct const& pid   = pids.at(i);

      if ( _verbosity > 1 && _eventid <= _maxPrint ){
        cout << "   PIDProduct: "
             << i            << " "
             << pid.GetTrkID() << " "
             << pid.GetResidualsSlope() << " "
             << pid.GetResidualsSlopeError() << " "
             << pid.GetLogEProb() << " "
             << pid.GetLogMProb() << " "
             << endl;
      }

    }

  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ParticleIDRead;
DEFINE_ART_MODULE(ParticleIDRead)
