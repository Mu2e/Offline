//
// Plugin to read virtual detectors data and create ntuples
//
//
// Original author Ivan Logashenko
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/Mu2eUtilities/inc/fromStrings.hh"
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;
namespace mu2e {
  class MakeTree : public art::EDAnalyzer {
    typedef SimParticleCollection::key_type key_type;
    art::InputTag _vdInputTag = "compressDetStepMCsSTM:virtualdetector";
    art::InputTag _simpInputTag = "compressDetStepMCsSTM:";

    TTree* _ntvd;
    Float_t time = 0.0;
    Int_t pdgId = 0;
    Float_t x = 0.0;
    Float_t y = 0.0;
    Float_t E = 0.0;
    key_type trackId;
  public:

    explicit MakeTree(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);

  };

  MakeTree::MakeTree(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset)
  {
    art::ServiceHandle<art::TFileService> tfs;
    _ntvd = tfs->make<TTree>( "ntvd", "Virtual Detectors ntuple");
    _ntvd->Branch("time", &time, "time/F");
    _ntvd->Branch("pdgId", &pdgId, "pdgId/I");
    _ntvd->Branch("x", &x, "x/F");
    _ntvd->Branch("y", &y, "y/F");
    _ntvd->Branch("E", &E, "E/F");

  }

  void MakeTree::analyze(const art::Event& event) {
    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_vdInputTag,hits);
    bool haveStepPart = hits.isValid();

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_simpInputTag, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());
    if (!(haveSimPart && haveStepPart))
      std::cout << "Event wierd"<< std::endl;

    // Loop over all VD hits.
    if( hits.isValid() ) for ( size_t i=0; i<hits->size(); ++i ){
        const StepPointMC& hit = (*hits)[i];
        if(hit.volumeId() != 101)
            continue;
        // Get the associated particle
        trackId = hit.trackId();
        SimParticle const& sim = simParticles->at(trackId);
        time = hit.time();
        pdgId = sim.pdgId();
        x = hit.position().x();
        y = hit.position().y();
        E = sim.endKineticEnergy();
        _ntvd->Fill();
      }
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MakeTree)
