//
// simple module to read stepPoints in PTM
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
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
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadPTM : public art::EDAnalyzer {
    public:

      typedef vector<int> Vint;
      typedef SimParticleCollection::key_type key_type;

      explicit ReadPTM(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _ptmStepPoints(pset.get<string>("ptmStepPoints","PTM")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _ntPTargetMon(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run"))
    {

      Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
      if( pdg_ids.size()>0 ) {
        cout << "ReadPTM: save following particle types in the ntuple: ";
        for( size_t i=0; i<pdg_ids.size(); ++i ) {
          pdg_save.insert(pdg_ids[i]);
          cout << pdg_ids[i] << ", ";
        }
        cout << endl;
      }

      nt = new float[1000];

    } // ReadPTM()

    virtual ~ReadPTM() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Name of the StepPoint collections
    std::string  _ptmStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntPTargetMon;

    float *nt; // Need this buffer to fill TTree ntPTargetMon

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;
  }; // class ReadPTM

  void ReadPTM::beginJob(){

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntPTargetMon = tfs->make<TNtuple>( "ntPTM", "PTM ntuple",
                                "run:evt:volId:trk:pdg:time:x:y:z:px:py:pz:iedep:totedep:"
                                "gtime");
  }

  void ReadPTM::beginRun(art::Run const& run){

  }

  void ReadPTM::analyze(const art::Event& event) {

    ++_nAnalyzed;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_ptmStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    // Loop over all hits.
    int validCount = 0;
    if( hits.isValid() ) for ( size_t i=0; i<hits->size(); ++i ){
      validCount++;
      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get track info
      key_type trackId = hit.trackId();
      int pdgId = 0;
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
        }
      }

      // Fill the ntuple.
      nt[0]  = event.id().run();
      nt[1]  = event.id().event();
      nt[2]  = hit.volumeId();
      nt[3]  = trackId.asInt();
      nt[4]  = pdgId;
      nt[5]  = hit.time();
      nt[6]  = pos.x();
      nt[7]  = pos.y();
      nt[8]  = pos.z();
      nt[9]  = mom.x();
      nt[10] = mom.y();
      nt[11] = mom.z();
      nt[12] = hit.ionizingEdep();
      nt[13] = hit.totalEDep();
      nt[14] = hit.properTime();

      _ntPTargetMon->Fill(nt);
      if ( _nAnalyzed < _maxPrint){
        cout << "PTM hit: "
             << event.id().run()   << " | "
             << event.id().event() << " | "
             << hit.volumeId()     << " "
             << pdgId              << " | "
             << hit.time()         << " "
             << mom.mag()
             << endl;

      }

    } // end loop over hits.

  } // ReadPTM::analyze

} // namespace mu2e

//using mu2e::ReadPTM;
DEFINE_ART_MODULE(mu2e::ReadPTM)
