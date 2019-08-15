// Ntuple dumper for SimParticles with compact output
//
// Zhengyun You, 2013-12-01

#include <string>
#include <vector>
#include <set>
#include <limits>
#include <cmath>

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

using namespace std;

namespace mu2e {

  class Mu2eCoordinates {};

  //================================================================
  struct VDHit {
    float x;
    float y;
    float z;
    float time;
    float px;
    float py;
    float pz;
    int   pdgId;

    VDHit() : x(std::numeric_limits<double>::quiet_NaN())
            , y(std::numeric_limits<double>::quiet_NaN())
            , z(std::numeric_limits<double>::quiet_NaN())
            , time(std::numeric_limits<double>::quiet_NaN())
            , px(std::numeric_limits<double>::quiet_NaN())
            , py(std::numeric_limits<double>::quiet_NaN())
            , pz(std::numeric_limits<double>::quiet_NaN())
            , pdgId(0)
    {}

    //----------------------------------------------------------------
    // Constructor using mu2e coords
    VDHit(const Mu2eCoordinates& , const SimParticle& particle)
      : x(particle.startPosition().x())
      , y(particle.startPosition().y())
      , z(particle.startPosition().z())
      , time(particle.startGlobalTime())
      , px(particle.startMomentum().x())
      , py(particle.startMomentum().y())
      , pz(particle.startMomentum().z())
      , pdgId(particle.pdgId())
    {}

  }; // struct VDHit

  //================================================================
  class SimParticleDumperCompact : public art::EDAnalyzer {
    art::InputTag particlesInputTag_;

    typedef vector<int> Vint;
    // List of particles of interest for the particles ntuple
    set<int> pdg_save;
    string inProcessCodeDrop_;
    ProcessCode inProcessCodeToLook_;

    // Members needed to write the ntuple
    TTree *nt_;
    VDHit particle_;

  public:
    explicit SimParticleDumperCompact(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  SimParticleDumperCompact::SimParticleDumperCompact(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , particlesInputTag_(pset.get<std::string>("particlesInputTag"))
    , inProcessCodeDrop_(pset.get<string>("inputProcessCodeDrop"))
    , nt_(0)
  {
    Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
    if( pdg_ids.size()>0 ) {
      cout << "SimParticleDumperCompact: save following particle types in the ntuple: ";
      for( size_t i=0; i<pdg_ids.size(); ++i ) {
        pdg_save.insert(pdg_ids[i]);
        cout << pdg_ids[i] << ", ";
      }
      cout << endl;
    }

    inProcessCodeToLook_ = ProcessCode::findByName(inProcessCodeDrop_);
  }

  //================================================================
  void SimParticleDumperCompact::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    static const char branchDesc[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pdgId/I";
    nt_ = tfs->make<TTree>("nt", "particle ntuple");
    nt_->Branch("particles", &particle_, branchDesc);
  }

  //================================================================
  void SimParticleDumperCompact::analyze(const art::Event& event) {
    const auto& ih = event.getValidHandle<SimParticleCollection>(particlesInputTag_);
    for(const auto& i : *ih) {
      const SimParticle &aParticle = i.second;
      particle_ = VDHit(Mu2eCoordinates(), aParticle);

      //Check the Pdg of the SimParticle
      if (pdg_save.find(aParticle.pdgId()) == pdg_save.end()) continue;

      //Check the SimParticle's creationCode
      if (inProcessCodeToLook_ == aParticle.creationCode()) continue;

      nt_->Fill();
    }

  } // analyze(event)

  //================================================================


} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleDumperCompact);
