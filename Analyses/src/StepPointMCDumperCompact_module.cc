// Ntuple dumper for StepPointMC with compact output
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

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

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
    VDHit(const Mu2eCoordinates& , const StepPointMC& hit)
      : x(hit.position().x())
      , y(hit.position().y())
      , z(hit.position().z())
      , time(hit.time())
      , px(hit.momentum().x())
      , py(hit.momentum().y())
      , pz(hit.momentum().z())
      , pdgId(hit.simParticle()->pdgId())
    {}

  }; // struct VDHit

  //================================================================
  class StepPointMCDumperCompact : public art::EDAnalyzer {
    typedef std::vector<art::InputTag> InputTags;
    InputTags stepInputs_; 

    typedef vector<int> Vint;
    // List of particles of interest for the particles ntuple
    set<int> pdg_save;
    string inProcessCodeDrop_;
    ProcessCode inProcessCodeToLook_;

    // Members needed to write the ntuple
    TTree *nt_;
    VDHit particle_;

  public:
    explicit StepPointMCDumperCompact(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  StepPointMCDumperCompact::StepPointMCDumperCompact(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inProcessCodeDrop_(pset.get<string>("inputProcessCodeDrop"))
    , nt_(0)
  {
    typedef std::vector<std::string> VS;
    const VS stepStrings(pset.get<VS>("stepInputs"));
    for(const auto& i : stepStrings) {
      stepInputs_.emplace_back(i);
    }

    Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
    if( pdg_ids.size()>0 ) {
      cout << "StepPointMCDumperCompact: save following particle types in the ntuple: ";
      for( size_t i=0; i<pdg_ids.size(); ++i ) {
        pdg_save.insert(pdg_ids[i]);
        cout << pdg_ids[i] << ", ";
      }
      cout << endl;
    }

    inProcessCodeToLook_ = ProcessCode::findByName(inProcessCodeDrop_);
  }

  //================================================================
  void StepPointMCDumperCompact::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    static const char branchDesc[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pdgId/I";
    nt_ = tfs->make<TTree>("nt", "particle ntuple");
    nt_->Branch("particles", &particle_, branchDesc);
  }

  //================================================================
  void StepPointMCDumperCompact::analyze(const art::Event& event) {

    for(const auto& tag : stepInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      for(const auto& i : *ih) {
        const art::Ptr<SimParticle>& particle(i.simParticle());
        //Check the Pdg of the StepPointMC
        if (pdg_save.find(particle->pdgId()) == pdg_save.end()) continue;

        //Check the StepPointMC's creationCode
        if (inProcessCodeToLook_ == particle->creationCode()) continue;

        particle_ = VDHit(Mu2eCoordinates(), i);
        nt_->Fill();
      }
    }

  } // analyze(event)

  //================================================================


} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StepPointMCDumperCompact);
