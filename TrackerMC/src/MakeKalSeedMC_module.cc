//
//  Module to create KalSeed objects from MC truth  These can be used to test analyses or algorithms against truth
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
// utiliities
// persistent data
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include <vector>
using namespace std;
namespace mu2e {
  namespace TrackerMC {
    class MakeKalSeedMC : public art::EDProducer {
      struct SDMCSimCol {
        art::Ptr<SimParticle> simp_; //
        std::vector<size_t> indices_;
        SDMCSimCol(art::Ptr<SimParticle> simp,size_t index) : simp_(simp) {
          indices_.push_back(index);
        }
      };
      public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag Level"), 0};
        fhicl::Atom<int> print{ Name("printLevel"), Comment("Print Level"), 0};
        fhicl::Atom<int> pdg{ Name("pdgCode"), Comment("PDG code to use")};
        fhicl::Atom<int> proc{ Name("processCode"), Comment("Process code to use")};
        fhicl::Atom<size_t> mindigi{ Name("minStrawDigiMC"), Comment("Minimum number of StrawDigiMCs")};
        fhicl::Atom<art::InputTag> strawDigiMCs { Name("StrawDigiMCs"), Comment("StrawDigiMC Collection")};
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeKalSeedMC(const Parameters& config);

      private:
      void produce(art::Event& e) override;
      int debug_, diag_, printLevel_;
      int pdg_, proc_;
      size_t mindigi_;
      art::InputTag sdColTag_;
      void makeSeeds(SDMCSimCol const& sdmcsim, KalSeedCollection& kseeds, KalSeedMCCollection& kseedmcs);
    };
    MakeKalSeedMC::MakeKalSeedMC(const Parameters& config) :
      EDProducer(config),
      debug_(config().debug()),
      diag_(config().diag()),
      printLevel_(config().print()),
      pdg_(config().pdg()),
      proc_(config().proc()),
      mindigi_(config().mindigi()),
      sdColTag_(config().strawDigiMCs()){
        consumes<StrawDigiMCCollection>(sdColTag_);
        produces<KalSeedCollection>();
        produces<KalSeedMCCollection>();
      }


    void MakeKalSeedMC::produce(art::Event& event) {
      art::Handle<StrawDigiMCCollection> sdColHandle;
      event.getByLabel(sdColTag_, sdColHandle);
      const StrawDigiMCCollection& sdCol(*sdColHandle);
      unique_ptr<KalSeedCollection> kseeds(new KalSeedCollection);
      unique_ptr<KalSeedMCCollection> kseedmcs(new KalSeedMCCollection);

      // sort StrawDigiMCs by particle
      std::vector<SDMCSimCol> sdmcsims;
      for(size_t index=0;index < sdCol.size(); ++index){
        auto const& sdmc = sdCol[index];
        bool found(false);
        for(auto& sdmcsim : sdmcsims){
          if(sdmcsim.simp_ == sdmc.earlyStrawGasStep()->simParticle()){
            sdmcsim.indices_.push_back(index);
            found = true;
            break;
          }
        }
        if(!found)sdmcsims.emplace_back(sdmc.earlyStrawGasStep()->simParticle(),index);
      }
      // find collections of StrawDigiMCs that meet requirements
      for(auto const& sdmcsim : sdmcsims) {
        if(sdmcsim.simp_->pdgId() == pdg_ && sdmcsim.simp_->creationCode() == proc_  && sdmcsim.indices_.size() > mindigi_){
          makeSeeds(sdmcsim,*kseeds,*kseedmcs);
        }
      }
      event.put(move(kseeds));
      event.put(move(kseedmcs));

    }
    void MakeKalSeedMC::makeSeeds(SDMCSimCol const& sdmcsim, KalSeedCollection& kseeds, KalSeedMCCollection& kseedmcs) {

    }
  }

}
using mu2e::TrackerMC::MakeKalSeedMC;
DEFINE_ART_MODULE(MakeKalSeedMC);
