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
// utiliities
// persistent data
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include <vector>
#include <limits>
using namespace std;
namespace mu2e {
  namespace TrackerMC {
    using MCTrajIter = MCTrajectoryCollection::const_iterator;
    class MakeMCKalSeed : public art::EDProducer {
      struct SDMCSimCol {
        MCTrajIter mctraj_;
        auto const& simParticle() const { return mctraj_->first; }
        auto const& mcTraj() const { return mctraj_->second; }
        std::vector<size_t> indices_;
        SDMCSimCol(MCTrajIter const& mctraj,size_t index) : mctraj_(mctraj) {
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
        fhicl::Atom<art::InputTag> mcTrajectoryTag{Name("mcTrajectoryTag"), Comment("InputTag for the MCTrajectory")};
        fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream") };
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeMCKalSeed(const Parameters& config);

      private:
      void produce(art::Event& e) override;
      int debug_, diag_, printLevel_;
      int pdg_, proc_;
      size_t mindigi_;
      art::InputTag sdColTag_, mcTrajColTag_;
      TrkFitDirection tdir_;
      // helper functions
      KalSeed makeKalSeed(SDMCSimCol const& sdmcsim) const;
      double MCT0(MCTrajectory const& mctraj) const;
   };

    MakeMCKalSeed::MakeMCKalSeed(const Parameters& config) :
      EDProducer(config),
      debug_(config().debug()),
      diag_(config().diag()),
      printLevel_(config().print()),
      pdg_(config().pdg()),
      proc_(config().proc()),
      mindigi_(config().mindigi()),
      sdColTag_(config().strawDigiMCs()),
      mcTrajColTag_(config().mcTrajectoryTag()),
      tdir_(TrkFitDirection::FitDirection(config().fitDirection()))
    {
      consumes<StrawDigiMCCollection>(sdColTag_);
      consumes<MCTrajectoryCollection>(mcTrajColTag_);
      produces<KalSeedCollection>();
    }

    void MakeMCKalSeed::produce(art::Event& event) {
      art::Handle<StrawDigiMCCollection> sdColHandle;
      event.getByLabel(sdColTag_, sdColHandle);
      const StrawDigiMCCollection& sdCol(*sdColHandle);

      art::Handle<MCTrajectoryCollection> mcTrajColHandle;
      event.getByLabel(mcTrajColTag_, mcTrajColHandle);
      const MCTrajectoryCollection& mcTrajCol(*mcTrajColHandle);

      // create output
      unique_ptr<KalSeedCollection> kseeds(new KalSeedCollection);

      // sort StrawDigiMCs by particle
      std::vector<SDMCSimCol> sdmcsims;
      for(size_t index=0;index < sdCol.size(); ++index){
        auto const& sdmc = sdCol[index];
        bool found(false);
        for(auto& sdmcsim : sdmcsims){
          if(sdmcsim.simParticle() == sdmc.earlyStrawGasStep()->simParticle()){
            sdmcsim.indices_.push_back(index);
            found = true;
            break;
          }
        }
        if(!found){
          auto mctraj = mcTrajCol.find(sdmc.earlyStrawGasStep()->simParticle());
          if(mctraj != mcTrajCol.end()) sdmcsims.emplace_back(mctraj,index);
        }
      }
      // find collections of StrawDigiMCs that meet requirements
      for(auto const& sdmcsim : sdmcsims) {
        if(sdmcsim.simParticle()->pdgId() == pdg_ && sdmcsim.simParticle()->creationCode() == proc_  && sdmcsim.indices_.size() > mindigi_){
          // should check particle direction for consistency with tdir TODO
          kseeds->push_back(makeKalSeed(sdmcsim));
        }
      }
      event.put(move(kseeds));
    }

     KalSeed MakeMCKalSeed::makeKalSeed(SDMCSimCol const& sdmcsim) const {
       // create seed with basic information
      KalSeed kseed(sdmcsim.simParticle()->pdgId(),TrkFitFlag(TrkFitFlag::MCSeed));
      // fill in segment, hit, and straw information TODO
      return kseed;
    }

    double MakeMCKalSeed::MCT0(MCTrajectory const& mctraj) const {
      // find the t0 from the MCTrajectory
      auto zminpoint = mctraj.points().begin();
      double zmin = fabs(zminpoint->z());
      for(auto ipoint = mctraj.points().begin(); ipoint != mctraj.points().end(); ++ipoint) {
        if(fabs(ipoint->z()) < zmin){
          zminpoint = ipoint;
          zmin = fabs(ipoint->z());
        }
      }
      // interpolate or extrapolate to find t0.
      double t0;
      if(zminpoint->z() < 0.0){
        auto nextpoint = zminpoint; ++nextpoint;
        if(nextpoint != mctraj.points().end()){
          t0 = zminpoint->t() - zminpoint->z()*(nextpoint->t()-zminpoint->t())/(nextpoint->z()-zminpoint->z());
        } else {
          auto prevpoint = zminpoint; --prevpoint;
          t0 = zminpoint->t() - zminpoint->z()*(zminpoint->t()-prevpoint->t())/(zminpoint->z()-prevpoint->z());
        }
      } else {
        if(zminpoint != mctraj.points().begin()){
          auto prevpoint = zminpoint; --prevpoint;
          t0 = zminpoint->t() - zminpoint->z()*(zminpoint->t()-prevpoint->t())/(zminpoint->z()-prevpoint->z());
        } else {
          auto nextpoint = zminpoint; ++nextpoint;
          t0 = zminpoint->t() - zminpoint->z()*(nextpoint->t()-zminpoint->t())/(nextpoint->z()-zminpoint->z());
        }
      }
      return t0;
    }
  }
}
using mu2e::TrackerMC::MakeMCKalSeed;
DEFINE_ART_MODULE(MakeMCKalSeed)
