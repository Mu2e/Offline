// Select cosmic events for mixing: this requires a 'signal-like' particle in the tracker
//  Original author; Dave Brown (LBNL) 2022

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include <map>
namespace mu2e {

  class CosmicMixingFilter : public art::EDFilter {
    public:

      struct TrkHitCount {
        unsigned nhits_;
        uint16_t minplane_, maxplane_;
        std::set<unsigned> planes_;
      };

      struct TimeCutConfig {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<double> maxTime { Name("MaximumTime"), Comment("Maximum time for good step (ns since POT)")};
        fhicl::Atom<double> minTime { Name("MinimumTime"), Comment("Minimum time for good step (ns since POT)")};
      };

      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Atom<art::InputTag> trkSteps { Name("StrawGasSteps"), Comment("StrawGasStep collection") };
        fhicl::Atom<int> diagLevel { Name("DiagLevel"), Comment("Diag printout level"),0 };
        fhicl::Atom<double> minTrkStepEnergy { Name("MinimumTrkStepEnergy"), Comment("Minimum Trk step energy"), 0.0 };
        fhicl::Atom<double> minPartMom { Name("MinimumPartMom"), Comment("Minimum particle momentum") };
        fhicl::Atom<double> maxPartMom { Name("MaximumPartMom"), Comment("Maximum particle momentum")};
        fhicl::Atom<unsigned> minNTrkSteps { Name("MinimumTrkSteps"), Comment("Minimum number of good tracker steps"), 10};
        fhicl::Atom<unsigned> minNTrkPlanes { Name("MinimumTrkPlanes"), Comment("Minimum number of good tracker planes")};
        fhicl::Atom<unsigned> minTrkPlaneSpan { Name("MinimumTrkPlaneSpan"), Comment("Minimum gap between good hittracker planes")};
        fhicl::Atom<double> maxImpact { Name("MaxImpact"), Comment("Maximum impact parameter, based on tracker info")};
        fhicl::OptionalTable<TimeCutConfig> timeCutConfig { fhicl::Name("TimeCutConfig") };

      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit CosmicMixingFilter(const Parameters& conf);
      bool filter(art::Event& event) override;
      void endJob() override;

    private:
      bool timeCut(double time) const; // particle time
      art::InputTag trkStepCol_;
      int diagLevel_;
      double minTrkE_, minPartM_, maxPartM_;
      unsigned minNTrk_, minNPlanes_, minPlaneSpan_;
      double maxImpact_;
      bool timecut_;
      double minTime_, maxTime_;
      unsigned nEvt_, nPassed_;
  };

  //================================================================
  CosmicMixingFilter::CosmicMixingFilter(const Parameters& conf) : art::EDFilter{conf}
  , trkStepCol_(conf().trkSteps())
    , diagLevel_(conf().diagLevel())
    , minTrkE_(conf().minTrkStepEnergy())
    , minPartM_(conf().minPartMom())
    , maxPartM_(conf().maxPartMom())
    , minNTrk_(conf().minNTrkSteps())
    , minNPlanes_(conf().minNTrkPlanes())
    , minPlaneSpan_(conf().minTrkPlaneSpan())
    , maxImpact_(conf().maxImpact())
    , timecut_(false)
    , minTime_(0.0), maxTime_(0.0)
    , nEvt_(0)
    , nPassed_(0)
    {
      consumes<StrawGasStepCollection>(trkStepCol_);
      TimeCutConfig tc;
      if(conf().timeCutConfig(tc)) {
        timecut_ = true;
        minTime_ = tc.minTime();
        maxTime_ = tc.maxTime();
      }
    }

  // input must be a physical time!
  bool CosmicMixingFilter::timeCut(double ptime) const { return (!timecut_) || (ptime > minTime_ && ptime < maxTime_); }

  bool CosmicMixingFilter::filter(art::Event& event) {
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double mbtime = accPar->deBuncherPeriod;
    GlobalConstantsHandle<ParticleDataList> pdt;
    GeomHandle<DetectorSystem> det;
    GeomHandle<BFieldManager> bfmgr;
    bool selecttrk(false);
    ++nEvt_;
    // Count Trk step from same particle
    using CT = std::map<const SimParticle*,TrkHitCount>;
    CT counttrk;
    auto sgscolH = event.getValidHandle<StrawGasStepCollection>(trkStepCol_);
    for(const auto& sgs : *sgscolH ) {
      double mom = sgs.momentum().R();
      if(sgs.ionizingEdep() > minTrkE_ &&
          mom > minPartM_ && mom < maxPartM_ &&
          timeCut(fmod(sgs.time(),mbtime))) {
        auto ifndhit = counttrk.find(sgs.simParticle().get());
        if(ifndhit == counttrk.end()){
          TrkHitCount thc;
          thc.nhits_ = 1;
          thc.planes_.insert(sgs.strawId().plane());
          thc.minplane_ = thc.maxplane_ = sgs.strawId().plane();
          counttrk.insert(CT::value_type(sgs.simParticle().get(),thc));
        } else {
          ifndhit->second.nhits_++;
          ifndhit->second.planes_.insert(sgs.strawId().plane());
          ifndhit->second.minplane_ = std::min(ifndhit->second.minplane_,sgs.strawId().plane());
          ifndhit->second.maxplane_ = std::max(ifndhit->second.maxplane_,sgs.strawId().plane());
        }
      }
    }
    for(auto itrk=counttrk.begin();itrk != counttrk.end();itrk++){
      if(diagLevel_ > 0)std::cout << "NHits " << itrk->second.nhits_ << " Nplanes " << itrk->second.planes_.size() << " plane span " << itrk->second.maxplane_-itrk->second.minplane_ << std::endl;
      if(itrk->second.nhits_ >= minNTrk_
          && itrk->second.planes_.size() >= minNPlanes_
          && unsigned(abs(itrk->second.maxplane_-itrk->second.minplane_)) >= minPlaneSpan_){
        // test the impact parameter.  Use info from the first sgs
        auto sp = itrk->first;
        double charge = pdt->particle(sp->pdgId()).charge();
        RobustHelix rh;
        for(const auto& sgs : *sgscolH ) {
          if(sgs.simParticle().get() == sp){
            CLHEP::Hep3Vector pos = sgs.position();
            CLHEP::Hep3Vector mom = sgs.momvec();
            double bz = bfmgr->getBField(det->toMu2e(pos)).z();
            TrkUtilities::RobustHelixFromMom(
                pos, mom,
                charge,
                bz,
                rh);
            if(diagLevel_ > 0)std::cout << "Impact Parameter " << fabs(rh.rcent() - rh.radius()) << std::endl;
            if(fabs(rh.rcent()-rh.radius()) < maxImpact_) {
              selecttrk = true;
              break;
            }
          }
        }
      }
      if(selecttrk)break;
    }
    if(selecttrk)nPassed_++;
    return selecttrk;
  }

  void CosmicMixingFilter::endJob() {
    mf::LogInfo("Summary")
      <<"CosmicMixingFilter_module: passed "
      <<nPassed_<<" / "<<nEvt_<<" events\n";
  }

}
DEFINE_ART_MODULE(mu2e::CosmicMixingFilter)
