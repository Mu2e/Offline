// Pass events when particle of specified types and energies create min detector signals (Calo, Trk, or CRV)
//  Original author; Dave Brown (LBNL) 2020

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
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include <map>
namespace mu2e {

  class DetectorStepFilter : public art::EDFilter {
    public:
      static constexpr double maxE_ = 1.0e6; // 1 TeV
      static constexpr unsigned maxN_ = 1000000;

      struct TimeCutConfig {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<double> maxTime { Name("MaximumTime"), Comment("Maximum time for good step (ns since POT)")};
        fhicl::Atom<double> minTime { Name("MinimumTime"), Comment("Minimum time for good step (ns since POT)")};
      };

      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Sequence<art::InputTag> trkSteps { Name("StrawGasSteps"), Comment("StrawGasStep collections") };
        fhicl::Sequence<art::InputTag> caloSteps { Name("CaloShowerSteps"), Comment("CaloShowerStep collections") };
        fhicl::Sequence<art::InputTag> crvSteps { Name("CrvSteps"), Comment("CrvStep collections") };

        fhicl::Atom<double> minTrkStepEnergy { Name("MinimumTrkStepEnergy"), Comment("Minimum Trk step energy"), 0.0 };
        fhicl::Atom<double> minCaloStepEnergy { Name("MinimumCaloStepEnergy"), Comment("Minimum Calo step energy"), 0.0 };
        fhicl::Atom<double> minCrvStepEnergy { Name("MinimumCrvStepEnergy"), Comment("Minimum Crv step energy"), 0.0 };

        fhicl::Atom<double> minPartMom { Name("MinimumPartMom"), Comment("Minimum particle momentum"), 0.0 };
        fhicl::Atom<double> maxPartMom { Name("MaximumPartMom"), Comment("Maximum particle momentum"), maxE_ };

        fhicl::Atom<bool> orRequirements { Name("ORRequirements"), Comment("Take the logical OR of requirements, otherwise take the AND"), true };

        fhicl::Atom<unsigned> minNTrkSteps { Name("MinimumTrkSteps"), Comment("Minimum number of good tracker steps"), 10};
        fhicl::Atom<double> minSumCaloStepE { Name("MinimumSumCaloStepE"), Comment("Minimum E sum of good calorimeter steps (MeV) "), 50.0 };
        fhicl::Atom<unsigned> minNCrvSteps { Name("MinimumCrvSteps"), Comment("Minimum number of good CRV steps"), 3 };

        fhicl::Atom<unsigned> maxNTrkSteps { Name("MaximumTrkSteps"), Comment("Maximum number of good tracker steps"), maxN_ };
        fhicl::Atom<double> maxSumCaloStepE { Name("MaximumSumCaloStepE"), Comment("Maximum E sum of good calorimeter steps (MeV) "),  maxE_ };
        fhicl::Atom<unsigned> maxNCrvSteps { Name("MaximumCrvSteps"), Comment("Maximum number of good CRV steps"), maxN_ };

        fhicl::Sequence<int> keepPDG { Name("KeepPDG"), Comment("PDG particle codes to keep") };

        fhicl::OptionalTable<TimeCutConfig> timeCutConfig { fhicl::Name("TimeCutConfig") };

      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit DetectorStepFilter(const Parameters& conf);
      bool filter(art::Event& event) override;
      void endJob() override;

    private:
      bool goodParticle(SimParticle const& simp) const; // select particles whose steps to count
      bool timeCut(double time) const; // particle time
      double minTrkE_, minCaloE_, minCrvE_;
      double minPartM_, maxPartM_;
      bool or_;
      std::vector<PDGCode::type> pdgToKeep_;
      std::vector<art::InputTag> trkStepCols_, caloStepCols_, crvStepCols_;
      unsigned minNTrk_, minNCrv_;
      double minSumCaloE_;
      unsigned maxNTrk_, maxNCrv_;
      double maxSumCaloE_;
      bool timecut_;
      double minTime_, maxTime_;
      unsigned nEvt_, nPassed_;
  };

  //================================================================
  DetectorStepFilter::DetectorStepFilter(const Parameters& conf)
    : art::EDFilter{conf}
  , minTrkE_(conf().minTrkStepEnergy())
    , minCaloE_(conf().minCaloStepEnergy())
    , minCrvE_(conf().minCrvStepEnergy())
    , minPartM_(conf().minPartMom())
    , maxPartM_(conf().maxPartMom())
    , or_(conf().orRequirements())
    , minNTrk_(conf().minNTrkSteps())
    , minNCrv_(conf().minNCrvSteps())
    , minSumCaloE_(conf().minSumCaloStepE())
    , maxNTrk_(conf().maxNTrkSteps())
    , maxNCrv_(conf().maxNCrvSteps())
    , maxSumCaloE_(conf().maxSumCaloStepE())
    , timecut_(false)
    , minTime_(0.0), maxTime_(0.0)
    , nEvt_(0)
    , nPassed_(0)
    {
      for(const auto ikeep : conf().keepPDG()) { pdgToKeep_.emplace_back(PDGCode::type(ikeep)); }
      for(const auto& trktag : conf().trkSteps()) { trkStepCols_.emplace_back(trktag); consumes<StrawGasStepCollection>(trktag); }
      for(const auto& calotag : conf().caloSteps()) { caloStepCols_.emplace_back(calotag); consumes<CaloShowerStepCollection>(calotag); }
      for(const auto& crvtag : conf().crvSteps()) { crvStepCols_.emplace_back(crvtag);  consumes<CrvStepCollection>(crvtag); }
      TimeCutConfig tc;
      if(conf().timeCutConfig(tc)) {
         timecut_ = true;
         minTime_ = tc.minTime();
         maxTime_ = tc.maxTime();
      }
    }

  // input must be a physical time!
  bool DetectorStepFilter::timeCut(double ptime) const { return (!timecut_) || (ptime > minTime_ || ptime < maxTime_); } // maxtime is 1 cycle around!!!

  bool DetectorStepFilter::filter(art::Event& event) {
    double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    bool selecttrk(false), selectcalo(false), selectcrv(false);
    ++nEvt_;
    // Count Trk step from same particle
    using CT = std::map<const SimParticle*,unsigned>;
    for(const auto& trkcoltag : trkStepCols_) {
      CT counttrk;
      auto sgscolH = event.getValidHandle<StrawGasStepCollection>(trkcoltag);
      for(const auto& sgs : *sgscolH ) {
        double mom = sgs.momentum().R();
        if(sgs.ionizingEdep() > minTrkE_ &&
            mom > minPartM_ && mom < maxPartM_ &&
            goodParticle(*sgs.simParticle()) &&
            timeCut(fmod(sgs.time(),mbtime))) {
          auto ifnd = counttrk.find(sgs.simParticle().get());
          if(ifnd == counttrk.end())
            counttrk.insert(CT::value_type(sgs.simParticle().get(),1));
          else
            ifnd->second++;
        }
      }
      for(auto itrk=counttrk.begin();itrk != counttrk.end();itrk++){
        if(itrk->second >= minNTrk_ && itrk->second <= maxNTrk_ ){
          selecttrk = true;
          break;
        }
      }
      if(selecttrk)break;
    }
    // sum Calo energy from same particle
    using CES = std::map<const SimParticle*,float>;
    for(const auto& calocoltag : caloStepCols_) {
      CES caloESum;
      auto csscolH = event.getValidHandle<CaloShowerStepCollection>(calocoltag);
      for(const auto& css : *csscolH ) {
        if(css.energyDepBirks() > minCaloE_ &&
            css.momentumIn() > minPartM_ && css.momentumIn() < maxPartM_ &&
            goodParticle(*css.simParticle()) &&
            timeCut(fmod(css.time(),mbtime))) {
          auto ifnd = caloESum.find(css.simParticle().get());
          if(ifnd == caloESum.end())
            caloESum.insert(CES::value_type(css.simParticle().get(),css.energyDepBirks()));
          else
            ifnd->second += css.energyDepBirks();
        }
      }
      for(auto icalo=caloESum.begin();icalo != caloESum.end();icalo++){
        if(icalo->second >= minSumCaloE_ && icalo->second <= maxSumCaloE_ ){
          selectcalo = true;
          break;
        }
      }
      if(selectcalo)break;
    }
    // Count good CRV from same particle
    using CC = std::map<const SimParticle*,unsigned>;
    for(const auto& crvcoltag : crvStepCols_) {
      CC countcrv;
      auto crvscolH = event.getValidHandle<CrvStepCollection>(crvcoltag);
      for(const auto& crvs : *crvscolH ) {
        double mom = crvs.startMom().R();
        if(crvs.visibleEDep() > minCrvE_ &&
            mom > minPartM_ && mom < maxPartM_ &&
            goodParticle(*crvs.simParticle()) &&
            timeCut(fmod(crvs.startTime(),mbtime))) {
          auto ifnd = countcrv.find(crvs.simParticle().get());
          if(ifnd == countcrv.end())
            countcrv.insert(CC::value_type(crvs.simParticle().get(),1));
          else
            ifnd->second++;
        }
      }
      for(auto icrv=countcrv.begin();icrv != countcrv.end();icrv++){
        if(icrv->second >= minNCrv_ && icrv->second <= maxNCrv_){
          selectcrv = true;
          break;
        }
      }
      if(selectcrv)break;
    }
    // global selection
    bool retval =( (or_ && (selecttrk || selectcalo || selectcrv)) ||
        ((!or_) && ( selecttrk && selectcalo && selectcrv)) );
    if(retval)nPassed_++;
    return retval;
  }

  void DetectorStepFilter::endJob() {
    mf::LogInfo("Summary")
      <<"DetectorStepFilter_module: passed "
      <<nPassed_<<" / "<<nEvt_<<" events\n";
  }

  bool DetectorStepFilter::goodParticle(SimParticle const& simp) const {
    bool retval = pdgToKeep_.size() > 0 ? false : true;
    for(auto pdg: pdgToKeep_){
      if(pdg == simp.pdgId()){
        retval = true;
        break;
      }
    }
    return retval;
  }

}

DEFINE_ART_MODULE(mu2e::DetectorStepFilter)
