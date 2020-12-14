// Pass events when particle of specified types and energies create min detector signals (Calo, Trk, or CRV)
//  Original author; Dave Brown (LBNL) 2020

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include <map>
namespace mu2e {

  class DetectorStepFilter : public art::EDFilter {
    public:
      struct Config {
	using Name=fhicl::Name;
	using Comment=fhicl::Comment;

	fhicl::Sequence<art::InputTag> trkSteps { Name("StrawGasSteps"), Comment("StrawGasStep collections") };
	fhicl::Sequence<art::InputTag> caloSteps { Name("CaloShowerSteps"), Comment("CaloShowerStep collections") };
	fhicl::Sequence<art::InputTag> crvSteps { Name("CrvSteps"), Comment("Crv StepPointMC collections") };

	fhicl::Atom<double> minTrkStepEnergy { Name("MinimumTrkStepEnergy"), Comment("Minimum Trk step energy"), -std::numeric_limits<double>::max() };
	fhicl::Atom<double> minCaloStepEnergy { Name("MinimumCaloStepEnergy"), Comment("Minimum Calo step energy"), -std::numeric_limits<double>::max() };
	fhicl::Atom<double> minCrvStepEnergy { Name("MinimumCrvStepEnergy"), Comment("Minimum Crv step energy"), -std::numeric_limits<double>::max() };

	fhicl::Atom<double> minPartMom { Name("MinimumPartMom"), Comment("Minimum particle momentum"), -std::numeric_limits<double>::max() };
	fhicl::Atom<double> maxPartMom { Name("MaximumPartMom"), Comment("Maximum particle momentum"), std::numeric_limits<double>::max() };

	fhicl::Atom<bool> testTrk { Name("TestTrk"), Comment("Test Tracker steps"), true};
	fhicl::Atom<bool> testCalo { Name("TestCalo"), Comment("Test calorimeter steps"), true };
	fhicl::Atom<bool> testCrv { Name("TestCrv"), Comment("Test CRV steps"), false };

	fhicl::Atom<unsigned> minNTrkSteps { Name("MinimumTrkSteps"), Comment("Minimum number of good tracker steps"), 10};
	fhicl::Atom<double> minSumCaloStepE { Name("MinimumSumCaloStepE"), Comment("Minimum E sum of good calorimeter steps (MeV) "), 50.0 };
	fhicl::Atom<unsigned> minNCrvSteps { Name("MinimumCrvSteps"), Comment("Minimum number of good CRV steps"), 3 };

	fhicl::Sequence<int> keepPDG { Name("KeepPDG"), Comment("PDG particle codes to keep") };

      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit DetectorStepFilter(const Parameters& conf);
      bool filter(art::Event& event) override;
      void endJob() override;

    private:
      bool goodParticle(SimParticle const& simp) const;
      double minTrkE_, minCaloE_, minCrvE_;
      double minPartM_, maxPartM_;
      bool testTrk_, testCalo_, testCrv_;
      std::vector<PDGCode::type> pdgToKeep_;
      std::vector<art::InputTag> trkStepCols_, caloStepCols_, crvStepCols_;
      unsigned minNTrk_, minNCrv_;
      double minSumCaloE_;
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
    , testTrk_(conf().testTrk())
    , testCalo_(conf().testCalo())
    , testCrv_(conf().testCrv())
    , minNTrk_(conf().minNTrkSteps())
    , minNCrv_(conf().minNCrvSteps())
    , minSumCaloE_(conf().minSumCaloStepE())
    , nEvt_(0)
    , nPassed_(0)
    {
      for(const auto ikeep : conf().keepPDG()) { pdgToKeep_.emplace_back(PDGCode::type(ikeep)); }
      for(const auto& trktag : conf().trkSteps()) { trkStepCols_.emplace_back(trktag); }
      for(const auto& calotag : conf().caloSteps()) { caloStepCols_.emplace_back(calotag); }
      for(const auto& crvtag : conf().crvSteps()) { crvStepCols_.emplace_back(crvtag); }
    }

  bool DetectorStepFilter::filter(art::Event& event) {
    bool retval(false);
    ++nEvt_;
    // Count Trk step from same particle
    if(testTrk_ ) {
      using CT = std::map<const SimParticle*,unsigned>;
      CT counttrk;
      for(const auto& trkcoltag : trkStepCols_) {
	auto sgscolH = event.getValidHandle<StrawGasStepCollection>(trkcoltag);
	for(const auto& sgs : *sgscolH ) {
	  double mom = sgs.momentum().R();
	  if(sgs.ionizingEdep() > minTrkE_ && 
	      mom > minPartM_ && mom < maxPartM_ &&
	      goodParticle(*sgs.simParticle())){
	    auto ifnd = counttrk.find(sgs.simParticle().get());
	    if(ifnd == counttrk.end())
	      counttrk.insert(CT::value_type(sgs.simParticle().get(),1));
	    else
	      ifnd->second++;
	  }
	}
      }
      for(auto itrk=counttrk.begin();itrk != counttrk.end();itrk++){
	if(itrk->second >= minNTrk_){
	  retval = true;
	  break;
	}
      }
    }
    // sum Calo energy from same particle
    if(testCalo_ && !retval){
      using CES = std::map<const SimParticle*,float>;
      CES caloESum;
      for(const auto& calocoltag : caloStepCols_) {
	auto csscolH = event.getValidHandle<CaloShowerStepCollection>(calocoltag);
	for(const auto& css : *csscolH ) {
	  if(css.energyDepBirks() > minCaloE_ && 
	      css.momentumIn() > minPartM_ && css.momentumIn() < maxPartM_ &&
	      goodParticle(*css.simParticle())){
	    auto ifnd = caloESum.find(css.simParticle().get());
	    if(ifnd == caloESum.end())
	      caloESum.insert(CES::value_type(css.simParticle().get(),css.energyDepBirks()));
	    else
	      ifnd->second += css.energyDepBirks(); // not sure Birks or G4 FIXME!
	  }
	}
	for(auto icalo=caloESum.begin();icalo != caloESum.end();icalo++){
	  if(icalo->second >= minSumCaloE_){
	    retval = true;
	    break;
	  }
	}
      }
    }
    // Count good CRV from same particle; update this when CrvStep is written FIXME!!!
    if(testCrv_ && !retval ){
      using CC = std::map<const SimParticle*,unsigned>;
      CC countcrv;
      for(const auto& crvcoltag : crvStepCols_) {
	auto spmccolH = event.getValidHandle<StepPointMCCollection>(crvcoltag);
	for(const auto& spmc : *spmccolH ) {
	  double mom = spmc.momentum().mag();
	  if(spmc.ionizingEdep() > minCrvE_ && 
	      mom > minPartM_ && mom < maxPartM_ &&
	      goodParticle(*spmc.simParticle())) {
	    auto ifnd = countcrv.find(spmc.simParticle().get());
	    if(ifnd == countcrv.end())
	      countcrv.insert(CC::value_type(spmc.simParticle().get(),1));
	    else
	      ifnd->second++;
	  }
	}
      }
      for(auto icrv=countcrv.begin();icrv != countcrv.end();icrv++){
	if(icrv->second >= minNCrv_){
	  retval = true;
	  break;
	}
      }
    }

    if(retval) ++nPassed_;
    return retval;
  }

  void DetectorStepFilter::endJob() {
    mf::LogInfo("Summary")
      <<"DetectorStepFilter_module: passed "
      <<nPassed_<<" / "<<nEvt_<<" events\n";
  }

  bool DetectorStepFilter::goodParticle(SimParticle const& simp) const {
    bool retval(false);
    for(auto pdg: pdgToKeep_){
      if(pdg == simp.pdgId()){
	retval = true;
	break;
      }
    }
    return retval;
  }

}

DEFINE_ART_MODULE(mu2e::DetectorStepFilter);
