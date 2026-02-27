// Accept events with calo E(primary + nearby pileup) is above a given threshold
//  Original author; Michael MacKenzie, 2026

// Framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// cet
#include "cetlib_except/exception.h"

// Offline
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"

// C++
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <cmath>

namespace mu2e {

  class CaloPileupDtsFilter : public art::EDFilter {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> caloSteps         { Name("CaloShowerSteps")      , Comment("CaloShowerStep collections") };
      fhicl::Atom<art::InputTag>     primary           { Name("PrimaryParticle")      , Comment("PrimaryParticle tag") };
      fhicl::Atom<double>            minCaloStepEnergy { Name("MinimumStepEnergy")    , Comment("Minimum Calo step energy"), 0.0 };
      fhicl::Atom<double>            minPrimaryEnergy  { Name("MinimumPrimaryEnergy") , Comment("Minimum Primary Calo energy deposited"), 0.};
      fhicl::Atom<double>            maxTotalEnergy    { Name("MaximumTotalEnergy")   , Comment("Maximum Primary + Pileup energy")};
      fhicl::Atom<double>            minTotalEnergy    { Name("MinimumTotalEnergy")   , Comment("Minimum Primary + Pileup energy")};
      fhicl::Atom<double>            timeWindow        { Name("TimeWindow")           , Comment("Time window around the Primary deposit")};
      fhicl::Atom<double>            spaceWindow       { Name("SpaceWindow")          , Comment("Space window around the Primary deposit")};
      fhicl::Atom<bool>              nullFilter        { Name("NullFilter")           , Comment("If true, do not apply any filtering and accept all events"), false };
      fhicl::Atom<int>               diagLevel         { Name("DiagLevel")            , Comment("Diagnostic output level"), 0 };
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit CaloPileupDtsFilter(const Parameters& conf);
    bool filter(art::Event& event) override;
    bool beginRun(art::Run& aRun) override;
    void endJob() override;

  private:
    bool goodParticle(SimParticle const& simp) const;
    CLHEP::Hep3Vector stepPosition(const CaloShowerStep& step) const;
    bool checkSimParticle(const SimParticle& sim, const art::Event& event) const;

    // Input data
    std::vector<art::InputTag> caloStepsTags_;
    art::InputTag              primaryTag_;
    double                     minCaloStepEnergy_;
    double                     minPrimaryEnergy_;
    double                     maxTotalEnergy_;
    double                     minTotalEnergy_;
    double                     timeWindow_;
    double                     spaceWindow_;
    bool                       nullFilter_;
    int                        diagLevel_;

    // Data
    const Calorimeter*         calorimeter_;
    long unsigned              nEvents_;
    long unsigned              nPassed_;
  };

  //================================================================
  CaloPileupDtsFilter::CaloPileupDtsFilter(const Parameters& conf)
    : art::EDFilter{conf}
    , caloStepsTags_(conf().caloSteps())
    , primaryTag_(conf().primary())
    , minCaloStepEnergy_(conf().minCaloStepEnergy())
    , minPrimaryEnergy_(conf().minPrimaryEnergy())
    , maxTotalEnergy_(conf().maxTotalEnergy())
    , minTotalEnergy_(conf().minTotalEnergy())
    , timeWindow_(conf().timeWindow())
    , spaceWindow_(conf().spaceWindow())
    , nullFilter_(conf().nullFilter())
    , diagLevel_(conf().diagLevel())
    , calorimeter_(nullptr)
    , nEvents_(0)
    , nPassed_(0)
  {
    if(caloStepsTags_.size() == 0)         throw cet::exception("BADCONFIG") << "At least one CaloShowerStep collection must be specified\n";
    if(minPrimaryEnergy_ <= 0.)            throw cet::exception("BADCONFIG") << "Minimum primary energy must be > 0\n";
    if(maxTotalEnergy_ <= 0.)              throw cet::exception("BADCONFIG") << "Maximum total energy must be > 0\n";
    if(minTotalEnergy_ < 0.)               throw cet::exception("BADCONFIG") << "Minimum total energy must be >= 0\n";
    if(maxTotalEnergy_ <= minTotalEnergy_) throw cet::exception("BADCONFIG") << "Maximum total energy must be > minimum total energy\n";
    if(timeWindow_ < 0.)                   throw cet::exception("BADCONFIG") << "Time window must be >= 0\n";
    if(spaceWindow_ < 0.)                  throw cet::exception("BADCONFIG") << "Space window must be >= 0\n";

    // Declare data dependencies to avoid implicit consumption.
    for (auto const& tag : caloStepsTags_) {
      consumes<CaloShowerStepCollection>(tag);
    }
    consumes<PrimaryParticle>(primaryTag_);
  }

  //--------------------------------------------------------------------------------
  bool CaloPileupDtsFilter::goodParticle(SimParticle const& simp) const {
    // Check if the SimParticle is a valid particle to consider for pileup filtering
    return true; // For now, accept all particles
  }

  //--------------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloPileupDtsFilter::stepPosition(const CaloShowerStep& step) const {
    // Get the position of the CaloShowerStep in Mu2e coordinates
    const CLHEP::Hep3Vector& pos = step.position();
    const auto& calo_geom = calorimeter_->geomUtil();
    const auto pos_in_mu2e = calo_geom.crystalToMu2e(step.volumeG4ID(), pos);
    return pos_in_mu2e;
  }

  //--------------------------------------------------------------------------------
  bool CaloPileupDtsFilter::checkSimParticle(const SimParticle& sim, const art::Event& event) const {
    // for time rolling
    const double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

    if(diagLevel_ > 1) {
      std::cout << "[CaloPileupDtsFilter::"
                << moduleDescription().moduleLabel()
                << "::" << __func__ << "]"
                << " Checking SimParticle with PDG = " << sim.pdgId()
                << " E = " << sim.startMomentum().e()
                << " t = " << sim.startGlobalTime()
                << std::endl;
    }

    // First get the edep information for the primary particle
    double primary_edep(0.), primary_time(0.);
    CLHEP::Hep3Vector primary_pos(0.,0.,0.);
    for(const auto& tag : caloStepsTags_) {
      const auto handle = event.getValidHandle<CaloShowerStepCollection>(tag);
      for(const auto& css : *handle) {
        const auto sim_step = css.simParticle();
        if(sim_step.isNull()) continue;
        if(&(*sim_step) == &sim) {
          const double edep_step = css.energyDepBirks();
          if(diagLevel_ > 2) {
            std::cout << "  Primary particle step: PDG = " << sim.pdgId()
                      << " E = " << sim.startMomentum().e()
                      << " step E birks = " << edep_step
                      << " time = " << css.time()
                      << " volumeG4ID = " << css.volumeG4ID()
                      << " position = " << stepPosition(css)
                      << std::endl;
          }
          primary_edep += edep_step;
          primary_time += edep_step*std::fmod(css.time(), mbtime);
          primary_pos  += edep_step*stepPosition(css);
        }
      }
    }

    if(primary_edep > 0.) {
      // Get the average time and position
      primary_time /= primary_edep;
      primary_pos  /= primary_edep;
    } else {
      if(diagLevel_ > 1) std::cout << "  No CaloShowerSteps found for primary particle with PDG = " << sim.pdgId()
                                   << " E = " << sim.startMomentum().e()
                                   << std::endl;
      return false; // no energy deposit from the primary particle, reject
    }

    // Check if the primary particle deposit passes the minimum energy requirement
    if(primary_edep < minPrimaryEnergy_) {
      if(diagLevel_ > 1) std::cout << "  Primary particle with PDG = " << sim.pdgId()
                                   << " E = " << sim.startMomentum().e()
                                   << " has energy deposit = " << primary_edep << " MeV, below the minimum of "
                                   << minPrimaryEnergy_ << " MeV, rejecting\n";
      return false;
    }

    // Check if the primary already exceeds the maximum total energy (in which case we can reject without even looking for pileup)
    if(primary_edep > maxTotalEnergy_) {
      if(diagLevel_ > 1) std::cout << "  Primary particle with PDG = " << sim.pdgId()
                                   << " E = " << sim.startMomentum().e()
                                   << " has energy deposit = " << primary_edep << " MeV, above the maximum total energy of "
                                   << maxTotalEnergy_ << " MeV, rejecting\n";
      return false;
    }

    if(diagLevel_ > 1) {
      std::cout << "  Primary edep = " << primary_edep
                << " MeV, time = " << primary_time
                << " ns, position = " << primary_pos
                << std::endl;
    }

    // Next check for pileup deposits near the primary particle
    double pileup_edep = 0.;

    for(const auto& tag : caloStepsTags_) {
      const auto handle = event.getValidHandle<CaloShowerStepCollection>(tag);
      for(const auto& css : *handle) {
        const auto sim_step = css.simParticle();
        if(sim_step.isNull()) continue;
        if(&(*sim_step) == &sim) continue; // skip steps from the primary particle itself
        if(css.energyDepBirks() < minCaloStepEnergy_) continue; // skip low energy steps
        if(!goodParticle(*sim_step)) continue; // skip steps from particles we don't care about
        double time_diff = std::fabs(std::fmod(css.time(), mbtime) - primary_time);
        if(time_diff < timeWindow_) {
          auto step_pos   = stepPosition(css);
          double space_diff = (step_pos - primary_pos).mag();
          if(space_diff < spaceWindow_) {
            pileup_edep += css.energyDepBirks();
            if(diagLevel_ > 2) {
              std::cout << "  Pileup step: PDG = " << sim_step->pdgId()
                        << " E = " << sim_step->startMomentum().e()
                        << " Process = " << sim_step->creationCode()
                        << " step E birks = " << css.energyDepBirks()
                        << " time = " << css.time()
                        << " position = " << step_pos
                        << " time diff = " << time_diff
                        << " space diff = " << space_diff
                        << std::endl;
            }
          }
        }
      }
    }

    double total_edep = pileup_edep + primary_edep;
    const bool result = (total_edep >= minTotalEnergy_ && total_edep <= maxTotalEnergy_);
    if(diagLevel_ > 3 || (result && diagLevel_ > 0))
      std::cout << "[CaloPileupDtsFilter::"
                << moduleDescription().moduleLabel()
                << "::" << __func__ << "] "
                << event.id()
                << " : primary edep = " << primary_edep
                << " pileup edep = " << pileup_edep
                << " total edep = " << total_edep
                << " result = " << result
                << std::endl;
    return result;
  }

  //--------------------------------------------------------------------------------
  bool CaloPileupDtsFilter::filter(art::Event& event) {
    ++nEvents_;
    if(nullFilter_) {
      ++nPassed_;
      return true;
    }

    // Get the input data
    auto primaryH = event.getValidHandle<PrimaryParticle>(primaryTag_);
    const auto& primary = *primaryH;

    // Diagnostic printout
    if(diagLevel_ > 1) {
      std::cout << "[CaloPileupDtsFilter::"
                << moduleDescription().moduleLabel()
                << "::" << __func__ << "]"
                << " Event " << event.id()
                << " N(primaries) = " << primary.primarySimParticles().size()
                << std::endl;
      for(const auto& p : primary.primarySimParticles()) {
        std::cout << "  Primary particle: PDG = " << p->pdgId()
                  << " E = " << p->startMomentum().e()
                  << " t = " << p->startGlobalTime()
                  << std::endl;
      }
    }

    // Require a primary particle
    if(primary.primarySimParticles().empty()) {
      if(diagLevel_ > 0) std::cout << "[CaloPileupDtsFilter::"
                                   << moduleDescription().moduleLabel()
                                   << "::" << __func__ << "]"
                                   << " Event " << event.id()
                                   << " has no PrimaryParticle, rejecting\n";
      return false;
    }

    bool result = false;
    for(const auto& sim : primary.primarySimParticles()) {
      if(sim.isNull()) continue;
      result = checkSimParticle(*sim, event);
      if(result) break;
    }
    if(result) ++nPassed_;
    return result;
  }

  //--------------------------------------------------------------------------------
  bool CaloPileupDtsFilter::beginRun(art::Run& run) {
    // Get the calorimeter geometry
    GeomHandle<Calorimeter> ch;
    calorimeter_ = ch.get();
    return true;
  }

  //--------------------------------------------------------------------------------
  void CaloPileupDtsFilter::endJob() {
    mf::LogInfo("Summary")
      << "[CaloPileupDtsFilter::"
      << moduleDescription().moduleLabel()
      << "::" << __func__ << "]"
      << " Passed "
      << nPassed_ << " / " << nEvents_ << " events\n";
  }

}

DEFINE_ART_MODULE(mu2e::CaloPileupDtsFilter)
