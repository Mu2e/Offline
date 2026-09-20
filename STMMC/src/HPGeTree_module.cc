// Adapted from ReadVirtualDetector_module.cc
// For StepPointMCs in STMDet, generates a TTree with energy in branch "E", time in branch "time",
// and true incoming energy in "incidentE".
//  - Iterate over the StepPointMCs, determine the associated SimParticle
//  - Trace the SimParticle lineage up until reaching a particle that originated outside the chosen detector volume
//  - Increment the energy associated with that incident particle and track its minimum interaction time
//  - Record the true start energy of the incident particle for purity calculations
//  - Print the PDG ID counts at the end of the job
// Input Parameters
//  - Detector - either "HPGe" or "LaBr" - applies a position cut to calculate the energy deposited by a particle going through the chosen detector
//  - StepPointMCsTag - tag of data product containing the StepPoints for STMDet
//  - SimParticlemvTag - tag of data product containing the SimParticles for STMDet

// stdlib includes
#include <limits>
#include <algorithm>
#include <tuple>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

// Mu2e type definitions
typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;

namespace mu2e {
  class HPGeTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<std::string> detector{ Name("Detector"), Comment("Which detector to generate energy histograms for, either 'HPGe' or 'LaBr'")};
        fhicl::Atom<art::InputTag> stepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> simParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit HPGeTree(const Parameters& conf);
      std::tuple<key_type, int, double> topParent(const SimParticle& particle);
      void analyze(const art::Event& event);
      void endJob();
    private:
      std::string detector = "";
      std::vector<std::string> detectors{"HPGe", "LaBr"};

      SimParticle stepParticle;

      key_type topParentId;
      std::vector<key_type> topParentIds;

      std::map<int, int> pdgIds; // <ID, count>
      std::map<key_type, double> EDeps, times, incidentEnergies; // Tracking variables

      TTree* ttree = nullptr;
      double xBeamCentre = -3904.0;
      int pdgId = 0;
      double E = 0.0, time = 0.0, incidentE = 0.0;

      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
  };

  HPGeTree::HPGeTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    detector(conf().detector()),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().simParticlemvTag())) {

      if (std::find(detectors.begin(), detectors.end(), detector) == detectors.end())
        throw cet::exception("Configuration") << "'detector' must be one of 'HPGe' or 'LaBr'";

      if (detector == "LaBr")
        throw cet::exception("Configuration") << "Currently this code only works for HPGe, exiting.\n";

      // Set up TTree
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "Detector ttree");
      ttree->Branch("E", &E, "E/D");
      ttree->Branch("time", &time, "time/D");
      ttree->Branch("incidentE", &incidentE, "incidentE/D"); // New branch for true origin energy
  };

  std::tuple<key_type, int, double> HPGeTree::topParent(const SimParticle& particle) {
    const SimParticle* current = &particle;

    // Trace lineage upward until we find a particle that originated OUTSIDE the detector volume
    while (current->parent().isNonnull()) {
      bool originatedInside = false;
      if (detector == "HPGe" && current->startPosition().x() <= xBeamCentre) {
        originatedInside = true;
      } else if (detector == "LaBr" && current->startPosition().x() >= xBeamCentre) {
        originatedInside = true;
      }

      // If the current particle originated outside the boundary, it is the incident parent
      if (!originatedInside) {
        break;
      }

      current = current->parent().get();
    }

    // Return the ID, PDG ID, and the true kinetic/total energy of the particle as it crossed the boundary
    return std::make_tuple(current->id(), current->pdgId(), current->startMomentum().e());
  };

  void HPGeTree::analyze(const art::Event& event) {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);

    if ((StepPointMCs.size() == 0) || (SimParticles.size() == 0))
      return;

    double parentEnergy = 0.0;

    for (const StepPointMC& step : StepPointMCs) {
      if ((detector == "HPGe") && (step.position().x() > xBeamCentre))
        continue;
      else if ((detector == "LaBr") && (step.position().x() < xBeamCentre))
        continue;

      stepParticle = SimParticles.at(step.trackId());
      std::tie(topParentId, pdgId, parentEnergy) = topParent(stepParticle);

      // Collate the data for this incident particle
      if (std::find(topParentIds.begin(), topParentIds.end(), topParentId) != topParentIds.end()) {
        EDeps[topParentId] += step.ionizingEdep();

        // Track the earliest hit time for the incident particle
        if (step.time() < times[topParentId]) {
          times[topParentId] = step.time();
        }
      }
      else {
        topParentIds.emplace_back(topParentId);
        EDeps.emplace(std::make_pair(topParentId, step.ionizingEdep()));
        times.emplace(std::make_pair(topParentId, step.time()));
        incidentEnergies.emplace(std::make_pair(topParentId, parentEnergy));

        // Log the incident PDG ID only once per parent particle
        if (pdgIds.find(pdgId) != pdgIds.end())
          pdgIds[pdgId] += 1;
        else
          pdgIds.emplace(std::make_pair(pdgId, 1));
      };
    };

    // Collect the data to the TTree
    for (size_t i = 0; i < topParentIds.size(); i++) {
      topParentId = topParentIds[i];
      E = EDeps[topParentId];
      time = times[topParentId];
      incidentE = incidentEnergies[topParentId];
      ttree->Fill();
    };

    EDeps.clear();
    times.clear();
    incidentEnergies.clear();
    topParentIds.clear();

    return;
  };

  void HPGeTree::endJob() {
    mf::LogInfo log("Detector tree");
    log << "==========Data summary==========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeTree)
