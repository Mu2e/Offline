// Adapted from ReadVirtualDetector_module.cc
// For StepPointMCs in STMDet, generates a TTree with energy in branch "E" and time in branch "time". For a chosen STM detector
//  - Iterate over the StepPointMCs, determine the associated SimParticle
//  - Find the most parent particle of the associated SimParticle within the associated STMDet
//  - Increment the energy associated with the parent particle
//    - If the parent particle ID doesn't exist in the collection, add a new entry to the parent particle ID vector, energy deposited map, and time map
//  - Determine the PDG ID of the top particles, increment the mapped counter
//    - If the PDG ID entry does not exist in the map, create a new entry.
//  - Print the PDG ID counts at the end of the job
// Input Parameters
//  - Detector - either "HPGe" or "LaBr" - applies a position cut to calculate the energy deposited by a particle going through the chosen detector
//  - StepPointMCsTag - tag of data product containing the StepPoints for STMDet
//  - SimParticlemvTag - tag of data product containing the SimParticles for STMDet
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <limits>

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
      std::tuple<key_type, int> topParent(std::set<key_type>& SimParticleIDs, const SimParticle particle);
      double parentTime(const art::Event& event, key_type parentId);
      void analyze(const art::Event& event);
      void endJob();
    private:
      std::string detector = "";
      std::vector<std::string> detectors{"HPGe", "LaBr"};

      SimParticle stepParticle;
      std::set<key_type> SimParticleIds;

      key_type topParentId;
      std::vector<key_type> topParentIds;
      std::vector<key_type>::iterator topParentIdsIt;

      std::map<int, int> pdgIds; // <ID, count>
      std::map<int, int>::iterator pdgIdsIt;

      std::map<key_type, double> EDeps, times; // EDeps = <ID, deposited energy>, times = <ID, time>

      TTree* ttree = nullptr;
      double xBeamCentre = -3904.0;
      int pdgId = 0;
      double E = 0.0, time = 0.0;

      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      art::Ptr<SimParticle> parent;
  };

  HPGeTree::HPGeTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    detector(conf().detector()),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().simParticlemvTag())) {
      // Check if detector is one of the allowed types
      if (std::find(detectors.begin(), detectors.end(), detector) == detectors.end())
        throw cet::exception("Configuration") << "'detector' must be one of 'HPGe' or 'LaBr'";
      // TODO - remove this when the LaBr shower energy is defined
      if (detector == "LaBr")
        throw cet::exception("Configuration") << "Currently this code only works for HPGe, exiting.\n";

      // Set up TTree
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "Detector ttree");
      ttree->Branch("E", &E, "E/D");
      ttree->Branch("time", &time, "time/D");
  };

  std::tuple<key_type, int> HPGeTree::topParent(std::set<key_type>& SimParticleIds, const SimParticle particle) {
    // Get the particle parent
    parent = particle.parent();

    // If the passed particle has no parent in STMDet, return its ID and PDG ID
    if (std::find(SimParticleIds.begin(), SimParticleIds.end(), parent->id()) == SimParticleIds.end())
      return std::make_tuple(particle.id(), particle.pdgId());

    // If the particle has a parent in STMDet, update the particle parent
    while (std::find(SimParticleIds.begin(), SimParticleIds.end(), parent->parent()->id()) != SimParticleIds.end())
      parent = parent->parent();

    // Return the ID and PDG ID
    return std::make_tuple(parent->id(), parent->pdgId());
  };

  double HPGeTree::parentTime(const art::Event& event, key_type parentId) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);

    // Set up a variable to track the maximum particle time
    time = std::numeric_limits<double>::max();

    // Loop over the StepPointMCs
    for (const StepPointMC& Step : StepPointMCs) {
      // Get the step time
      stepParticle = SimParticles.at(Step.trackId());
      // If the step ID (track ID) is equal to the topmost parent and its time is lower than those previously selected, update the time
      if ((stepParticle.id() == parentId) && (Step.time() < time))
          time = Step.time();
    };

    // If the time has not been updated, throw, otherwise return the updated time
    if (time > std::numeric_limits<double>::max() * 0.99)
      throw cet::exception("LogicError") << "Time has not been updated! \n";
    return time;
  };

  void HPGeTree::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);

    // Validate that these data products exist
    if ((StepPointMCs.size() == 0) || (SimParticles.size() == 0))
      return;

    // Collect the particle IDs for StepPointMCs in STMDet
    // Note - IDs are from the StepPoints as the SimParticle contains the full geneaology, this is required to keep the volume constrianed to STMDet
    for (const StepPointMC& step : StepPointMCs)
      SimParticleIds.insert(SimParticles.at(step.trackId()).id());

    // Loop over all steps
    for (const StepPointMC& step : StepPointMCs) {
      // Select the appropriate steps
      if ((detector == "HPGe") && (step.position().x() > xBeamCentre))
        continue;
      else if ((detector == "LaBr") && (step.position().x() < xBeamCentre))
        continue;

      // Get the associated top particle
      stepParticle = SimParticles.at(step.trackId());
      std::tie(topParentId, pdgId) = topParent(SimParticleIds, stepParticle);
      // Precautionary check if the SimParticle ID is in STMDet
      if (std::find(SimParticleIds.begin(), SimParticleIds.end(), topParentId) == SimParticleIds.end())
        throw cet::exception("LogicError") << "The found parent ID is not a member of the data product\n";

      // Collate the data
      topParentIdsIt = std::find(topParentIds.begin(), topParentIds.end(), topParentId);
      if (topParentIdsIt != topParentIds.end()) {
        if (detector == "HPGe")
          EDeps[topParentId] += step.ionizingEdep();
        // TODO - include the LaBr response here
      }
      else {
        topParentIds.emplace_back(topParentId);
        if (detector == "HPGe")
          EDeps.emplace(std::make_pair(topParentId, step.ionizingEdep()));
        // TODO - include the LaBr response here
        times.emplace(std::make_pair(topParentId, parentTime(event, topParentId)));
      };

      // Generate the data summary
      pdgIdsIt = pdgIds.find(pdgId);
      if (pdgIdsIt != pdgIds.end())
        pdgIds[pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(pdgId, 1));
    }; // end for step

    // Collect the data to the TTree
    for (size_t i = 0; i < topParentIds.size(); i++) {
      topParentId = topParentIds[i];
      E = EDeps[topParentId];
      time = times[topParentId];
      ttree->Fill();
    };

    // Set up the data products for collection from the next event
    SimParticleIds.clear();
    EDeps.clear();
    times.clear();
    topParentIds.clear();

    return;
  }; // end analyze

  void HPGeTree::endJob() {
    mf::LogInfo log("Detector tree");
    log << "==========Data summary==========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeTree)
