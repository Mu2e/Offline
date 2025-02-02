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

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// exception handling
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "TTree.h"

// stdlib includes
#include <iostream>

typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;

namespace mu2e {
  class EDepFromShower : public art::EDAnalyzer {
    // Initialization
    std::string Detector = "";
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    art::ProductToken<SimParticleCollection> SimParticlemvToken;

    // Set up
    std::vector<std::string> Detectors{"HPGe", "LaBr"};
    TTree* ttree = nullptr;

    // Data collection
    int pdgId = 0;
    double E = 0.0, time = 0.0;
    art::Ptr<SimParticle> parent;

    // Data selection
    double xBeamCentre = -3904.0;
    key_type topParentId;
    SimParticle stepParticle;
    std::set<key_type> SimParticleIds;
    std::vector<key_type> topParentIds;
    std::vector<key_type>::iterator topParentIdsIt;

    // Data storage
    std::map<int, int> pdgIds; // <ID, count>
    std::map<int, int>::iterator pdgIdsIt;
    std::map<key_type, double> EDeps, times; // EDeps = <ID, deposited energy>, times = <ID, time>

  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<std::string> Detector{ Name("Detector"), Comment("Which detector to generate energy histograms for, either 'HPGe' or 'LaBr'")};
      fhicl::Atom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
      fhicl::Atom<art::InputTag> SimParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit EDepFromShower(const Parameters& conf);
    std::tuple<key_type, int> topParent(std::set<key_type>& SimParticleIDs, const SimParticle particle);
    double parentTime(const art::Event& event, key_type parentId);
    void analyze(const art::Event& event);
    void endJob();
  };

  EDepFromShower::EDepFromShower(const Parameters& conf) :
    art::EDAnalyzer(conf),
    Detector(conf().Detector()),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag()))
  {
    // TODO - remove this when the LaBr shower energy is defined
    if (Detector == "LaBr")
      throw cet::exception("Configuration") << "Currently this code only works for LaBr, exiting.\n";

    // Check if "Detector" is one of the allowed types
    if (std::find(Detectors.begin(), Detectors.end(), Detector) == Detectors.end())
      throw cet::exception("Configuration") << "'Detector' must be one of 'HPGe' or 'LaBr'";

    // Set up TTree
    art::ServiceHandle<art::TFileService> tfs;
    ttree = tfs->make<TTree>( "ttree", "Detector ttree");
    ttree->Branch("E", &E, "E/D");
    ttree->Branch("time", &time, "time/D");
  };

  std::tuple<key_type, int> EDepFromShower::topParent(std::set<key_type>& SimParticleIds, const SimParticle particle) {
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

  double EDepFromShower::parentTime(const art::Event& event, key_type parentId) {
    // Get the data products from the event
    auto const& Steps = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);

    // Set up a variable to track the maximum particle time
    time = 1.79769e+308; // max value of double

    // Loop over the StepPointMCs
    for (const StepPointMC& Step : Steps) {
      // Get the step time
      stepParticle = SimParticles.at(Step.trackId());
      // If the step ID (track ID) is equal to the topmost parent and its time is lower than those previously selected, update the time
      if ((stepParticle.id() == parentId) && (Step.time() < time))
          time = Step.time();
    };

    // If the time has not been updated, throw, otherwise return the updated time
    if (time > 1.7976e+308)
      throw cet::exception("LogicError") << "Time has not been updated! \n";
    return time;
  };

  void EDepFromShower::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPoints = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);

    // Validate that these data products exist
    if ((StepPoints.size() == 0) || (SimParticles.size() == 0))
      return;

    // Set up the data products for collection from this event
    EDeps.clear();
    times.clear();
    topParentIds.clear();

    // Collect the particle IDs for StepPointMCs in STMDet
    // Note - IDs are from the StepPoints as the SimParticle contains the full geneaology, this is required to keep the volume constrianed to STMDet
    SimParticleIds.clear();
    for (const StepPointMC& step : StepPoints)
      SimParticleIds.insert(SimParticles.at(step.trackId()).id());

    // Loop over all steps
    for (const StepPointMC& step : StepPoints) {
      // Select the appropriate steps
      if ((Detector == "HPGe") && (step.position().x() > xBeamCentre))
    continue;
      else if ((Detector == "LaBr") && (step.position().x() < xBeamCentre))
    continue;

      // Get the associated top particle
      stepParticle = SimParticles.at(step.trackId());
      std::tie(topParentId, pdgId) = topParent(SimParticleIds, stepParticle);
      if (std::find(SimParticleIds.begin(), SimParticleIds.end(), topParentId) == SimParticleIds.end())
    throw cet::exception("LogicError") << "The found parent ID is not a member of the data product. Exiting\n";

      // Save the associated top particle metadata
      topParentIdsIt = std::find(topParentIds.begin(), topParentIds.end(), topParentId);
      // If we've already found the parent particle
      if (topParentIdsIt != topParentIds.end()) {
    // Increment the energy by the ionizing energy deposited for this detector
    if (Detector == "HPGe")
      EDeps[topParentId] += step.ionizingEdep();
    // TODO - include the LaBr response here
      }
      // If we have not found the parent particle before
      else {
    // Add the parent particle ID to the vector of parent particle IDs
    topParentIds.emplace_back(topParentId);
    // Allocate the energy deposited as a new entry to the map of parent particle deposited energies
    if (Detector == "HPGe")
      EDeps.emplace(std::make_pair(topParentId, step.ionizingEdep()));
    // TODO - include the LaBr response here
    // Allocate the time of the top parent particle's first step as a new entry to the map of top parent particle times
    times.emplace(std::make_pair(topParentId, parentTime(event, topParentId)));
      };

      // Generate the data summary
      // Find the particle in the map of PDG IDs
      pdgIdsIt = pdgIds.find(pdgId);
      // If it is found, increment the count of these particles
      if (pdgIdsIt != pdgIds.end())
    pdgIds[pdgId] += 1;
      // If it has not been found, allocate a new entry of the particle type to the map of PDG IDs
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
    return;
  }; // end analyze

  void EDepFromShower::endJob() {
    std::cout << "==========Data summary==========" << std::endl;
    for (auto _pdgId : pdgIds)
      std::cout << "PDGID " << _pdgId.first << ": " << _pdgId.second << std::endl;
    std::cout << "================================" << std::endl;
  };

}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::EDepFromShower)
