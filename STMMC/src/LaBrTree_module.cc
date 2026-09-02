// Adapted from HPGeTree_module.cc
// For StepPointMCs in STMDet, generates a TTree with energy in branch "E" and time in branch "time", for the LaBr detector
//  - Iterate over the StepPointMCs, determine the associated SimParticle
//  - Find the most parent particle of the associated SimParticle within the associated STMDet
//  - Increment the energy associated with the parent particle
//    - If the parent particle ID doesn't exist in the collection, add a new entry to the parent particle ID vector, energy deposited map, and time map
//  - Determine the PDG ID of the top particles, increment the mapped counter
//    - If the PDG ID entry does not exist in the map, create a new entry.
//  - Print the PDG ID counts at the end of the job
//
// WARNING - this module reports the RAW IONIZING ENERGY DEPOSITED in the LaBr crystal, summed
// per top parent particle in exactly the same way as HPGeTree does for the HPGe crystal. It is
// NOT a LaBr detector response. In particular the following are NOT applied:
//  - the scintillation light yield and its non-proportionality with energy, which for LaBr3:Ce
//    is a significant, energy dependent effect
//  - the photosensor quantum efficiency and gain
//  - the intrinsic energy resolution, which is far worse than that of a HPGe crystal
// The deposited energy is still physically meaningful - the scintillation light produced is
// driven by it - so this tree is useful for relative studies and for comparing the flux and
// spectrum reaching the two crystals. Do not read the "E" branch as a measured LaBr spectrum.
//
// FIXME - replace the raw ionizing energy sum below with a proper LaBr response once the light
// yield, non-proportionality and resolution model is defined in Offline. At that point this
// module and HPGeTree_module.cc should probably be merged back into one detector-aware module,
// which is what the "Detector" parameter of HPGeTree was originally intended for - it currently
// throws for "LaBr" because that response is undefined, and this module exists to fill that gap
// without changing HPGeTree's behaviour.
//
// Input Parameters
//  - StepPointMCsTag - tag of data product containing the StepPoints for STMDet
//  - SimParticlemvTag - tag of data product containing the SimParticles for STMDet
// Original author: Ivan Logashenko, Pawel Plesniak
// Adapted by: Yongyi Wu, Aug. 2026

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
  class LaBrTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> simParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit LaBrTree(const Parameters& conf);
      std::tuple<key_type, int> topParent(std::set<key_type>& SimParticleIDs, const SimParticle particle);
      double parentTime(const art::Event& event, key_type parentId);
      void analyze(const art::Event& event);
      void endJob();
    private:
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
      unsigned long nNoParent = 0; // steps whose SimParticle had no reachable parent

      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      art::Ptr<SimParticle> parent;
  };

  LaBrTree::LaBrTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().simParticlemvTag())) {
      // Make it obvious in the job log that this is not a detector response
      mf::LogWarning("LaBrTree")
        << "The 'E' branch of this tree is the raw ionizing energy deposited in the LaBr crystal.\n"
        << "No scintillation light yield, non-proportionality, photosensor response or energy\n"
        << "resolution is applied, so this is NOT a simulated LaBr spectrum. See the FIXME at the\n"
        << "top of LaBrTree_module.cc.\n";

      // Set up TTree
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "Detector ttree");
      ttree->Branch("E", &E, "E/D");
      ttree->Branch("time", &time, "time/D");
  };

  std::tuple<key_type, int> LaBrTree::topParent(std::set<key_type>& SimParticleIds, const SimParticle particle) {
    // Get the particle parent
    parent = particle.parent();

    // A primary has no parent, and a compressed SimParticleCollection can leave a non-null Ptr
    // whose target was dropped. Either way the walk stops here and the particle is its own top
    // parent - the same answer as the STMDet test below, which cannot be reached without a
    // dereferenceable Ptr. Without this, art throws ProductNotFound (ProductID 0).
    if (parent.isNull() || !parent.isAvailable()) {
      ++nNoParent;
      return std::make_tuple(particle.id(), particle.pdgId());
    }

    // If the passed particle has no parent in STMDet, return its ID and PDG ID
    if (std::find(SimParticleIds.begin(), SimParticleIds.end(), parent->id()) == SimParticleIds.end())
      return std::make_tuple(particle.id(), particle.pdgId());

    // If the particle has a parent in STMDet, update the particle parent. The genealogy can run
    // out before the STMDet boundary does, so the Ptr is checked before each dereference.
    while (parent->parent().isNonnull() && parent->parent().isAvailable() &&
           std::find(SimParticleIds.begin(), SimParticleIds.end(), parent->parent()->id()) != SimParticleIds.end())
      parent = parent->parent();

    // Return the ID and PDG ID
    return std::make_tuple(parent->id(), parent->pdgId());
  };

  double LaBrTree::parentTime(const art::Event& event, key_type parentId) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken); // TODO - resture this so we don't access the SimParticles from the data product directly, but through the parent particle

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

  void LaBrTree::analyze(const art::Event& event) {
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
      // Select the steps in the LaBr crystal, which sits on the +x side of the beam centre
      if (step.position().x() < xBeamCentre)
        continue;

      // Get the associated top particle
      stepParticle = SimParticles.at(step.trackId());
      std::tie(topParentId, pdgId) = topParent(SimParticleIds, stepParticle);
      // Precautionary check if the SimParticle ID is in STMDet
      if (std::find(SimParticleIds.begin(), SimParticleIds.end(), topParentId) == SimParticleIds.end())
        throw cet::exception("LogicError") << "The found parent ID is not a member of the data product\n";

      // Collate the data
      // FIXME - this is the raw ionizing energy deposited, not the LaBr response. See the top of this file.
      topParentIdsIt = std::find(topParentIds.begin(), topParentIds.end(), topParentId);
      if (topParentIdsIt != topParentIds.end()) {
        EDeps[topParentId] += step.ionizingEdep();
      }
      else {
        topParentIds.emplace_back(topParentId);
        EDeps.emplace(std::make_pair(topParentId, step.ionizingEdep()));
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

  void LaBrTree::endJob() {
    mf::LogInfo log("Detector tree");
    log << "==========Data summary==========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
    // Not an error - primaries legitimately have no parent - but a large count means the
    // genealogy is truncated, so the "top parent" is only the top of what was kept.
    if (nNoParent > 0)
      log << "Steps whose SimParticle had no reachable parent (treated as their own top parent): "
          << nNoParent << "\n";
    mf::LogWarning("LaBrTree")
      << "Reminder - the 'E' branch is raw ionizing energy deposited, not a LaBr detector response.\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrTree)
