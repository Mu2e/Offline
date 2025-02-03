// Adapted from ReadVirtualDetector_module.cc
// For StepPointMCs in virtualdetectors, generates a TTree with hit time, PDG ID, virtualdetector ID, kinetic energy, positions x, y, and z in branches "time", "virtualdetectorId", "pdgId", "E", "x", "y", and "z" respectively.
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <iostream>

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

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;

namespace mu2e {
  class VirtualDetectorTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VirtualDetectorTree(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0;
      double x = 0.0, y = 0.0, z = 0.0, mass = 0.0, E = 0.0, time = 0.0;
      VolumeId_type virtualdetectorId = 0;
      TTree* ttree;
      std::map<int, int> pdgIds; // <id, count>
  };

  VirtualDetectorTree::VirtualDetectorTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())) {
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "Virtual Detectors ttree");
      ttree->Branch("time", &time, "time/D"); // ns
      ttree->Branch("virtualdetectorId", &virtualdetectorId, "virtualdetectorId/l");
      ttree->Branch("pdgId", &pdgId, "pdgId/I");
      ttree->Branch("x", &x, "x/D"); // mm
      ttree->Branch("y", &y, "y/D"); // mm
      ttree->Branch("z", &z, "z/D"); // mm
      ttree->Branch("E", &E, "E/D"); // MeV
    };

  void VirtualDetectorTree::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.empty())
      throw cet::exception("DataError", "Requested data product not found");
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    if (SimParticles.empty())
      throw cet::exception("DataError", "Requested data product not found");

    // Loop over all VD hits
    for (const StepPointMC& step : StepPointMCs) {
      // Get the associated particle
      const SimParticle& particle = SimParticles.at(step.trackId());

      // Extract the parameters
      time = step.time();
      virtualdetectorId = step.virtualDetectorId();
      pdgId = particle.pdgId();
      x = step.position().x();
      y = step.position().y();
      z = step.position().z();
      mass = pdt->particle(pdgId).mass();
      E = step.momentum().mag2()/(2 * mass);
      if (E < 0)
        throw cet::exception("LogicError", "Energy is negative");
      ttree->Fill();

      // Generate the data summary
      if (pdgIds.find(pdgId) != pdgIds.end())
        pdgIds[pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(pdgId, 1));
    };
    return;
  };

  void VirtualDetectorTree::endJob() {
    mf::LogInfo log("Virtual detector tree summary");
    log << "========= Data summary =========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VirtualDetectorTree)
