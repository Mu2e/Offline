// Adapted from ReadVirtualDetector_module.cc
// Extracts the end position (stopping vertex) of particles from a SimParticleCollection
// Generates a TTree with stop time, PDG ID, and positions x, y, and z.
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <map>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/SimParticle.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class TSCTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> SimParticleTag{
            Name("SimParticleTag"),
            Comment("Tag identifying the SimParticle collection (e.g. TargetStopFilter)")
        };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit TSCTree(const Parameters& conf);
      void analyze(const art::Event& e) override;
      void endJob() override;

    private:
      art::ProductToken<SimParticleCollection> simParticleToken_;
      int pdgId = 0;
      double x = 0.0, y = 0.0, z = 0.0, time = 0.0;
      TTree* ttree;
      std::map<int, int> pdgIds; // <id, count>
  };

  TSCTree::TSCTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    simParticleToken_(consumes<SimParticleCollection>(conf().SimParticleTag())) {
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>("ttree", "Target Stops ttree");
      ttree->Branch("time", &time, "time/D"); // ns
      ttree->Branch("pdgId", &pdgId, "pdgId/I");
      ttree->Branch("x", &x, "x/D"); // mm
      ttree->Branch("y", &y, "y/D"); // mm
      ttree->Branch("z", &z, "z/D"); // mm
  }

  void TSCTree::analyze(const art::Event& event) {
    auto simHandle = event.getHandle<SimParticleCollection>(simParticleToken_);
    if (!simHandle || simHandle->empty()) {
      return;
    }

    auto const& SimParticles = *simHandle;

    // SimParticleCollection is a cet::map_vector, iterate over the pairs
    for (auto const& iter : SimParticles) {
      const SimParticle& particle = iter.second;

      // Extract initial filtering parameters
      pdgId = particle.pdgId();
      z = particle.endPosition().z();

      // Filter for Muons (13) stopping within the ST (5470 mm to 6271 mm)
      if (pdgId == 13 && z >= 5470.0 && z <= 6271.0) {
          x = particle.endPosition().x();
          y = particle.endPosition().y();
          time = particle.endGlobalTime();

          ttree->Fill();
          pdgIds[pdgId]++;
      }
    }
  }

  void TSCTree::endJob() {
    mf::LogInfo log("Target Stops tree summary");
    log << "========= Data summary =========\n";
    for (auto const& part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "================================\n";
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TSCTree)
