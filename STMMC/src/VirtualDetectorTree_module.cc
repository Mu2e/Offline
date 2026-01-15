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
#include "fhiclcpp/types/OptionalAtom.h"

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
        fhicl::OptionalAtom<int> consecutiveEmptyFileThreshold{Name("consecutiveEmptyFileThreshold"), Comment("Number of consecutive empty files before stopping the job")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VirtualDetectorTree(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0, consecutiveEmptyFileCounter = 0, consecutiveEmptyFileThreshold = 0, simParticleId=0;
      double x = 0.0, y = 0.0, z = 0.0, mass = 0.0, Ekin = 0.0, Etot = 0.0, time = 0.0, p = 0.0, p2 = 0.0, px = 0.0, py = 0.0, pz = 0.0;
      VolumeId_type virtualdetectorId = 0;
      TTree* ttree;
      std::map<int, int> pdgIds; // <id, count>
      uint simParticleIdKey = 0;
  };

  VirtualDetectorTree::VirtualDetectorTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())) {
      consecutiveEmptyFileThreshold = conf().consecutiveEmptyFileThreshold() ? *(conf().consecutiveEmptyFileThreshold()) : 10;
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "Virtual Detectors ttree");
      ttree->Branch("time", &time, "time/D"); // ns
      ttree->Branch("virtualdetectorId", &virtualdetectorId, "virtualdetectorId/l");
      ttree->Branch("pdgId", &pdgId, "pdgId/I");
      ttree->Branch("x", &x, "x/D"); // mm
      ttree->Branch("y", &y, "y/D"); // mm
      ttree->Branch("z", &z, "z/D"); // mm
      ttree->Branch("px", &px, "px/D"); // mm
      ttree->Branch("py", &py, "py/D"); // mm
      ttree->Branch("pz", &pz, "pz/D"); // mm
      ttree->Branch("p", &p, "p/D"); // mm
      ttree->Branch("p2", &p2, "p2/D"); // mm
      ttree->Branch("mass", &mass, "mass/D"); // MeV/c^2
      ttree->Branch("Ekin", &Ekin, "Ekin/D"); // MeV
      ttree->Branch("Etot", &Etot, "Etot/D"); // MeV
      ttree->Branch("SimParticleId", &simParticleId, "SimParticleId/i");
    };

  void VirtualDetectorTree::analyze(const art::Event& event) {
    auto stepHandle = event.getHandle< std::vector<StepPointMC> >(StepPointMCsToken);
    if (!stepHandle || stepHandle->empty()) {
      consecutiveEmptyFileCounter++;
      return;
    }

    // Try to get SimParticles (assuming it's a map or vector)
    auto simHandle = event.getHandle< SimParticleCollection >(SimParticlemvToken);
    if (!simHandle || simHandle->empty()) {
      consecutiveEmptyFileCounter++;
      return;
    }
    if (consecutiveEmptyFileCounter > consecutiveEmptyFileThreshold) {
      throw cet::exception("LogicError", "Too many consecutive empty files, stopping the job");
    }
    auto const& StepPointMCs = *stepHandle;
    auto const& SimParticles = *simHandle;
    consecutiveEmptyFileCounter = 0;
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
      px = step.momentum().x();
      py = step.momentum().y();
      pz = step.momentum().z();
      p2 = step.momentum().mag2();
      p = std::sqrt(p2);
      mass = pdt->particle(pdgId).mass();
      Etot = std::sqrt(p2 + mass * mass); // Total energy
      Ekin = Etot - mass; // Subtract the rest mass
      simParticleId = particle.id().asInt();
      if (Ekin < 0)
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
