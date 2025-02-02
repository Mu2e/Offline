// Adapted from ReadVirtualDetector_module.cc
// For StepPointMCs in virtualdetectors, generates a TTree with hit time, PDG ID, virtualdetector ID, energy, positions x, y, and z in branches "time", "virtualdetectorId", "pdgId", "E", "x", "y", and "z" respectively.
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <iostream>

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
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "TTree.h"

typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;
namespace mu2e {
  class VirtualDetectorTree : public art::EDAnalyzer {
    // Initialization
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    art::ProductToken<SimParticleCollection> SimParticlemvToken;

    // Data collection
    int pdgId = 0;
    double x = 0.0, y = 0.0, z = 0.0, mass = 0.0, KE = 0.0, time = 0.0;
    VolumeId_type virtualdetectorId = 0;

    // Data storage
    TTree* _ttree;
    std::map<int,int> pdgIds; // <id, count>

  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
      fhicl::Atom<art::InputTag> SimParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit VirtualDetectorTree(const Parameters& conf);
    void analyze(const art::Event& e);
    void endJob();
  };

  VirtualDetectorTree::VirtualDetectorTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag()))
  {
    art::ServiceHandle<art::TFileService> tfs;
    _ttree = tfs->make<TTree>( "ttree", "Virtual Detectors ttree");
    _ttree->Branch("time", &time, "time/D");
    _ttree->Branch("virtualdetectorId", &virtualdetectorId, "virtualdetectorId/l");
    _ttree->Branch("pdgId", &pdgId, "pdgId/I");
    _ttree->Branch("x", &x, "x/D");
    _ttree->Branch("y", &y, "y/D");
    _ttree->Branch("z", &z, "z/D");
    _ttree->Branch("KE", &KE, "KE/D");
  };

  void VirtualDetectorTree::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPoints = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    GlobalConstantsHandle<ParticleDataList> pdt;

    // Loop over all VD hits
    for (const StepPointMC& step : StepPoints) {
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
      // E = std::sqrt(step.momentum().mag2() + mass*mass);
      KE = std::sqrt(step.momentum().mag2() + mass*mass) - mass;
      _ttree->Fill();

      // Generate the data summary
      if (pdgIds.find(pdgId) != pdgIds.end())
        pdgIds[pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(pdgId, 1));
    };
    return;
  };

  void VirtualDetectorTree::endJob() {
    std::cout << "========= Data summary =========" << std::endl;
    std::cout << " PDG ID: count" << std::endl;
    for (auto part : pdgIds)
      std::cout << part.first << ": " << part.second << std::endl;
    std::cout << "================================" << std::endl;
  };

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VirtualDetectorTree)
