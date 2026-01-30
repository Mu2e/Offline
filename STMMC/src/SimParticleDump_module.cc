// Adapted from ReadVirtualDetector_module.cc
// For StepPointMCs in virtualdetectors, generates a TTree showing the provenance of the generated SimParticles with the following branches
//  - pdgId - PDG ID
//  - creationCode - creation code. See Offline/MCDataProducts/inc/ProcessCode.hh for details
//  - x [mm] - start x position of the SimParticle
//  - y [mm] - start y position of the SimParticle
//  - z [mm] - start z position of the SimParticle
//  - px [MeV/c] - start x momentum of the SimParticle
//  - py [MeV/c] - start y momentum of the SimParticle
//  - pz [MeV/c] - start z momentum of the SimParticle
//  - time [ns] - start global time of the SimParticle
//  - Ekin [MeV] - end kinetic energy of the SimParticle
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// TODO before upload - validate that the code works and generates the same particle trace as the SimParticleAndVDBacktrace module

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

namespace mu2e {
  class SimParticleDump : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
        fhicl::Atom<unsigned> FilterVirtualDetectorId{Name("FilterVirtualDetectorId"), Comment("ID of the virtual detector to filter")};
        fhicl::OptionalAtom<int> consecutiveEmptyFileThreshold{Name("consecutiveEmptyFileThreshold"), Comment("Number of consecutive empty files before stopping the job")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit SimParticleDump(const Parameters& conf);
      void analyze(const art::Event& e);
      void addToTree(const SimParticle& particle, TTree* ttree);
      void addParentToTree(art::Ptr<SimParticle> &_particle, TTree* ttree);
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0, creationCode = 0, consecutiveEmptyFileCounter = 0, consecutiveEmptyFileThreshold = 0;
      double x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, time = 0.0, Ekin = 0.0;
      unsigned virtualdetectorId = 0, FilterVirtualDetectorId = 0;
      bool forward = false;
      uint16_t parentCounter = 0;
      cet::map_vector_key particleId;
      TTree* ttree;
      CLHEP::Hep3Vector startPosition, startMomentum;
      art::Ptr<SimParticle> parent;
  };

  SimParticleDump::SimParticleDump(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    FilterVirtualDetectorId(conf().FilterVirtualDetectorId()) {
      consecutiveEmptyFileThreshold = conf().consecutiveEmptyFileThreshold() ? *(conf().consecutiveEmptyFileThreshold()) : 10;
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "SimParticle ttree");
      ttree->Branch("parentCounter", &parentCounter, "parentCounter/s");
      ttree->Branch("particleId", &particleId, "particleId/l");
      ttree->Branch("pdgId", &pdgId, "pdgId/I");
      ttree->Branch("creationCode", &creationCode, "creationCode/I");
      ttree->Branch("x", &x, "x/D"); // mm
      ttree->Branch("y", &y, "y/D"); // mm
      ttree->Branch("z", &z, "z/D"); // mm
      ttree->Branch("px", &px, "px/D"); // mm
      ttree->Branch("py", &py, "py/D"); // mm
      ttree->Branch("pz", &pz, "pz/D"); // mm
      ttree->Branch("time", &time, "time/D"); // ns
      ttree->Branch("Ekin", &Ekin, "Ekin/D"); // MeV
    };

  void SimParticleDump::addToTree(const SimParticle& particle, TTree* ttree) {
    parentCounter = 0;
    particleId = particle.id();
    pdgId = particle.pdgId();
    creationCode = particle.creationCode();
    startPosition = particle.startPosition();
    x = startPosition.x();
    y = startPosition.y();
    z = startPosition.z();
    startMomentum = particle.startMomentum();
    px = startMomentum.x();
    py = startMomentum.y();
    pz = startMomentum.z();
    time = particle.startGlobalTime();
    Ekin = particle.endKineticEnergy();
    ttree->Fill();

    if (!particle.hasParent()) {
      return;
    };
    parentCounter++;
    parent = particle.parent();
    addParentToTree(parent, ttree);
    while (parent->hasParent()) {
      parentCounter++;
      parent = parent->parent();
      addParentToTree(parent, ttree);
    };
  };

  void SimParticleDump::addParentToTree(art::Ptr<SimParticle> &_particle, TTree* ttree) {
    particleId = _particle->id();
    pdgId = _particle->pdgId();
    creationCode = _particle->creationCode();
    startPosition = _particle->startPosition();
    x = startPosition.x();
    y = startPosition.y();
    z = startPosition.z();
    startMomentum = _particle->startMomentum();
    px = startMomentum.x();
    py = startMomentum.y();
    pz = startMomentum.z();
    time = _particle->startGlobalTime();
    Ekin = _particle->endKineticEnergy();
    ttree->Fill();
  };

  void SimParticleDump::analyze(const art::Event& event) {
    auto stepHandle = event.getHandle< std::vector<StepPointMC> >(StepPointMCsToken);
    if (!stepHandle || stepHandle->empty()) {
      consecutiveEmptyFileCounter++;
      return;
    };
    auto simHandle = event.getHandle< SimParticleCollection >(SimParticlemvToken);
    if (!simHandle || simHandle->empty()) {
      consecutiveEmptyFileCounter++;
      return;
    };
    if (consecutiveEmptyFileCounter > consecutiveEmptyFileThreshold) {
      throw cet::exception("LogicError", "Too many consecutive empty files, stopping the job");
    };

    auto const& StepPointMCs = *stepHandle;
    auto const& SimParticles = *simHandle;
    consecutiveEmptyFileCounter = 0;

    for (const StepPointMC& step : StepPointMCs) {
      // Filter by virtual detector ID if requested
      if (step.virtualDetectorId() != FilterVirtualDetectorId) continue;

      // Get the associated particle
      const SimParticle& particle = SimParticles.at(step.trackId());
      addToTree(particle, ttree);
    };
    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleDump)
