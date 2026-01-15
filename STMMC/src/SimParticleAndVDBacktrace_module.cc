// Adapted from ReadVirtualDetector_module.cc
// Generates a large ROOT TTree with SimParticle and StepPointMC information. For each SimParticle that has a StepPointMC in the specified virtual detector, gets its
//  - parentCounter - 0 for the original particle, 1 for its parent, 2 for its grandparent, etc. This is traced all the way back to the POT
//  - particleId - SimParticle ID (from GEANT)
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
// For each chosen virtual detector ID to trace back to, also gets the following StepPointMC information (0 if the SimParticle did not have a StepPointMC in that virtual detector):
//  - VD_<virtualdetectorId>_x [mm] - x position of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_y [mm] - y position of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_z [mm] - z position of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_px [MeV/c] - x momentum of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_py [MeV/c] - y momentum of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_pz [MeV/c] - z momentum of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_time [ns] - time of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_Ekin [MeV] - kinetic energy of the StepPointMC in that virtual detector
//  - VD_<virtualdetectorId>_pdgId - PDG ID of the SimParticle at that StepPointMC
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

namespace mu2e {
  class SimParticleAndVDBacktrace : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<unsigned> StartVirtualDetectorId{Name("StartVirtualDetectorId"), Comment("ID of the virtual detector to filter")};
        fhicl::Sequence<unsigned> TraceVirtualdetectorIds{Name("TraceVirtualdetectorIds"), Comment("IDs of the virtual detectors to trace back to the origin")};
        fhicl::OptionalAtom<int> consecutiveEmptyFileThreshold{Name("consecutiveEmptyFileThreshold"), Comment("Number of consecutive empty files before stopping the job")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit SimParticleAndVDBacktrace(const Parameters& conf);
      void analyze(const art::Event& e);
      void addToTree(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps);
      void addParentToTree(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps);
      void traceSteps(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps);
      void clear_virtualdetector_vectors();
      void populate_virtualdetector_vector(const StepPointMC & step_it, const unsigned index);
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      unsigned StartVirtualDetectorId = 0;
      std::vector<unsigned> TraceVirtualdetectorIds;
      int consecutiveEmptyFileThreshold = 0;

      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0, creationCode = 0, consecutiveEmptyFileCounter = 0, nTraceVirtualdetectorIds = 0;
      double x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, time = 0.0, Ekin = 0.0;
      std::vector<double> VD_x, VD_y, VD_z, VD_px, VD_py, VD_pz, VD_time, VD_Ekin;
      std::vector<int> VD_pdgId;
      uint16_t parentCounter = 0;
      cet::map_vector_key particleId;
      TTree* ttree;
      CLHEP::Hep3Vector startPosition, startMomentum;
      art::Ptr<SimParticle> parent;
      std::vector<std::vector<double>*> VD_vectors;
  };

    SimParticleAndVDBacktrace::SimParticleAndVDBacktrace(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    StartVirtualDetectorId(conf().StartVirtualDetectorId()),
    TraceVirtualdetectorIds(conf().TraceVirtualdetectorIds()) {
      consecutiveEmptyFileThreshold = conf().consecutiveEmptyFileThreshold() ? *(conf().consecutiveEmptyFileThreshold()) : 10;
      if (std::find(TraceVirtualdetectorIds.begin(), TraceVirtualdetectorIds.end(), StartVirtualDetectorId) == TraceVirtualdetectorIds.end()) {
        TraceVirtualdetectorIds.push_back(StartVirtualDetectorId);
      };
      nTraceVirtualdetectorIds = TraceVirtualdetectorIds.size();
      VD_vectors = {&VD_x, &VD_y, &VD_z, &VD_px, &VD_py, &VD_pz, &VD_time, &VD_Ekin};
      for (auto vec : VD_vectors) {
        vec->resize(nTraceVirtualdetectorIds, 0.0);
      };
      VD_pdgId.resize(nTraceVirtualdetectorIds, 0);

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
      for (size_t i=0; i<TraceVirtualdetectorIds.size(); i++) {
        ttree->Branch(Form("VD_%u_x", TraceVirtualdetectorIds[i]), &VD_x[i], Form("VD_%u_x/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_y", TraceVirtualdetectorIds[i]), &VD_y[i], Form("VD_%u_y/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_z", TraceVirtualdetectorIds[i]), &VD_z[i], Form("VD_%u_z/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_px", TraceVirtualdetectorIds[i]), &VD_px[i], Form("VD_%u_px/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_py", TraceVirtualdetectorIds[i]), &VD_py[i], Form("VD_%u_py/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_pz", TraceVirtualdetectorIds[i]), &VD_pz[i], Form("VD_%u_pz/D", TraceVirtualdetectorIds[i])); // mm
        ttree->Branch(Form("VD_%u_time", TraceVirtualdetectorIds[i]), &VD_time[i], Form("VD_%u_time/D", TraceVirtualdetectorIds[i])); // ns
        ttree->Branch(Form("VD_%u_Ekin", TraceVirtualdetectorIds[i]), &VD_Ekin[i], Form("VD_%u_Ekin/D", TraceVirtualdetectorIds[i])); // MeV
        ttree->Branch(Form("VD_%u_pdgId", TraceVirtualdetectorIds[i]), &VD_pdgId[i], Form("VD_%u_pdgId/I", TraceVirtualdetectorIds[i]));
      };
    };

  void SimParticleAndVDBacktrace::populate_virtualdetector_vector(const StepPointMC & step_it, const unsigned index) {
    const CLHEP::Hep3Vector momentum = step_it.momentum();
    const CLHEP::Hep3Vector position = step_it.position();
    const double mass = pdt->particle(step_it.simParticle()->pdgId()).mass();
    const double EKin = std::sqrt(momentum.mag2() + mass * mass) - mass;
    VD_x[index] = position.x();
    VD_y[index] = position.y();
    VD_z[index] = position.z();
    VD_px[index] = momentum.x();
    VD_py[index] = momentum.y();
    VD_pz[index] = momentum.z();
    VD_time[index] = step_it.time();
    VD_Ekin[index] = EKin;
    VD_pdgId[index] = step_it.simParticle()->pdgId();
  };

  void SimParticleAndVDBacktrace::traceSteps(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps) {
    unsigned step_virtualdetectorId = 0, index = 0;
    for (const StepPointMC & step_it : steps) {
      if (step_it.simParticle()->id() != particle->id()) continue;
      step_virtualdetectorId = step_it.virtualDetectorId();
      auto it = std::find(TraceVirtualdetectorIds.begin(), TraceVirtualdetectorIds.end(), step_virtualdetectorId);
      if (it == TraceVirtualdetectorIds.end()) continue;
      index = std::distance(TraceVirtualdetectorIds.begin(), it);
      populate_virtualdetector_vector(step_it, index);
    };
  };

  void SimParticleAndVDBacktrace::clear_virtualdetector_vectors() {
    for (auto & vec : VD_vectors) {
      std::fill(vec->begin(), vec->end(), 0.0);
    };
    std::fill(VD_pdgId.begin(), VD_pdgId.end(), 0);
    return;
  };

  void SimParticleAndVDBacktrace::addToTree(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps) {
    // Assigne all the SimParticle variables to the tree variables
    parentCounter = 0;
    particleId = particle->id();
    pdgId = particle->pdgId();
    creationCode = particle->creationCode();
    startPosition = particle->startPosition();
    x = startPosition.x();
    y = startPosition.y();
    z = startPosition.z();
    startMomentum = particle->startMomentum();
    px = startMomentum.x();
    py = startMomentum.y();
    pz = startMomentum.z();
    time = particle->startGlobalTime();
    Ekin = particle->endKineticEnergy();

    // Trace the StepPointMCs in the relevant virtualdetectors
    clear_virtualdetector_vectors();
    traceSteps(particle, steps);
    ttree->Fill();

    // Now trace the parents
    if (!particle->hasParent()) return;

    // Require separate definitions of addParentToTree because of data accessors
    parent = particle->parent();
    addParentToTree(parent, steps);

    while (parent->hasParent()) {
      parent = parent->parent();
      addParentToTree(parent, steps);
    };
  };

  void SimParticleAndVDBacktrace::addParentToTree(const art::Ptr<SimParticle> &particle, const std::vector<StepPointMC> steps) {
    // Increment the parent counter
    parentCounter++;

    // Assign all the SimParticle variables to the tree variables
    particleId = particle->id();
    pdgId = particle->pdgId();
    creationCode = particle->creationCode();
    startPosition = particle->startPosition();
    x = startPosition.x();
    y = startPosition.y();
    z = startPosition.z();
    startMomentum = particle->startMomentum();
    px = startMomentum.x();
    py = startMomentum.y();
    pz = startMomentum.z();
    time = particle->startGlobalTime();
    Ekin = particle->endKineticEnergy();

    // Trace the StepPointMCs in the relevant virtualdetectors
    clear_virtualdetector_vectors();
    traceSteps(particle, steps);
    ttree->Fill();
  };

  void SimParticleAndVDBacktrace::analyze(const art::Event& event) {
    // Get the data products
    auto stepHandle = event.getHandle< std::vector<StepPointMC> >(StepPointMCsToken);
    if (!stepHandle || stepHandle->empty()) {
      consecutiveEmptyFileCounter++;
      if (consecutiveEmptyFileCounter > consecutiveEmptyFileThreshold) {
        throw cet::exception("LogicError", "Too many consecutive empty files, are you sure you have the correct data product name?\n");
      };
      return;
    };
    const std::vector<StepPointMC> StepPointMCs = *stepHandle;
    consecutiveEmptyFileCounter = 0;

    for (const StepPointMC& step : StepPointMCs) {
      // Filter by virtual detector ID if requested
      if (step.virtualDetectorId() != StartVirtualDetectorId) continue;

      // Get the associated particle
      const art::Ptr<SimParticle> particle = step.simParticle();
      addToTree(particle, StepPointMCs);
    };
    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleAndVDBacktrace)
