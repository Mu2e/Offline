// stdlib includes
#include <limits>
#include <map>
#include <cmath>
#include <utility>

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
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
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
  class GeometricShiftTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> BeforeShiftTag{ Name("BeforeShiftTag"), Comment("Tag of the StepPointMCs before running the shift")};
        fhicl::Atom<art::InputTag> AfterShiftTag{ Name("AfterShiftTag"), Comment("Tag of the StepPointMCs after running the shift")};
        fhicl::Atom<long unsigned int> VDId{ Name("VDId"), Comment("Virtual Detector ID to filter on"), 101}; // Default to VD101
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit GeometricShiftTree(const Parameters& conf);
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      TTree* ttree = nullptr;

      // Kinematic variables
      double Ekin_before = 0.0, Ekin_after = 0.0;
      double px_before = 0.0, px_after = 0.0, py_before = 0.0, py_after = 0.0, pz_before = 0.0, pz_after = 0.0;

      // Spatial variables to plot the shift quantities
      double x_before = 0.0, x_after = 0.0, y_before = 0.0, y_after = 0.0, z_before = 0.0, z_after = 0.0;

      int processedSteps = 0;
      long unsigned int vdId_ = 101;
      art::ProductToken<StepPointMCCollection> BeforeShiftToken;
      art::ProductToken<StepPointMCCollection> AfterShiftToken;
  };

  GeometricShiftTree::GeometricShiftTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    vdId_(conf().VDId()),
    BeforeShiftToken(consumes<StepPointMCCollection>(conf().BeforeShiftTag())),
    AfterShiftToken(consumes<StepPointMCCollection>(conf().AfterShiftTag())) {
      // Set up TTree
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>( "ttree", "GeometricShiftTree ttree");

      // Energy and Momentum Branches
      ttree->Branch("Ekin_before",  &Ekin_before, "Ekin_before/D");
      ttree->Branch("Ekin_after",   &Ekin_after,  "Ekin_after/D");
      ttree->Branch("px_before",    &px_before,   "px_before/D");
      ttree->Branch("px_after",     &px_after,    "px_after/D");
      ttree->Branch("py_before",    &py_before,   "py_before/D");
      ttree->Branch("py_after",     &py_after,    "py_after/D");
      ttree->Branch("pz_before",    &pz_before,   "pz_before/D");
      ttree->Branch("pz_after",     &pz_after,    "pz_after/D");

      // Position Branches
      ttree->Branch("x_before",     &x_before,    "x_before/D");
      ttree->Branch("x_after",      &x_after,     "x_after/D");
      ttree->Branch("y_before",     &y_before,    "y_before/D");
      ttree->Branch("y_after",      &y_after,     "y_after/D");
      ttree->Branch("z_before",     &z_before,    "z_before/D");
      ttree->Branch("z_after",      &z_after,     "z_after/D");
  }

  void GeometricShiftTree::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& stepsBefore = event.getProduct(BeforeShiftToken);
    auto const& stepsAfter = event.getProduct(AfterShiftToken);

    // Validate that these data products exist
    if (stepsBefore.empty() || stepsAfter.empty()) {
      return;
    }

    // Access the Particle Data List for mass lookups to calculate kinetic energy
    GlobalConstantsHandle<ParticleDataList> pdt;

    // Collect the data products.
    // Map key is a pair of {SimParticleKey, Time} to guarantee 100% bulletproof matching
    std::map<std::pair<art::Ptr<SimParticle>::key_type, double>, const StepPointMC*> beforeMap;

    for (const auto& step : stepsBefore) {
        // STRICT FILTER: Only map hits from the target Virtual Detector BEFORE the shift
        if (step.volumeId() == vdId_) {
            beforeMap[{step.simParticle().key(), step.time()}] = &step;
        }
    }

    // Collect the data to the TTree
    for (const auto& stepAfter : stepsAfter) {
        // Match based on SimParticle ID and Time
        auto key = std::make_pair(stepAfter.simParticle().key(), stepAfter.time());
        auto it = beforeMap.find(key);

        // Make sure the exact step entry was present in the filtered beforeMap
        // (This inherently validates that the step originally came from vdId_)
        if (it != beforeMap.end()) {
            const auto* stepBefore = it->second;

            // Lookup particle mass dynamically based on its PDG ID
            double mass = pdt->particle(stepBefore->simParticle()->pdgId()).mass();

            // Fill Spatial variables for before and after shift
            x_before = stepBefore->position().x();
            y_before = stepBefore->position().y();
            z_before = stepBefore->position().z();

            x_after = stepAfter.position().x();
            y_after = stepAfter.position().y();
            z_after = stepAfter.position().z();

            // Fill Kinematic variables for before shift
            px_before = stepBefore->momentum().x();
            py_before = stepBefore->momentum().y();
            pz_before = stepBefore->momentum().z();
            Ekin_before = std::sqrt(stepBefore->momentum().mag2() + mass * mass) - mass;

            // Fill Kinematic variables for after shift
            px_after = stepAfter.momentum().x();
            py_after = stepAfter.momentum().y();
            pz_after = stepAfter.momentum().z();
            Ekin_after = std::sqrt(stepAfter.momentum().mag2() + mass * mass) - mass;

            ttree->Fill();
            processedSteps++;
        }
    }

    return;
  } // end analyze

  void GeometricShiftTree::endJob() {
    mf::LogInfo log("Detector tree");
    log << "==========Data summary==========\n";
    log << "Processed " << processedSteps << " StepPointMCs validated at VD" << vdId_ << " before shift\n";
    log << "================================\n";
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GeometricShiftTree)
