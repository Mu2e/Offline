// Sums the energy associated with a simulated POT in VD01
// Original author: Pawel Plesniak

// stdlib
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>

// art
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl
#include "canvas/Utilities/InputTag.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT
#include "art_root_io/TFileService.h"
#include "TNtuple.h"


namespace mu2e{
    class DetectorEfficiency : public art::EDAnalyzer {
      public:
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          struct Config {};
          using Parameters=art::EDAnalyzer::Table<Config>;

          explicit DetectorEfficiency(const Parameters& pset);
          void analyze(art::Event const& event);
          void endJob();
      private:
          bool verbose = false;
          TNtuple* _nt = nullptr;
          double E = 0.0;
          int acceptedSteps = 0, rejectedSteps = 0, countedEvents = 0;
          const unsigned long virtualdetectorId = 101;
    };

    DetectorEfficiency::DetectorEfficiency(const Parameters& config):
        art::EDAnalyzer{config} {
            art::ServiceHandle<art::TFileService> tfs;
            _nt = tfs->make<TNtuple>("nt", "Energy Deposit", "E");
            std::cout << "Initializing DetectorEfficiency TTree" << std::endl;
        };

    void DetectorEfficiency::analyze(art::Event const& event) {
        countedEvents++;
        art::Handle<StepPointMCCollection> StepPointMCs;
        event.getByLabel("g4run:virtualdetector", StepPointMCs);
        if (StepPointMCs.size() == 0)
            throw cet::exception("DataError") << "Requested data not found";

        // Collect the associated StepPointMC energy
        E = 0.0;
        for (const StepPointMC &step : *StepPointMCs) {
            if (step.virtualDetectorId() == virtualdetectorId) {
                E += step.ionizingEdep();
            acceptedSteps++;
        }
        else
            rejectedSteps++;
        };
        _nt->Fill(E);
        return;
    };

    void DetectorEfficiency::endJob() {
        mf::LogInfo log("Unfiltered virtual detector hits summary");
        log << "Accepted steps: " << acceptedSteps << "\n";
        log << "Rejected steps: " << rejectedSteps << "\n";
        log << "Counted events: " << countedEvents << "\n";
    };
};

DEFINE_ART_MODULE(mu2e::DetectorEfficiency)
