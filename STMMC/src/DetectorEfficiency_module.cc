// Sums the energy associated with a simulated POT in VD101
// Original author: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TNtuple.h"


namespace mu2e {
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
            const unsigned long virtualdetectorId = 101;
            double E = 0.0;
            unsigned long acceptedSteps = 0, rejectedSteps = 0, countedEvents = 0;
    };

    DetectorEfficiency::DetectorEfficiency(const Parameters& config):
        art::EDAnalyzer{config} {
            art::ServiceHandle<art::TFileService> tfs;
            _nt = tfs->make<TNtuple>("nt", "Energy Deposit", "E");
        };

    void DetectorEfficiency::analyze(art::Event const& event) {
        countedEvents++;
        art::Handle<StepPointMCCollection> StepPointMCs;
        event.getByLabel("g4run:virtualdetector", StepPointMCs);
        if (StepPointMCs->size() == 0)
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
        mf::LogInfo log("Detector Efficiency");
        log << "========= Data summary =========\n";
        log << "Accepted steps: " << acceptedSteps << "\n";
        log << "Rejected steps: " << rejectedSteps << "\n";
        log << "Counted events: " << countedEvents << "\n";
        log << "================================\n";
    };
};

DEFINE_ART_MODULE(mu2e::DetectorEfficiency)
