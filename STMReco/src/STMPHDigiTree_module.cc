// NTuple dumper for Detector pulse heights

// stdlib includes
#include <limits>
#include <map>

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

// Offline includes I added
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TH1D.h"


// Mu2e type definitions

namespace mu2e {
    class STMPHDigiTree : public art::EDAnalyzer {
        public:
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          struct Config {
            fhicl::Atom<art::InputTag> phDigiTag{ Name("phDigiTag"), Comment("Input Tag for STMPHDigiCollectionMap")};
            // Already seperates by detector name -> makeSTMDigis:phHPGe or phLaBr
          };
          using Parameters = art::EDAnalyzer::Table<Config>;
          explicit STMPHDigiTree(const Parameters& conf);

        private:
          void analyze(const art::Event& event) override;
          void beginJob() override;
          //void endJob() override;

          // New token for map
          art::ProductToken<STMPHDigiCollectionMap> _stmPHDigisMapToken;
          STMChannel _channel;

          // Store STM PH Digi information
          int16_t pulseHeight = 0;
          uint32_t uncalibratedTime = 0;

          // Store file information
          Int_t art_evt;
          Int_t run;
          Int_t subrun;

          // Store from EventHeader
          uint64_t ewt = 0;
          uint8_t evtMode = 0;
          uint64_t adcClock = 0;
          uint64_t dtcClock = 0;

          // Tree reference
          TTree* ttree = nullptr;

          // STM PH spectrum
          TH1D* _phSpectrum = nullptr;
    };

    STMPHDigiTree::STMPHDigiTree(const Parameters& config) :
        art::EDAnalyzer{config},
        _stmPHDigisMapToken(consumes<STMPHDigiCollectionMap>(config().phDigiTag())),
        _channel(STMUtils::getChannel(config().phDigiTag()))
        {}

    void STMPHDigiTree::beginJob(){
        // Set up TTree here
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "Detector ttree");
        ttree->Branch("pulse_height", &pulseHeight, "pulse_height/S");
        ttree->Branch("uncalibrated_time", &uncalibratedTime, "uncalibrated_time/i");

        // Event Information
        ttree->Branch("art_event", &art_evt, "art_event/I");
        ttree->Branch("run", &run, "run/I");
        ttree->Branch("subrun", &subrun, "subrun/I");

        // EventHeader information
        ttree->Branch("EWT", &ewt, "EWT/l");
        ttree->Branch("event_mode", &evtMode, "event_mode/b");
        ttree->Branch("adc_clock", &adcClock, "adc_clock/l");
        ttree->Branch("dtc_clock", &dtcClock, "dtc_clock/l");

        // Set up histogram
        std::string phSpectrumTitle = "PH Spectrum (" + _channel.name() + ")" ; // Builds title PH Spectrum (HPGe/ LaBr)

        _phSpectrum = tfs->make<TH1D>("phSpectrum",
            (phSpectrumTitle +";PulseHeight;Count").c_str(),
            10000, -10000, 0);
    };

    void STMPHDigiTree::analyze(const art::Event& event) {
        // We fill art based information here
        art_evt = event.event();
        run = event.run();
        subrun = event.subRun();

        // Get handle for PH Digis
        auto phDigisMapHandle = event.getValidHandle(_stmPHDigisMapToken);
        const auto& phDigiMap = *phDigisMapHandle;
        for (const auto& i_phDigiMap : phDigiMap){
            // get EventHeader
            const auto& header = i_phDigiMap.first;
            // store related info
            ewt = header.eventWindowTag();
            evtMode = header.eventMode();
            adcClock = header.adcClock();
            dtcClock = header.dtcClock();

            // get PH Digi Collection
            const auto& phDigis = i_phDigiMap.second;

            // Second loop for PH Digis
            for (const auto& phDigi : phDigis){
                pulseHeight = phDigi.energy();
                uncalibratedTime = phDigi.time();
                // fill tree
                ttree->Fill();
                // fill histogram
                _phSpectrum->Fill(pulseHeight);
            }
        }
    } // end of analyze

}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::STMPHDigiTree)
