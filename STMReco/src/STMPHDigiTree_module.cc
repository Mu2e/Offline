// NTuple dumper for Detector pulse heights

// stdlib includes
#include <limits>

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
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

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
            fhicl::Atom<art::InputTag> phDigiTag{ Name("phDigiTag"), Comment("Input Tag for STMPHDigiCollection")};
            // Already seperates by detector name -> makeSTMDigis:phHPGe or phLaBr
          };
          using Parameters = art::EDAnalyzer::Table<Config>;
          explicit STMPHDigiTree(const Parameters& conf);

        private:
          void analyze(const art::Event& event) override;
          void beginJob() override;
          //void endJob() override;

          // Token for STM PH Digis
          art::ProductToken<STMPHDigiCollection> _stmPHDigisToken;
          STMChannel _channel;

          // Store STM PH Digi information
          int16_t pulseHeight_ = 0;
          uint32_t uncalibratedTime_ = 0;

          // Store file information
          Int_t evt_;
          Int_t run_;
          Int_t subrun_;

          // Tree reference
          TTree* ttree = nullptr;

          // STM PH spectrum
          TH1D* _phSpectrum = nullptr;
    };

    STMPHDigiTree::STMPHDigiTree(const Parameters& config) :
        art::EDAnalyzer{config},
        _stmPHDigisToken(consumes<STMPHDigiCollection>(config().phDigiTag())),
        _channel(STMUtils::getChannel(config().phDigiTag()))
        {}

    void STMPHDigiTree::beginJob(){
        // Set up TTree here
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "Detector ttree");
        ttree->Branch("pulse_height", &pulseHeight_, "pulse_height/S");
        ttree->Branch("uncalibrated_time", &uncalibratedTime_, "uncalibrated_time/i");

        // Event Information
        ttree->Branch("event", &evt_, "event/I");
        ttree->Branch("run", &run_, "run/I");
        ttree->Branch("subrun", &subrun_, "subrun/I");

        // Set up histogram
        std::string phSpectrumTitle = "PH Spectrum (" + _channel.name() + ")" ; // Builds title PH Spectrum (HPGe/ LaBr)

        _phSpectrum = tfs->make<TH1D>("phSpectrum",
            (phSpectrumTitle +";PulseHeight;Count").c_str(),
            10000, -10000, 0);
    };

    void STMPHDigiTree::analyze(const art::Event& event) {

        // We fill here information
        // you can reference event and get ifno
        evt_ = event.event();
        run_ = event.run();
        subrun_=event.subRun();

        // Get handle for PH Digis
        auto phDigisHandle = event.getValidHandle(_stmPHDigisToken);
        const auto& phDigis = *phDigisHandle;
        for (const auto& phDigi : phDigis){
            pulseHeight_= phDigi.energy();
            uncalibratedTime_ = phDigi.time();

            // Fill one row per PH Digi
            ttree->Fill();
            _phSpectrum->Fill(pulseHeight_);

        }

    } // end of analyze

}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::STMPHDigiTree)
