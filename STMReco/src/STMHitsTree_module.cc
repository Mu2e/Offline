// NTuple dumper for Detector pulse heights - calibrated

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
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TH1D.h"


// Mu2e type definitions

namespace mu2e {
    class STMHitsTree : public art::EDAnalyzer {
        public:
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          struct Config {
            fhicl::Atom<art::InputTag> stmHitsMapTag{ Name("stmHitsMapTag"), Comment("Input Tag for STMHitCollectionMap")};
            // reading from the makeSTMHits output
          };
          using Parameters = art::EDAnalyzer::Table<Config>;
          explicit STMHitsTree(const Parameters& conf);

        private:
          void analyze(const art::Event& event) override;
          void beginJob() override;
          //void endJob() override;

          art::ProductToken<STMHitCollectionMap> _stmHitCollectionMapToken; // map token
          STMChannel _channel;
          ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;

          // Store STM Hit information
          float energy   {0};
          float time     {0};

          // Store file information
          Int_t art_evt  {0};
          Int_t run      {0};
          Int_t subrun   {0};

          // Store from EventHeader
          uint64_t ewt   {0};

          // Tree reference
          TTree* ttree = nullptr;

          // STM PH spectrum
          TH1D* _energySpectrum = nullptr;
    };

    STMHitsTree::STMHitsTree(const Parameters& config) :
        art::EDAnalyzer{config},
        _stmHitCollectionMapToken(consumes<STMHitCollectionMap>(config().stmHitsMapTag())),
        _channel(STMUtils::getChannel(config().stmHitsMapTag()))
        {}

    void STMHitsTree::beginJob(){
        // Set up TTree here
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("ttree", "Detector ttree");
        ttree->Branch("energy", &energy, "energy/F");
        ttree->Branch("time", &time, "time/F");

        // Event Information
        ttree->Branch("art_event", &art_evt, "art_event/I");
        ttree->Branch("run", &run, "run/I");
        ttree->Branch("subrun", &subrun, "subrun/I");

        // EventHeader information
        ttree->Branch("EWT", &ewt, "EWT/l");

        // Set up energy spectrum
        std::string energySpectrumTitle = "Energy Spectrum (" + _channel.name() + ")" ;
        double min_energy = 0;
        double max_energy = 10;
        double energy_bin_width = 0.001;
        int n_bins = (max_energy - min_energy) / energy_bin_width;
        _energySpectrum = tfs->make<TH1D>(
            "energySpectrum",
            (energySpectrumTitle +";Energy;Count").c_str(),
            n_bins, min_energy, max_energy);
    };

    void STMHitsTree::analyze(const art::Event& event) {
        // We fill art based information here
        art_evt  = event.event();
        run      = event.run();
        subrun   = event.subRun();

        // Get handle for Hit Collection Map
        auto stmHitCollectionMapHandle = event.getValidHandle(_stmHitCollectionMapToken);
        for (const auto& mu2e_evt : *stmHitCollectionMapHandle) {
            // get EventHeader
            const auto& header = mu2e_evt.first;

            // store related info
            ewt = header.eventWindowTag();

            // get Hit Collection
            const auto& stmHits = mu2e_evt.second;

            // Second loop for Hit Collection
            for (const auto& stmHit : stmHits){
                // get energy and time from hit collection
                energy = stmHit.energy();
                time   = stmHit.time();
                // Fill tree with EWT and calibrated (energy,time)
                ttree->Fill();
                // fill histogram
                _energySpectrum->Fill(energy);
            }
        }
    } // end of analyze
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::STMHitsTree)
