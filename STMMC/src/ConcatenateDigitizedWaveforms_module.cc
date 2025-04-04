// Concatenates STMWaveformDigis generated with HPGeWaveformsFromStepPointMCs
// Generates a summary of the number of generated concatenated waveforms and the number discarded
// Original author: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <utility>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


namespace mu2e {
  class ConcatenateDigitizedWaveforms : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> STMWaveformDigisTag{ Name("STMWaveformDigisTag"), Comment("InputTag for StepPointMCs")};
      fhicl::Atom<int> nMerge{ Name("nMerge"), Comment("Number of digitized waveforms")};
      fhicl::OptionalAtom<bool> makeTTree{ Name("makeTTree"), Comment{"Controls whether to make a TTree"}};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit ConcatenateDigitizedWaveforms(const Parameters& conf);
  private:
    void produce(art::Event& event) override;
    void endJob() override;

    // fhicl variables
    art::ProductToken<STMWaveformDigiCollection> STMWaveformDigisToken;
    int nMerge = 0; // Number of digized waveforms to concatenate

    // data product generation variables
    uint32_t outputTime = 0;
    std::vector<int16_t> outputADCs;

    // TTree variables
    bool makeTTree = false;
    TTree* ttree;
    size_t nOutputADCs = 0, i = 0;
    int16_t ADC = 0;
    uint32_t time = 0;

    // message facility variables
    uint inputEvents = 0, outputEvents = 0;
  };

  ConcatenateDigitizedWaveforms::ConcatenateDigitizedWaveforms(const Parameters& conf):
    art::EDProducer{conf},
    STMWaveformDigisToken(consumes<STMWaveformDigiCollection>(conf().STMWaveformDigisTag())),
    nMerge(conf().nMerge()) {
        produces<STMWaveformDigiCollection>();
        // Assign TTrees
        makeTTree = conf().makeTTree() ? *(conf().makeTTree()) : false;
        if (makeTTree) {
          art::ServiceHandle<art::TFileService> tfs;
          ttree = tfs->make<TTree>("ttree", "Concatenated ttree");
          ttree->Branch("time", &time, "time/i");
          ttree->Branch("ADC", &ADC, "ADC/S");
        };
      };

  void ConcatenateDigitizedWaveforms::produce(art::Event& event) {
    // Get the hits in the detector
    std::vector<STMWaveformDigi> waveforms = event.getProduct(STMWaveformDigisToken);
    inputEvents++;
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    if ((inputEvents % nMerge) == 1)
        outputTime = waveforms[0].trigTimeOffset();
    outputADCs.insert(outputADCs.end(), waveforms[0].adcs().begin(), waveforms[0].adcs().end());
    if ((inputEvents % nMerge) == 0) {
        STMWaveformDigi waveformDigi(outputTime, outputADCs);
        outputDigis->emplace_back(waveformDigi);
        outputEvents++;
        nOutputADCs = outputADCs.size();
        if (makeTTree) {
            time = outputTime;
            for (i = 0; i < nOutputADCs; i++) {
                ADC = outputADCs[i];
                ttree->Fill();
                time++;
            };
        };
        outputADCs.clear();
    };
    event.put(std::move(outputDigis));
    return;
  };

  void ConcatenateDigitizedWaveforms::endJob() {
    mf::LogInfo log("ConcatenateDigitizedWaveforms summary");
    log << "=====ConcatenateDigitizedWaveforms summary=====\n";
    log << std::left << std::setw(25) << "\tNo. input events:  " << inputEvents  << "\n";
    log << std::left << std::setw(25) << "\tNo. output events: " << outputEvents << "\n";
    log << "===============================================\n";
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::ConcatenateDigitizedWaveforms)
