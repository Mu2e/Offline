// module to write out CrvDigis and CrvStatus ntuple

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "TTree.h"

namespace mu2e {

  class CrvDigiDump : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> crvDigiCollection { Name("crvDigiTag"), Comment("CRV Digi producer name"), "CrvDigi" };
        fhicl::Atom<art::InputTag> crvStatusCollection { Name("crvStatusTag"), Comment("CRV Status producer name"), "CrvDigi" };
        fhicl::Atom<int> diagLevel { Name("diagLevel"), Comment("Diagnostic level"), 0 };
      };

      explicit CrvDigiDump(const art::EDAnalyzer::Table<Config>& config);
      virtual ~CrvDigiDump() {}

      virtual void beginJob() override;
      virtual void analyze(const art::Event& event) override;

    private:
      art::InputTag crvDigiTag_;
      art::InputTag crvStatusTag_;
      int diagLevel_;
      //int nADCsamples_;

      TTree* digiTree_;
      TTree* statusTree_;

      // Digi ntuple variables - one row per digi (flattened structure)
      Int_t evt_;
      Int_t run_;
      Int_t subrun_;
      Int_t barIndex_;      // Bar index (can be large)
      static constexpr int nADC_ = 12;  // 12 ADC samples per digi
      Short_t ADC_[nADC_]; // ADC waveform (12-bit values, use Short_t) - 2 bytes each
      UShort_t time_;       // Time in TDC units (0-65535) - 2 bytes
      UChar_t sipmNumber_;  // SiPM number (0-3) - 1 byte
      UChar_t roc_;         // ROC (1-18) - 1 byte
      UChar_t feb_;         // FEB (1-24) - 1 byte
      UChar_t febChannel_;  // FEB channel (0-63) - 1 byte
      Bool_t NZS_;          // Non-zero-suppressed flag - 1 byte
      Bool_t oddTimestamp_; // Odd timestamp flag - 1 byte

      // Status ntuple variables - one row per event, arrays up to 18 ROCs
      static constexpr int maxROCs_ = 18;

      Int_t status_evt_;
      Int_t status_run_;
      Int_t status_subrun_;
      Int_t nStatus_;                  // Number of status/ROC entries (0-18)
      Long64_t status_eventWindowTag_[maxROCs_];         // Event window tag (from data header)
      Long64_t status_subeventEWT_[maxROCs_];            // Event window tag (from subevent header)
      Long64_t status_roc_eventWindowTag_[maxROCs_];     // ROC event window tag
      Long64_t status_roc_microBunchStatus_[maxROCs_];   // Micro bunch status
      Int_t status_roc_activeFEBFlags_[maxROCs_];        // 24-bit flags

      UShort_t status_roc_controllerEventWordCount_[maxROCs_];  // Word count
      UShort_t status_roc_triggerCount_[maxROCs_];              // Trigger count

      UChar_t status_linkID_[maxROCs_];  // Link ID (0-5 typically)
      UChar_t status_dtcID_[maxROCs_];   // DTC ID (small number)
      UChar_t status_statusByte_[maxROCs_]; // Status byte (0-255)
      UChar_t status_linkStatus_[maxROCs_]; // Link status from subevent header
      UChar_t status_linkLatency_[maxROCs_]; // Link latency from subevent header
      Bool_t status_valid_[maxROCs_];    // Valid flag
  };

  CrvDigiDump::CrvDigiDump(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    crvDigiTag_(config().crvDigiCollection()),
    crvStatusTag_(config().crvStatusCollection()),
    diagLevel_(config().diagLevel()),
    //nADCsamples_(config().ADCsamples()),
    digiTree_(nullptr),
    statusTree_(nullptr)
    {}

  void CrvDigiDump::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    // Create Digi tree - flattened structure (one row per digi)
    digiTree_ = tfs->make<TTree>("CrvDigis", "CRV Digis");

    // Event information
    digiTree_->Branch("EWT", &evt_, "EWT/I");
    digiTree_->Branch("run", &run_, "run/I");
    digiTree_->Branch("subrun", &subrun_, "subrun/I");

    // Digi information (one per row)
    digiTree_->Branch("barIndex", &barIndex_, "barIndex/I");
    digiTree_->Branch("sipmNumber", &sipmNumber_, "sipmNumber/b");   // UChar_t
    digiTree_->Branch("roc", &roc_, "roc/b");                        // UChar_t
    digiTree_->Branch("feb", &feb_, "feb/b");                        // UChar_t
    digiTree_->Branch("febChannel", &febChannel_, "febChannel/b");   // UChar_t
    digiTree_->Branch("time", &time_, "time/s");                     // UShort_t
    digiTree_->Branch("NZS", &NZS_, "NZS/O");                        // Bool_t
    digiTree_->Branch("oddTimestamp", &oddTimestamp_, "oddTimestamp/O"); // Bool_t

    // ADC waveform
    digiTree_->Branch("ADC", ADC_, "ADC[12]/S");  // Short_t

    // Create Status tree - one row per event with arrays (up to 18 ROCs)
    statusTree_ = tfs->make<TTree>("CrvStatus", "CRV Status");

    statusTree_->Branch("EWT", &status_evt_, "EWT/I");
    statusTree_->Branch("run", &status_run_, "run/I");
    statusTree_->Branch("subrun", &status_subrun_, "subrun/I");
    statusTree_->Branch("nStatus", &nStatus_, "nStatus/I");  // Int_t

    // Status information (arrays indexed by ROC)
    statusTree_->Branch("valid", status_valid_, "valid[nStatus]/O");          // Bool_t
    statusTree_->Branch("linkID", status_linkID_, "linkID[nStatus]/b");       // UChar_t
    statusTree_->Branch("headerEWT", status_eventWindowTag_, "headerEWT[nStatus]/L");
    statusTree_->Branch("subeventEWT", status_subeventEWT_, "subeventEWT[nStatus]/L");
    statusTree_->Branch("dtcID", status_dtcID_, "dtcID[nStatus]/b");          // UChar_t
    statusTree_->Branch("headerStatus", status_statusByte_, "headerStatus[nStatus]/b"); // UChar_t
    statusTree_->Branch("linkStatus", status_linkStatus_, "linkStatus[nStatus]/b"); // UChar_t
    statusTree_->Branch("linkLatency", status_linkLatency_, "linkLatency[nStatus]/b"); // UChar_t

    // ROC header information (same indexing as status)
    statusTree_->Branch("wordCount", status_roc_controllerEventWordCount_, "wordCount[nStatus]/s"); // UShort_t
    statusTree_->Branch("drCount", status_roc_triggerCount_, "drCount[nStatus]/s");  // UShort_t
    statusTree_->Branch("activeFEB", status_roc_activeFEBFlags_, "activeFEB[nStatus]/I");
    statusTree_->Branch("rocStatus", status_roc_microBunchStatus_, "rocStatus[nStatus]/L");
    statusTree_->Branch("rocEWT", status_roc_eventWindowTag_, "rocEWT[nStatus]/L");

    if (diagLevel_ > 0) {
      std::cout << "CrvDigiDump: Initialized with the following settings:" << std::endl;
      std::cout << "  CrvDigi collection: " << crvDigiTag_ << std::endl;
      std::cout << "  CrvStatus collection: " << crvStatusTag_ << std::endl;
      std::cout << "  ADC samples: " << nADC_ << std::endl;
    }
  }

  void CrvDigiDump::analyze(const art::Event& event) {

    if (diagLevel_ > 1) {
      std::cout << "\nCrvDigiDump: Processing Run " << event.run()
        << ", SubRun " << event.subRun()
        << ", Event " << event.event() << std::endl;
    }

    // Get CrvDigis
    art::Handle<CrvDigiCollection> crvDigisHandle;
    event.getByLabel(crvDigiTag_, crvDigisHandle);

    int nDigis = 0;
    if (crvDigisHandle.isValid() && !crvDigisHandle->empty()) {
      const CrvDigiCollection& crvDigis = *crvDigisHandle;
      nDigis = crvDigis.size();

      if (diagLevel_ > 1) {
        std::cout << "  Found " << nDigis << " CRV digis" << std::endl;
      }

      // Fill one row per digi
      for (size_t i = 0; i < crvDigis.size(); ++i) {
        const CrvDigi& digi = crvDigis[i];

        // Event info (repeated for each digi)
        evt_ = event.event();
        run_ = event.run();
        subrun_ = event.subRun();

        // Digi info
        barIndex_ = digi.GetScintillatorBarIndex().asInt();
        sipmNumber_ = digi.GetSiPMNumber();
        roc_ = digi.GetROC();
        feb_ = digi.GetFEB();
        febChannel_ = digi.GetFEBchannel();
        time_ = digi.GetStartTDC();
        NZS_ = digi.IsNZS();
        oddTimestamp_ = digi.HasOddTimestamp();

        const std::vector<int16_t>& adcs = digi.GetADCs();
        int nSamples = std::min(static_cast<int>(adcs.size()), nADC_);
        std::copy_n(adcs.begin(), nSamples, ADC_);

        // Zero-fill remaining slots only if needed
        if (nSamples < nADC_) {
          std::fill(ADC_ + nSamples, ADC_ + nADC_, 0);
        }

        if (diagLevel_ > 2) {
          std::cout << "    Digi " << i << ": bar=" << barIndex_
            << " sipm=" << sipmNumber_
            << " roc=" << roc_
            << " feb=" << feb_
            << " febChannel=" << febChannel_
            << " TDC=" << time_
            << " nADC=" << nSamples << std::endl;
        }

        // Fill one row for this digi
        digiTree_->Fill();
      }
    } else {
      if (diagLevel_ > 1) {
        std::cout << "  No CRV digis found" << std::endl;
      }
    }

    // Get CrvStatus - one row per event with arrays
    art::Handle<CrvStatusCollection> crvStatusHandle;
    event.getByLabel(crvStatusTag_, crvStatusHandle);

    // Event info
    status_evt_ = event.event();
    status_run_ = event.run();
    status_subrun_ = event.subRun();
    nStatus_ = 0;

    if (crvStatusHandle.isValid() && !crvStatusHandle->empty()) {
      const CrvStatusCollection& crvStatus = *crvStatusHandle;
      nStatus_ = std::min(static_cast<int>(crvStatus.size()), maxROCs_);

      if (diagLevel_ > 1) {
        std::cout << "  Found " << crvStatus.size() << " CRV status entries";
        if (crvStatus.size() > static_cast<size_t>(maxROCs_)) {
          std::cout << " (truncating to " << maxROCs_ << ")";
        }
        std::cout << std::endl;
      }

      // Fill arrays
      for (int i = 0; i < nStatus_; ++i) {
        const CrvStatus& status = crvStatus[i];

        // Status info
        status_valid_[i] = status.IsValid();
        status_linkID_[i] = status.GetLinkID();
        status_eventWindowTag_[i] = status.GetEventWindowTag();
        status_subeventEWT_[i] = status.GetSubeventEventWindowTag();
        status_dtcID_[i] = status.GetDTCID();
        status_statusByte_[i] = status.GetStatus();
        status_linkStatus_[i] = status.GetLinkStatus();
        status_linkLatency_[i] = status.GetLinkLatency();

        // Get ROC header (should be one per status entry)
        const std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII>& rocHeaders =
          const_cast<CrvStatus&>(status).GetROCHeader();

        // Fill ROC header if available, otherwise set to -1
        if (!rocHeaders.empty()) {
          const CRVDataDecoder::CRVROCStatusPacketFEBII& rocHeader = rocHeaders[0];

          status_roc_controllerEventWordCount_[i] = rocHeader.ControllerEventWordCount;
          status_roc_triggerCount_[i] = rocHeader.TriggerCount;
          status_roc_activeFEBFlags_[i] = rocHeader.GetActiveFEBFlags().to_ulong();
          status_roc_microBunchStatus_[i] = rocHeader.GetMicroBunchStatus();
          status_roc_eventWindowTag_[i] = rocHeader.GetEventWindowTag();
        } else {
          // No ROC header - set sentinel values
          status_roc_controllerEventWordCount_[i] = -1;
          status_roc_triggerCount_[i] = -1;
          status_roc_activeFEBFlags_[i] = -1;
          status_roc_microBunchStatus_[i] = -1;
          status_roc_eventWindowTag_[i] = -1;
        }

        if (diagLevel_ > 2) {
          std::cout << "    Status " << i << ": linkID=" << status_linkID_[i]
            << " dtcID=" << status_dtcID_[i]
            << " EWT=" << status_eventWindowTag_[i] << std::endl;
        }
      }
    } else {
      if (diagLevel_ > 1) {
        std::cout << "  No CRV status found" << std::endl;
      }
    }

    // Fill one row per event
    statusTree_->Fill();

    if (diagLevel_ > 0 && event.event() % 100 == 0) {
      std::cout << "CrvDigiDump: Processed " << event.event() << " events" << std::endl;
    }
  }

} // namespace mu2e

using mu2e::CrvDigiDump;
DEFINE_ART_MODULE(CrvDigiDump)
