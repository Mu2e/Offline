// =====================================================================
//
// STMBinaryDigisFromFragments: Binary File Writing
//
// =====================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <sstream>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

namespace art {
  class STMBinaryDigisFromFragments;
}

using art::STMBinaryDigisFromFragments;

class art::STMBinaryDigisFromFragments : public EDProducer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> stmTag{
      fhicl::Name("stmTag"),
      fhicl::Comment("Input module")
    };

    // For all file names below: an empty string disables that output.
    fhicl::Atom<std::string> rawHPGeFile{
      fhicl::Name("rawHPGeFile"), "rawHPGe.bin"
    };
    fhicl::Atom<std::string> zsHPGeFile{
      fhicl::Name("zsHPGeFile"), "zsHPGe.bin"
    };
    fhicl::Atom<std::string> phHPGeFile{
      fhicl::Name("phHPGeFile"), "phHPGe.bin"
    };
    fhicl::Atom<std::string> rawHeaderHPGeFile{
      fhicl::Name("rawHeaderHPGeFile"), "rawWithHeaderHPGe.bin"
    };

    fhicl::Atom<std::string> rawLaBrFile{
      fhicl::Name("rawLaBrFile"), "rawLaBr.bin"
    };
    fhicl::Atom<std::string> zsLaBrFile{
      fhicl::Name("zsLaBrFile"), "zsLaBr.bin"
    };
    fhicl::Atom<std::string> phLaBrFile{
      fhicl::Name("phLaBrFile"), "phLaBr.bin"
    };
    fhicl::Atom<std::string> rawHeaderLaBrFile{
      fhicl::Name("rawHeaderLaBrFile"), "rawWithHeaderLaBr.bin"
    };

    fhicl::Atom<std::string> eventFile{
      fhicl::Name("eventFile"), "event.bin"
    };

    fhicl::OptionalAtom<int> verbosityLevel{
      fhicl::Name("verbosityLevel"),
      fhicl::Comment("Verbosity level")
    };
  };

  explicit STMBinaryDigisFromFragments(
    const art::EDProducer::Table<Config>& config);

  ~STMBinaryDigisFromFragments() override;

  void produce(Event&) override;
  void endJob() override;

private:
  struct FragRecord {
    uint64_t seqID;
    uint16_t dataset;
    std::vector<int16_t> data;
  };

  void writeSortedSingle(std::vector<FragRecord>& buf, std::ofstream& out);
  void writeSortedCombined(std::vector<FragRecord>& buf, std::ofstream& out);
  void flushAndClearBuffers();

  art::InputTag _stmFragmentsTag;
  int _verbosityLevel{0};

  bool _filesInitialized{false};
  std::string _baseName;

  std::vector<FragRecord> _rawHPGeBuf;
  std::vector<FragRecord> _zsHPGeBuf;
  std::vector<FragRecord> _phHPGeBuf;
  std::vector<FragRecord> _rawHeaderHPGeBuf;

  std::vector<FragRecord> _rawLaBrBuf;
  std::vector<FragRecord> _zsLaBrBuf;
  std::vector<FragRecord> _phLaBrBuf;
  std::vector<FragRecord> _rawHeaderLaBrBuf;

  std::vector<FragRecord> _eventBuf;

  std::string _rawHPGeFile;
  std::string _rawHeaderHPGeFile;
  std::string _zsHPGeFile;
  std::string _phHPGeFile;

  std::string _rawLaBrFile;
  std::string _rawHeaderLaBrFile;
  std::string _zsLaBrFile;
  std::string _phLaBrFile;

  std::string _eventFile;

  std::ofstream _rawHPGeOut;
  std::ofstream _rawHeaderHPGeOut;
  std::ofstream _zsHPGeOut;
  std::ofstream _phHPGeOut;

  std::ofstream _rawLaBrOut;
  std::ofstream _rawHeaderLaBrOut;
  std::ofstream _zsLaBrOut;
  std::ofstream _phLaBrOut;

  std::ofstream _eventOut;

  size_t _totalEvents{0};
  size_t _totalFragments{0};
  size_t _totalContainers{0};
  size_t _totalNonContainers{0};
  size_t _totalUnreadInnerFrags{0};
  size_t _totalInnerFrags{0};

  size_t _totalRaw{0};
  size_t _totalZS{0};
  size_t _totalPH{0};

  size_t _totalRawHPGe{0};
  size_t _totalZSHPGe{0};
  size_t _totalPHHPGe{0};

  size_t _totalRawLaBr{0};
  size_t _totalZSLaBr{0};
  size_t _totalPHLaBr{0};

  // Counters for unknown cases
  size_t _totalUnknownRawHeaderDetector{0};
  size_t _totalUnknownRawPayloadDetector{0};
  size_t _totalUnknownZSDetector{0};
  size_t _totalUnknownPHDetector{0};

  // Raw fragments that are physically too short to contain a RawHeader
  size_t _totalShortRawHeader{0};
};

// =====================================================================

STMBinaryDigisFromFragments::STMBinaryDigisFromFragments(
  const art::EDProducer::Table<Config>& config)
  : art::EDProducer{config}
  , _stmFragmentsTag(config().stmTag())
  , _verbosityLevel(config().verbosityLevel()? *(config().verbosityLevel()):0)
  , _rawHPGeFile(config().rawHPGeFile())
  , _rawHeaderHPGeFile(config().rawHeaderHPGeFile())
  , _zsHPGeFile(config().zsHPGeFile())
  , _phHPGeFile(config().phHPGeFile())
  , _rawLaBrFile(config().rawLaBrFile())
  , _rawHeaderLaBrFile(config().rawHeaderLaBrFile())
  , _zsLaBrFile(config().zsLaBrFile())
  , _phLaBrFile(config().phLaBrFile())
  , _eventFile(config().eventFile())
{
  if (_rawHPGeFile.empty() && _rawHeaderHPGeFile.empty() &&
      _zsHPGeFile.empty() && _phHPGeFile.empty() &&
      _rawLaBrFile.empty() && _rawHeaderLaBrFile.empty() &&
      _zsLaBrFile.empty() && _phLaBrFile.empty() &&
      _eventFile.empty()) {
    mf::LogWarning("STMBinaryDigisFromFragments")
      << "All output file names are empty: this module will not write any files";
  }
}

// =====================================================================

STMBinaryDigisFromFragments::~STMBinaryDigisFromFragments() {
  if (_rawHPGeOut.is_open())       {_rawHPGeOut.close();}
  if (_zsHPGeOut.is_open())        {_zsHPGeOut.close();}
  if (_phHPGeOut.is_open())        {_phHPGeOut.close();}
  if (_rawHeaderHPGeOut.is_open()) {_rawHeaderHPGeOut.close();}

  if (_rawLaBrOut.is_open())       {_rawLaBrOut.close();}
  if (_zsLaBrOut.is_open())        {_zsLaBrOut.close();}
  if (_phLaBrOut.is_open())        {_phLaBrOut.close();}
  if (_rawHeaderLaBrOut.is_open()) {_rawHeaderLaBrOut.close();}

  if (_eventOut.is_open())         {_eventOut.close();}
}

// =====================================================================

void STMBinaryDigisFromFragments::writeSortedSingle(
  std::vector<FragRecord>& buf,
  std::ofstream& out) {

  if (!out.is_open()) { return; }

  std::sort(
    buf.begin(),
    buf.end(),
    [](const FragRecord& a, const FragRecord& b) {
      return a.seqID < b.seqID;
    });

  for (const auto& rec : buf) {
    out.write(
      reinterpret_cast<const char*>(rec.data.data()),
      rec.data.size() * sizeof(int16_t));
  }
}

// =====================================================================

void STMBinaryDigisFromFragments::writeSortedCombined(
  std::vector<FragRecord>& buf,
  std::ofstream& out) {

  if (!out.is_open()) { return; }

  auto datasetOrder = [](uint16_t d) {
    if (d == static_cast<uint16_t>(stm::Dataset::RAW_HPGE)) return 0;
    if (d == static_cast<uint16_t>(stm::Dataset::ZS_HPGE))  return 1;
    if (d == static_cast<uint16_t>(stm::Dataset::PH_HPGE))  return 2;
    if (d == static_cast<uint16_t>(stm::Dataset::RAW_LABR)) return 3;
    if (d == static_cast<uint16_t>(stm::Dataset::ZS_LABR))  return 4;
    if (d == static_cast<uint16_t>(stm::Dataset::PH_LABR))  return 5;
    return 99;
  };

  std::sort(
    buf.begin(),
    buf.end(),
    [&](const FragRecord& a, const FragRecord& b) {
      if (a.seqID != b.seqID) {
        return a.seqID < b.seqID;
      }
      return datasetOrder(a.dataset) < datasetOrder(b.dataset);
    });

  for (const auto& rec : buf) {
    out.write(
      reinterpret_cast<const char*>(rec.data.data()),
      rec.data.size() * sizeof(int16_t));
  }
}

// =====================================================================

void STMBinaryDigisFromFragments::flushAndClearBuffers() {
  writeSortedSingle(_rawHPGeBuf,       _rawHPGeOut);
  writeSortedSingle(_rawHeaderHPGeBuf, _rawHeaderHPGeOut);
  writeSortedSingle(_zsHPGeBuf,        _zsHPGeOut);
  writeSortedSingle(_phHPGeBuf,        _phHPGeOut);

  writeSortedSingle(_rawLaBrBuf,       _rawLaBrOut);
  writeSortedSingle(_rawHeaderLaBrBuf, _rawHeaderLaBrOut);
  writeSortedSingle(_zsLaBrBuf,        _zsLaBrOut);
  writeSortedSingle(_phLaBrBuf,        _phLaBrOut);

  writeSortedCombined(_eventBuf, _eventOut);

  _rawHPGeBuf.clear();
  _rawHeaderHPGeBuf.clear();
  _zsHPGeBuf.clear();
  _phHPGeBuf.clear();

  _rawLaBrBuf.clear();
  _rawHeaderLaBrBuf.clear();
  _zsLaBrBuf.clear();
  _phLaBrBuf.clear();

  _eventBuf.clear();
}

// =====================================================================

void STMBinaryDigisFromFragments::produce(Event& event) {
  ++_totalEvents;

  if (!_filesInitialized) {
    std::ostringstream base;
    base << "artDump_run" << event.run()
         << "_subrun" << event.subRun() << "_";

    _baseName = base.str();

    auto openIf = [&](std::ofstream& out, const std::string& suffix) {
      if (suffix.empty()) { return; }   // empty name => output disabled
      const std::string name = _baseName + suffix;
      out.open(name, std::ios::binary);
      if (!out) {
        throw cet::exception("FILEOPEN")
          << "Failed to open STM binary output file " << name << "\n";
      }
      if (_verbosityLevel > 0) {
        mf::LogInfo("STMBinaryDigisFromFragments")
          << "Output file: " << name;
      }
    };

    openIf(_rawHPGeOut,       _rawHPGeFile);
    openIf(_zsHPGeOut,        _zsHPGeFile);
    openIf(_phHPGeOut,        _phHPGeFile);
    openIf(_rawHeaderHPGeOut, _rawHeaderHPGeFile);

    openIf(_rawLaBrOut,       _rawLaBrFile);
    openIf(_zsLaBrOut,        _zsLaBrFile);
    openIf(_phLaBrOut,        _phLaBrFile);
    openIf(_rawHeaderLaBrOut, _rawHeaderLaBrFile);

    openIf(_eventOut,         _eventFile);

    _filesInitialized = true;
  }

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMFragments = STMFragmentsH.product();

  for (const auto& frag : *STMFragments) {
    ++_totalFragments;

    if (_verbosityLevel >= 3) {
      mf::LogDebug("STMBinaryDigisFromFragments")
        << "Frag_ID: " << frag.fragmentID();
    }

    if (frag.type() == artdaq::Fragment::ContainerFragmentType) {
      artdaq::ContainerFragment cont_frag(frag);

      ++_totalContainers;

      const size_t blocks = cont_frag.block_count();
      _totalInnerFrags += blocks;

      for (size_t i = 0; i < cont_frag.block_count(); ++i) {
        auto inner_frag = cont_frag.at(i);
        mu2e::STMFragment stm_frag(*inner_frag);

        const uint64_t seqID = inner_frag->sequenceID();
        const uint16_t dataset =
          static_cast<uint16_t>(stm_frag.dataset());

        // Number of int16_t words physically present in this artdaq fragment
        const size_t physicalWords =
          inner_frag->dataSizeBytes() / sizeof(int16_t);

        // Note: FragRecords (and their data copies) are only built when the
        // corresponding output is enabled, so disabled outputs cost nothing

        if (stm_frag.isRaw()) {
          ++_totalRaw;

          auto const* dataPtr = stm_frag.dataBegin();
          const bool isHPGe = stm_frag.isHPGe();
          const bool isLaBr = stm_frag.isLaBr();

          // Header + payload stream
          if (_verbosityLevel > 1) {
            std::ostringstream msg;
            msg << "First 25 raw words: ";
            for (size_t j = 0;
                 j < std::min<size_t>(physicalWords, 25);
                 ++j) {
              msg << dataPtr[j] << " ";
            }

            mf::LogDebug("STMBinaryDigisFromFragments")
              << msg.str();
          }

          if (isHPGe) {
            ++_totalRawHPGe;
            if (_rawHeaderHPGeOut.is_open()) {
              _rawHeaderHPGeBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(dataPtr, dataPtr + physicalWords)});
            }
          } else if (isLaBr) {
            ++_totalRawLaBr;
            if (_rawHeaderLaBrOut.is_open()) {
              _rawHeaderLaBrBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(dataPtr, dataPtr + physicalWords)});
            }
          } else {
            ++_totalUnknownRawHeaderDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "Raw fragment is neither HPGe nor LaBr while "
              << "dispatching raw-with-header data. "
              << "sequenceID=" << seqID;
          }

          // Payload-only stream.
          if (physicalWords >= stm::RawHeader::WORDS) {
            auto const* payloadPtr = stm_frag.payloadBegin();

            const size_t physicalPayloadWords =
              physicalWords - stm::RawHeader::WORDS;
            const size_t claimedPayloadWords = stm_frag.payloadWords();
            const size_t wordsToWrite =
              std::min(claimedPayloadWords, physicalPayloadWords);

            const int16_t* hdr = dataPtr;

            const int64_t ewt =
              int64_t(hdr[stm::RawHeader::EWT_0]) |
              (int64_t(hdr[stm::RawHeader::EWT_1]) << 16) |
              (int64_t(hdr[stm::RawHeader::EWT_2]) << 32);

            const uint64_t raw_len =
              static_cast<uint64_t>(hdr[stm::RawHeader::RAW_LEN]);

            if (_verbosityLevel > 1) {
              mf::LogDebug("STMBinaryDigisFromFragments")
                << "RAW PAYLOAD DEBUG: "
                << "SeqID=" << seqID
                << " EWT=" << ewt
                << " RAW_LEN=" << raw_len
                << " physical payload=" << physicalPayloadWords
                << " claimed payload=" << claimedPayloadWords;
            }

            if (isHPGe) {
              if (_rawHPGeOut.is_open()) {
                _rawHPGeBuf.push_back(FragRecord{
                  seqID, dataset,
                  std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
              }
            } else if (isLaBr) {
              if (_rawLaBrOut.is_open()) {
                _rawLaBrBuf.push_back(FragRecord{
                  seqID, dataset,
                  std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
              }
            } else {
              ++_totalUnknownRawPayloadDetector;
              mf::LogWarning("STMBinaryDigisFromFragments")
                << "Raw fragment is neither HPGe nor LaBr while "
                << "dispatching raw payload data. "
                << "sequenceID=" << seqID;
            }
          } else {
            ++_totalShortRawHeader;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "Raw fragment is shorter than the expected RawHeader. "
              << "sequenceID=" << seqID
              << ", physicalWords=" << physicalWords
              << ", headerWords=" << stm::RawHeader::WORDS;
          }
        }
        else if (stm_frag.isZS()) {
          ++_totalZS;

          const bool isHPGe = stm_frag.isHPGe();
          const bool isLaBr = stm_frag.isLaBr();

          auto const* payloadPtr = stm_frag.payloadBegin();
          const size_t claimedWords = stm_frag.payloadWords();
          const size_t wordsToWrite =
            std::min(claimedWords, physicalWords);

          if (isHPGe) {
            ++_totalZSHPGe;
            if (_zsHPGeOut.is_open()) {
              _zsHPGeBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
            }
          } else if (isLaBr) {
            ++_totalZSLaBr;
            if (_zsLaBrOut.is_open()) {
              _zsLaBrBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
            }
          } else {
            ++_totalUnknownZSDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "ZS fragment is neither HPGe nor LaBr. "
              << "sequenceID=" << seqID;
          }
        }
        else if (stm_frag.isPH()) {
          ++_totalPH;

          const bool isHPGe = stm_frag.isHPGe();
          const bool isLaBr = stm_frag.isLaBr();

          auto const* payloadPtr = stm_frag.payloadBegin();
          const size_t claimedWords = stm_frag.payloadWords();
          const size_t wordsToWrite =
            std::min(claimedWords, physicalWords);

          if (isHPGe) {
            ++_totalPHHPGe;
            if (_phHPGeOut.is_open()) {
              _phHPGeBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
            }
          } else if (isLaBr) {
            ++_totalPHLaBr;
            if (_phLaBrOut.is_open()) {
              _phLaBrBuf.push_back(FragRecord{
                seqID, dataset,
                std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)});
            }
          } else {
            ++_totalUnknownPHDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "PH fragment is neither HPGe nor LaBr. "
              << "sequenceID=" << seqID;
          }
        }
        else {
          ++_totalUnreadInnerFrags;

          mf::LogWarning("STMBinaryDigisFromFragments")
            << "Encountered unknown STM inner fragment. "
            << "sequenceID=" << seqID;
        }

        // Combined stream. Preserve the existing dataset order for a
        // given sequence ID
        const int16_t* cptr = nullptr;
        size_t cwords = 0;

        if (stm_frag.isRaw()) {
          cptr = stm_frag.dataBegin();
          cwords = physicalWords;
        }
        else if (stm_frag.isZS()) {
          cptr = stm_frag.payloadBegin();
          cwords = std::min(stm_frag.payloadWords(), physicalWords);
        }
        else if (stm_frag.isPH()) {
          cptr = stm_frag.payloadBegin();
          cwords = std::min(stm_frag.payloadWords(), physicalWords);
        }
        else {
          mf::LogWarning("STMBinaryDigisFromFragments")
            << "Unknown STM dataset: " << dataset
            << ", sequenceID=" << seqID;
        }

        if (cptr && cwords > 0 && _eventOut.is_open()) {
          _eventBuf.push_back(FragRecord{
            seqID, dataset,
            std::vector<int16_t>(cptr, cptr + cwords)});
        }
      }
    }
    else {
      ++_totalNonContainers;

      mf::LogWarning("STMBinaryDigisFromFragments")
        << "Encountered non-container STM fragment. "
        << "fragmentID=" << frag.fragmentID();
    }
  }

  // Flush once per art event so the buffers do not grow for the full job.
  // Sorting is still done before each write.
  flushAndClearBuffers();
}

// =====================================================================

void STMBinaryDigisFromFragments::endJob() {
  if (_verbosityLevel > 0) {
    mf::LogInfo("STMBinaryDigisFromFragments")
      << "\n========== STM JOB SUMMARY - (Binary Module) ==========\n"
      << "Total events                         : " << _totalEvents << "\n"
      << "Total fragments                      : " << _totalFragments << "\n"
      << "Container frags                      : " << _totalContainers << "\n"
      << "Inner fragments                      : " << _totalInnerFrags << "\n"
      << "Total non-container fragments        : " << _totalNonContainers << "\n"
      << "Unknown inner fragments              : " << _totalUnreadInnerFrags << "\n"
      << "\n--- Data types read ---\n"
      << "RAW                                  : " << _totalRaw << "\n"
      << "ZS                                   : " << _totalZS << "\n"
      << "PH                                   : " << _totalPH << "\n"
      << "RAW HPGe                             : " << _totalRawHPGe << "\n"
      << "ZS  HPGe                             : " << _totalZSHPGe << "\n"
      << "PH  HPGe                             : " << _totalPHHPGe << "\n"
      << "RAW LaBr                             : " << _totalRawLaBr << "\n"
      << "ZS  LaBr                             : " << _totalZSLaBr << "\n"
      << "PH  LaBr                             : " << _totalPHLaBr << "\n"
      << "\n--- Binary dispatch checks ---\n"
      << "Unknown Raw header detector          : "
      << _totalUnknownRawHeaderDetector << "\n"
      << "Unknown Raw payload detector         : "
      << _totalUnknownRawPayloadDetector << "\n"
      << "Unknown ZS detector                  : "
      << _totalUnknownZSDetector << "\n"
      << "Unknown PH detector                  : "
      << _totalUnknownPHDetector << "\n"
      << "Raw fragments shorter than header    : "
      << _totalShortRawHeader << "\n"
      << "=======================================================\n";
  }
}

DEFINE_ART_MODULE(STMBinaryDigisFromFragments)
