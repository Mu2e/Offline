// =====================================================================
//
// STMBinaryDigisFromFragments: Binary File Writing
//
// ======================================================================

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
{}

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

    const std::string rawHPGeName       = _baseName + _rawHPGeFile;
    const std::string zsHPGeName        = _baseName + _zsHPGeFile;
    const std::string phHPGeName        = _baseName + _phHPGeFile;
    const std::string rawHeaderHPGeName = _baseName + _rawHeaderHPGeFile;

    const std::string rawLaBrName       = _baseName + _rawLaBrFile;
    const std::string zsLaBrName        = _baseName + _zsLaBrFile;
    const std::string phLaBrName        = _baseName + _phLaBrFile;
    const std::string rawHeaderLaBrName = _baseName + _rawHeaderLaBrFile;

    const std::string eventName         = _baseName + _eventFile;

    _rawHPGeOut.open(rawHPGeName, std::ios::binary);
    _zsHPGeOut.open(zsHPGeName, std::ios::binary);
    _phHPGeOut.open(phHPGeName, std::ios::binary);
    _rawHeaderHPGeOut.open(rawHeaderHPGeName, std::ios::binary);

    _rawLaBrOut.open(rawLaBrName, std::ios::binary);
    _zsLaBrOut.open(zsLaBrName, std::ios::binary);
    _phLaBrOut.open(phLaBrName, std::ios::binary);
    _rawHeaderLaBrOut.open(rawHeaderLaBrName, std::ios::binary);

    _eventOut.open(eventName, std::ios::binary);

    if (!_rawHPGeOut || !_zsHPGeOut || !_phHPGeOut ||
        !_rawHeaderHPGeOut || !_rawLaBrOut || !_zsLaBrOut ||
        !_phLaBrOut || !_rawHeaderLaBrOut || !_eventOut) {
      throw cet::exception("FILEOPEN")
        << "Failed to open one or more STM binary output files\n";
    }

    if (_verbosityLevel > 0) {
      mf::LogInfo("STMBinaryDigisFromFragments")
        << "Output files:\n"
        << "  " << rawHPGeName << "\n"
        << "  " << zsHPGeName << "\n"
        << "  " << phHPGeName << "\n"
        << "  " << rawHeaderHPGeName << "\n"
        << "  " << rawLaBrName << "\n"
        << "  " << zsLaBrName << "\n"
        << "  " << phLaBrName << "\n"
        << "  " << rawHeaderLaBrName << "\n"
        << "  " << eventName;
    }

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

        // Number of int16_t words physically present in this artdaq fragment.
        const size_t physicalWords =
          inner_frag->dataSizeBytes() / sizeof(int16_t);

        if (stm_frag.isRaw()) {
          ++_totalRaw;

          auto const* dataPtr = stm_frag.dataBegin();
          const bool isHPGe = stm_frag.isHPGe();
          const bool isLaBr = stm_frag.isLaBr();

          // Header + payload stream.
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

          FragRecord rawHeaderRec{
            inner_frag->sequenceID(),
            static_cast<uint16_t>(stm_frag.dataset()),
            std::vector<int16_t>(dataPtr, dataPtr + physicalWords)
          };

          if (isHPGe) {
            ++_totalRawHPGe;
            _rawHeaderHPGeBuf.push_back(std::move(rawHeaderRec));
          } else if (isLaBr) {
            ++_totalRawLaBr;
            _rawHeaderLaBrBuf.push_back(std::move(rawHeaderRec));
          } else {
            ++_totalUnknownRawHeaderDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "Raw fragment is neither HPGe nor LaBr while "
              << "dispatching raw-with-header data. "
              << "sequenceID=" << inner_frag->sequenceID();
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
                << "SeqID=" << inner_frag->sequenceID()
                << " EWT=" << ewt
                << " RAW_LEN=" << raw_len
                << " physical payload=" << physicalPayloadWords
                << " claimed payload=" << claimedPayloadWords;
            }

            FragRecord rawRec{
              inner_frag->sequenceID(),
              static_cast<uint16_t>(stm_frag.dataset()),
              std::vector<int16_t>(
                payloadPtr,
                payloadPtr + wordsToWrite)
            };

            if (isHPGe) {
              _rawHPGeBuf.push_back(std::move(rawRec));
            } else if (isLaBr) {
              _rawLaBrBuf.push_back(std::move(rawRec));
            } else {
              ++_totalUnknownRawPayloadDetector;
              mf::LogWarning("STMBinaryDigisFromFragments")
                << "Raw fragment is neither HPGe nor LaBr while "
                << "dispatching raw payload data. "
                << "sequenceID=" << inner_frag->sequenceID();
            }
          } else {
            ++_totalShortRawHeader;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "Raw fragment is shorter than the expected RawHeader. "
              << "sequenceID=" << inner_frag->sequenceID()
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

          FragRecord rec{
            inner_frag->sequenceID(),
            static_cast<uint16_t>(stm_frag.dataset()),
            std::vector<int16_t>(
              payloadPtr,
              payloadPtr + wordsToWrite)
          };

          if (isHPGe) {
            ++_totalZSHPGe;
            _zsHPGeBuf.push_back(std::move(rec));
          } else if (isLaBr) {
            ++_totalZSLaBr;
            _zsLaBrBuf.push_back(std::move(rec));
          } else {
            ++_totalUnknownZSDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "ZS fragment is neither HPGe nor LaBr. "
              << "sequenceID=" << inner_frag->sequenceID();
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

          FragRecord rec{
            inner_frag->sequenceID(),
            static_cast<uint16_t>(stm_frag.dataset()),
            std::vector<int16_t>(
              payloadPtr,
              payloadPtr + wordsToWrite)
          };

          if (isHPGe) {
            ++_totalPHHPGe;
            _phHPGeBuf.push_back(std::move(rec));
          } else if (isLaBr) {
            ++_totalPHLaBr;
            _phLaBrBuf.push_back(std::move(rec));
          } else {
            ++_totalUnknownPHDetector;
            mf::LogWarning("STMBinaryDigisFromFragments")
              << "PH fragment is neither HPGe nor LaBr. "
              << "sequenceID=" << inner_frag->sequenceID();
          }
        }
        else {
          ++_totalUnreadInnerFrags;

          mf::LogWarning("STMBinaryDigisFromFragments")
            << "Encountered unknown STM inner fragment. "
            << "sequenceID=" << inner_frag->sequenceID();
        }

        // Combined stream. Preserve the existing dataset order for a
        // given sequence ID.
        const int16_t* cptr = nullptr;
        size_t cwords = 0;
        const uint16_t dataset =
          static_cast<uint16_t>(stm_frag.dataset());

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
            << ", sequenceID=" << inner_frag->sequenceID();
        }

        if (cptr && cwords > 0) {
          FragRecord rec{
            inner_frag->sequenceID(),
            dataset,
            std::vector<int16_t>(cptr, cptr + cwords)
          };

          _eventBuf.push_back(std::move(rec));
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
