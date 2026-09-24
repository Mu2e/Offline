// =====================================================================
//
// STMBinaryDigisFromFragments: Binary File Writing
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

namespace art
{
  class STMBinaryDigisFromFragments;
}
using art::STMBinaryDigisFromFragments;
class art::STMBinaryDigisFromFragments : public EDProducer
{
public:
  struct Config {
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
    fhicl::Atom<std::string> rawHPGeFile {fhicl::Name("rawHPGeFile"), "rawHPGe.bin"};
    fhicl::Atom<std::string> zsHPGeFile {fhicl::Name("zsHPGeFile"), "zsHPGe.bin"};
    fhicl::Atom<std::string> phHPGeFile {fhicl::Name("phHPGeFile"), "phHPGe.bin"};
    fhicl::Atom<std::string> rawHeaderHPGeFile {fhicl::Name("rawHeaderHPGeFile"), "rawWithHeaderHPGe.bin"};
    fhicl::Atom<std::string> rawLaBrFile {fhicl::Name("rawLaBrFile"), "rawLaBr.bin"};
    fhicl::Atom<std::string> zsLaBrFile {fhicl::Name("zsLaBrFile"), "zsLaBr.bin"};
    fhicl::Atom<std::string> phLaBrFile {fhicl::Name("phLaBrFile"), "phLaBr.bin"};
    fhicl::Atom<std::string> rawHeaderLaBrFile {fhicl::Name("rawHeaderLaBrFile"), "rawWithHeaderLaBr.bin"};
    fhicl::Atom<std::string> eventFile {fhicl::Name("eventFile"), "event.bin"};
    fhicl::OptionalAtom<int> verbosityLevel{fhicl::Name("verbosityLevel"), fhicl::Comment("Verbosity level")};
  };

  explicit STMBinaryDigisFromFragments(const art::EDProducer::Table<Config>& config); // constructor created, config via fcl
  virtual ~STMBinaryDigisFromFragments(); //declares destructor

  virtual void produce(Event &) override;
  void endJob() override;

struct FragRecord{
  uint64_t seqID;
  uint16_t dataset;
  std::vector<int16_t> data;
};

private:
  void writeSortedSingle(std::vector<FragRecord>& buf, std::ofstream& out);
  void writeSortedCombined(std::vector<FragRecord>& buf, std::ofstream& out);
  std::string getBaseName(const std::string& fullPath);

  art::InputTag _stmFragmentsTag;
  int _verbosityLevel = 0;
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
};

// ======================================================================

STMBinaryDigisFromFragments::STMBinaryDigisFromFragments(const art::EDProducer::Table<Config>& config)
  : art::EDProducer{config}
  , _stmFragmentsTag(config().stmTag())
  , _verbosityLevel(config().verbosityLevel() ? *(config().verbosityLevel()) : 0)
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

STMBinaryDigisFromFragments::~STMBinaryDigisFromFragments(){
  if (_rawHPGeOut.is_open())  {_rawHPGeOut.close();}
  if (_zsHPGeOut.is_open()) {_zsHPGeOut.close();}
  if (_phHPGeOut.is_open()) {_phHPGeOut.close();}
  if (_rawHeaderHPGeOut.is_open()) {_rawHeaderHPGeOut.close();}

  if (_rawLaBrOut.is_open())  {_rawLaBrOut.close();}
  if (_zsLaBrOut.is_open()) {_zsLaBrOut.close();}
  if (_phLaBrOut.is_open()) {_phLaBrOut.close();}
  if (_rawHeaderLaBrOut.is_open()) {_rawHeaderLaBrOut.close();}

  if (_eventOut.is_open()) {_eventOut.close();}
} // Closing files

// ==================================================================
// Defined helper functions here

void STMBinaryDigisFromFragments::writeSortedSingle(std::vector<FragRecord>& buf,
                                              std::ofstream& out) {

  std::sort(buf.begin(), buf.end(),
            [](const FragRecord& a, const FragRecord& b) {
              return a.seqID < b.seqID;
            });

  for (const auto& rec : buf) {
    out.write(reinterpret_cast<const char*>(rec.data.data()),
              rec.data.size() * sizeof(int16_t));
  }
}

void STMBinaryDigisFromFragments::writeSortedCombined(std::vector<FragRecord>& buf,
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

  std::sort(buf.begin(), buf.end(),
            [&](const FragRecord& a, const FragRecord& b) {
              if (a.seqID != b.seqID)
                return a.seqID < b.seqID;
              return datasetOrder(a.dataset) < datasetOrder(b.dataset);
            });

  for (const auto& rec : buf) {
    out.write(reinterpret_cast<const char*>(rec.data.data()),
              rec.data.size() * sizeof(int16_t));
  }
}

std::string STMBinaryDigisFromFragments::getBaseName(const std::string& fullPath) {
  auto slash = fullPath.find_last_of("/\\");
  std::string name = (slash == std::string::npos) ? fullPath : fullPath.substr(slash + 1);

  // strip .art
  auto dot = name.rfind(".art");
  if (dot != std::string::npos) {
    name = name.substr(0, dot);
  }

  return name;
}

// ==================================================================

void STMBinaryDigisFromFragments::produce(Event& event)
{
  ++_totalEvents; //Increment Total Event Counter

  if (!_filesInitialized){
    std::ostringstream base;
    base << "artDump_run" << event.run()
         << "_subrun" << event.subRun() << "_";

    _baseName = base.str();

    // HPGe
    std::string rawHPGeName        = _baseName + _rawHPGeFile;
    std::string zsHPGeName         = _baseName + _zsHPGeFile;
    std::string phHPGeName         = _baseName + _phHPGeFile;
    std::string rawHeaderHPGeName  = _baseName + _rawHeaderHPGeFile;
    std::string eventName          = _baseName + _eventFile;

    // LaBr
    std::string rawLaBrName        = _baseName + _rawLaBrFile;
    std::string zsLaBrName         = _baseName + _zsLaBrFile;
    std::string phLaBrName         = _baseName + _phLaBrFile;
    std::string rawHeaderLaBrName  = _baseName + _rawHeaderLaBrFile;
    // HPGe
    _rawHPGeOut.open(rawHPGeName, std::ios::binary);
    _zsHPGeOut.open(zsHPGeName, std::ios::binary);
    _phHPGeOut.open(phHPGeName, std::ios::binary);
    _rawHeaderHPGeOut.open(rawHeaderHPGeName, std::ios::binary);
    _eventOut.open(eventName, std::ios::binary);
    // LaBr
    _rawLaBrOut.open(rawLaBrName, std::ios::binary);
    _zsLaBrOut.open(zsLaBrName, std::ios::binary);
    _phLaBrOut.open(phLaBrName, std::ios::binary);
    _rawHeaderLaBrOut.open(rawHeaderLaBrName, std::ios::binary);

    if(!_rawHPGeOut || !_zsHPGeOut || !_phHPGeOut || !_rawHeaderHPGeOut || !_eventOut
      ||!_rawLaBrOut || !_zsLaBrOut || !_phLaBrOut || !_rawHeaderLaBrOut ) {
        throw cet::exception("FILEOPEN") << "Failed to open output files\n";
  }
  if (_verbosityLevel > 0) {
    std::cout << "[INFO] Output files:\n"
              << "  " << rawHPGeName << "\n"
              << "  " << zsHPGeName << "\n"
              << "  " << phHPGeName << "\n"
              << "  " << rawHeaderHPGeName << "\n"
              << "  " << rawLaBrName << "\n"
              << "  " << zsLaBrName << "\n"
              << "  " << phLaBrName << "\n"
              << "  " << rawHeaderLaBrName << "\n"
              << "  " << eventName << "\n";
  }
  _filesInitialized = true;
  }

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMFragments = STMFragmentsH.product();

  for (const auto& frag : *STMFragments) {
    ++_totalFragments;

    if (_verbosityLevel >=3){ std::cout <<"\nFrag_ID : " << frag.fragmentID() << "\n";}

    //Check if this is a container fragment
    if (frag.type() == artdaq::Fragment::ContainerFragmentType) {

      mu2e::STMFragment container_frag(frag);
      artdaq::ContainerFragment cont_frag(frag);
      ++_totalContainers;
      size_t blocks = cont_frag.block_count();
      _totalInnerFrags += blocks;

      for (size_t i = 0; i < cont_frag.block_count(); ++i) {

        auto inner_frag = cont_frag.at(i);
        mu2e::STMFragment stm_frag(*inner_frag);
        const size_t physicalWords = inner_frag->dataSizeBytes() / sizeof(int16_t); // gets physical words stores in this frag

        if (stm_frag.isRaw()) {
          ++_totalRaw; //Increment job counter
          auto const* dataPtr = stm_frag.dataBegin();
          bool const isHPGe = stm_frag.isHPGe();
          bool const isLaBr = stm_frag.isLaBr();
          // header + payload
          {
            // -- optional : print first few raw words
            if (_verbosityLevel > 1) {
              std::cout << "  First 25 header words: ";
              for (size_t j = 0; j < std::min<size_t>(physicalWords,25); ++j) {
                std::cout << dataPtr[j] << " ";
              }
              std::cout << "\n";
            }

            FragRecord rec{
              inner_frag->sequenceID(),
              static_cast<uint16_t>(stm_frag.dataset()),
              std::vector<int16_t>(dataPtr, dataPtr + physicalWords)
            };

            if (isHPGe){
              ++_totalRawHPGe;
              _rawHeaderHPGeBuf.push_back(std::move(rec));
            } else if (isLaBr){
              ++_totalRawLaBr;
              _rawHeaderLaBrBuf.push_back(std::move(rec));
            }
          }

          // payload
          // interpret header
          if (physicalWords >= stm::RawHeader::WORDS) {
            auto const* payloadPtr = stm_frag.payloadBegin();
            const size_t physicalPayloadWords = physicalWords - stm::RawHeader::WORDS;
            const size_t claimedPayloadWords = stm_frag.payloadWords();
            const size_t wordsToWrite = std::min(claimedPayloadWords,physicalPayloadWords); // Pick smallest as bound
            const int16_t* hdr = dataPtr;

            int64_t ewt =
                int64_t(hdr[stm::RawHeader::EWT_0]) |
                (int64_t(hdr[stm::RawHeader::EWT_1]) << 16) |
                (int64_t(hdr[stm::RawHeader::EWT_2]) << 32);

            uint64_t raw_len = hdr[stm::RawHeader::RAW_LEN];

            // --- print debug ---
            if (_verbosityLevel > 1) {
              std::cout << "[RAW PAYLOAD DEBUG] "
                      << "SeqID=" << inner_frag->sequenceID()
                      << " EWT=" << ewt
                      << " RAW_LEN=" << raw_len
                      << " physical payload=" << physicalPayloadWords
                      << " claimed payload=" << claimedPayloadWords
                      << "\n";
            }

            FragRecord rec{
              inner_frag->sequenceID(),
              static_cast<uint16_t>(stm_frag.dataset()),
              std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)
            };

            if (isHPGe) {
              _rawHPGeBuf.push_back(std::move(rec));
            } else if (isLaBr) {
              _rawLaBrBuf.push_back(std::move(rec));
            }
          }

        }//End of isRaw
        else if (stm_frag.isZS()) {
          ++_totalZS;
          bool const isHPGe = stm_frag.isHPGe();
          bool const isLaBr = stm_frag.isLaBr();

          auto const* payloadPtr = stm_frag.payloadBegin();
          const size_t claimedWords = stm_frag.payloadWords();
          const size_t wordsToWrite = std::min(claimedWords, physicalWords);

          FragRecord rec{
            inner_frag->sequenceID(),
            static_cast<uint16_t>(stm_frag.dataset()),
            std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)
          };
          if (isHPGe) {
            ++_totalZSHPGe;
            _zsHPGeBuf.push_back(std::move(rec));
          } else if (isLaBr) {
            ++_totalZSLaBr;
            _zsLaBrBuf.push_back(std::move(rec));
          }

        } // End of isZS
        else if (stm_frag.isPH()) {
          ++_totalPH;
          bool const isHPGe = stm_frag.isHPGe();
          bool const isLaBr = stm_frag.isLaBr();

          auto const* payloadPtr = stm_frag.payloadBegin();
          const size_t claimedWords = stm_frag.payloadWords();
          const size_t wordsToWrite = std::min(claimedWords,physicalWords);

          FragRecord rec{
            inner_frag->sequenceID(),
            static_cast<uint16_t>(stm_frag.dataset()),
            std::vector<int16_t>(payloadPtr, payloadPtr + wordsToWrite)
          };

          if (isHPGe) {
            ++_totalPHHPGe;
            _phHPGeBuf.push_back(std::move(rec));
          } else if (isLaBr) {
            ++_totalPHLaBr;
            _phLaBrBuf.push_back(std::move(rec));
          }

      } else {
        // Non Raw/ZS/PH case
        ++_totalUnreadInnerFrags;
        if (_verbosityLevel > 1){
          std::cout << "\n[WARNING] Unknown Inner Frag" << std::endl;
        }
      } // End of isPH

        //---Combined stream write w. order preserved ----
        {
          const int16_t* cptr = nullptr;
          size_t cwords = 0;
          uint16_t dataset = static_cast<uint16_t>(stm_frag.dataset());

          if (stm_frag.isRaw()){
            cptr = stm_frag.dataBegin();
            cwords = physicalWords;
          }
          else if (stm_frag.isZS()){
            cptr = stm_frag.payloadBegin();
            cwords = std::min(stm_frag.payloadWords(),physicalWords);
          }
          else if (stm_frag.isPH()){
            cptr = stm_frag.payloadBegin();
            cwords = std::min(stm_frag.payloadWords(),physicalWords);
          }
          else {
            if(_verbosityLevel > 1){
              std::cout << "[WARNING] Unknown dataset : " << dataset << "\n";
            }
          }

          if (cptr && cwords > 0) {
            FragRecord rec;
            rec.seqID = inner_frag->sequenceID();
            rec.dataset = dataset;
            rec.data.assign(cptr, cptr + cwords);

            _eventBuf.push_back(std::move(rec));
          }
        }
        // end of inner fragment loop
      }
    } else {
      // Non Container type
      ++_totalNonContainers;
      if (_verbosityLevel > 1) {
        std::cout << "[WARNING] Non-Container STM Frag Type" << std::endl;
      }
    }
  }
} // produce

// ======================================================================

void STMBinaryDigisFromFragments::endJob() {
  if (_verbosityLevel > 0) {
    std::cout << "\nWriting sorted streams...\n";
  }
  writeSortedSingle(_rawHPGeBuf,       _rawHPGeOut);
  writeSortedSingle(_rawHeaderHPGeBuf, _rawHeaderHPGeOut);
  writeSortedSingle(_zsHPGeBuf,        _zsHPGeOut);
  writeSortedSingle(_phHPGeBuf,        _phHPGeOut);

  writeSortedSingle(_rawLaBrBuf,       _rawLaBrOut);
  writeSortedSingle(_rawHeaderLaBrBuf, _rawHeaderLaBrOut);
  writeSortedSingle(_zsLaBrBuf,        _zsLaBrOut);
  writeSortedSingle(_phLaBrBuf,        _phLaBrOut);

  writeSortedCombined(_eventBuf,   _eventOut);

  if(_verbosityLevel > 0){

    // Print out of summary
    std::cout << "\n========== STM JOB SUMMARY - (Binary Module) ==========\n";

    std::cout << "Total events                  : " << _totalEvents << "\n";
    std::cout << "Total fragments               : " << _totalFragments << "\n";
    std::cout << "Container frags               : " << _totalContainers << "\n";
    std::cout << "Inner fragments               : " << _totalInnerFrags << "\n";
    std::cout << "Total Non Container Fragments : " << _totalNonContainers << "\n";
    std::cout << "Unknown Inner fragments       : " << _totalUnreadInnerFrags << "\n";

    std::cout << "\n--- Data types read ---\n";
    std::cout << "RAW                           : " << _totalRaw << "\n";
    std::cout << "ZS                            : " << _totalZS << "\n";
    std::cout << "PH                            : " << _totalPH << "\n";

    std::cout << "RAW HPGE                      : " << _totalRawHPGe << "\n";
    std::cout << "ZS  HPGE                      : " << _totalZSHPGe << "\n";
    std::cout << "PH  HPGE                      : " << _totalPHHPGe << "\n";

    std::cout << "RAW LABR                      : " << _totalRawLaBr << "\n";
    std::cout << "ZS  LABR                      : " << _totalZSLaBr << "\n";
    std::cout << "PH  LABR                      : " << _totalPHLaBr << "\n";

    std::cout << "=================================\n";

  }
}

DEFINE_ART_MODULE(STMBinaryDigisFromFragments)
