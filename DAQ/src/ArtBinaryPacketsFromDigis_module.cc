//
// EDProducer module for converting tracker/calo/crv digis
// into DTC formatted packets
//
//
// Original author G. Pezzullo, E. Flumerfelt, and R. Ehrlich
// V2 version with Calo DMAP handling from S.Miscetti
//

// C++ includes.
#include <cmath>
#include <iostream>
#include <string>

#include <math.h>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// geometry
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

// artdaq-core-mu2e includes
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"

// pci_linux_kernel_module includes
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_Event.h"

// Mu2e includes.
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
// #include "Offline/DAQDataProducts/inc/DataBlockCollection.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/SeedService/inc/SeedService.hh"

#include "TRACE/tracemf.h"

#include <fstream>
#include <stdexcept>

// Typedefs needed for compatibility with HLS codeblock
using timestamp = uint64_t;
struct DataBlockHeader // From mu2e_pcie_utils mu2e_mmap_ioctl.h TODO: Use DTC_DataHeaderPacket
                       // instead!
{
  uint16_t TransferByteCount; ///< Block Byte count

  uint16_t Resv1 : 4;      ///< Reserved
  uint16_t PacketType : 4; ///< Type of packet
  uint16_t LinkID : 3;     ///< Link ID of packet
  uint16_t DTCErrors : 4;  ///< Subsystem ID
  uint16_t Valid : 1;      ///< Is the packet valid?

  uint16_t PacketCount : 11;    ///< Packet count requested
  uint16_t Resv2 : 2;           ///< Reserved
  uint16_t SubsystemID : 3;     ///< Subsystem of Data Block
  uint16_t ts10;                ///< Timestamp bytes 1 and 2 (Least significant)
  uint16_t ts32;                ///< Timestamp bytes 3 and 4
  uint16_t ts54;                ///< Timestamp bytes 5 and 6 (Most significant)
  uint16_t Status : 8;          ///< Status word
  uint16_t Version : 8;         ///< Data packet format version
  uint16_t DTCID : 8;           ///< ID of receiving DTC
  uint16_t EventWindowMode : 8; ///< Window mode byte from CFO
};

using TrackerDataPacket = mu2e::TrackerDataDecoder::TrackerDataPacket;
using TrackerADCPacket = mu2e::TrackerDataDecoder::TrackerADCPacket;
using adc_t = uint16_t;

// adding in new calorimter packet
using CalorimeterHitDataPacket = mu2e::CalorimeterDataDecoder::CalorimeterHitDataPacket;
using CalorimeterFooterPacket = mu2e::CalorimeterDataDecoder::CalorimeterFooterPacket;
using CRVROCStatusPacket = mu2e::CRVDataDecoder::CRVROCStatusPacket;
using CRVHitWaveformSample = mu2e::CRVDataDecoder::CRVHitWaveformSample;
using CRVHitInfo = mu2e::CRVDataDecoder::CRVHitInfo;
using CRVHit = mu2e::CRVDataDecoder::CRVHit;

// data struct for the calorimeter
struct CaloDataPacket {
  CalorimeterHitDataPacket dataPacket;
  std::vector<CalorimeterHitDataPacket> hitPacketVec;
  CalorimeterFooterPacket dataFooterPacket;
  std::vector<std::vector<adc_t>> waveformVec;
  std::vector<uint16_t> hitIndex;
};

// data struct for the tracker
struct TrackerFullHitFormat {
  TrackerDataPacket mainPacket;
  std::vector<TrackerADCPacket> adcPacketVec;
};

using tracker_data_block_t = std::pair<DataBlockHeader, std::vector<TrackerFullHitFormat>>;
using calo_data_block_t = std::pair<DataBlockHeader, CaloDataPacket>;

using tracker_data_block_list_t = std::deque<tracker_data_block_t>;
using calo_data_block_list_t = std::deque<calo_data_block_t>;

// data struct for the crv
struct CrvDataPacket {
  DataBlockHeader header;
  CRVROCStatusPacket rocStatus;
  std::vector<CRVHit> hits;

  CrvDataPacket() : rocStatus(), hits() { bzero(&header, sizeof(header)); }
};

using crv_data_block_list_t = std::map<int, CrvDataPacket>; // the map key is the CRV ROC ID

namespace mu2e {

constexpr int format_version = 1;
constexpr int NUM_PRESAMPLES = 4;
constexpr int START_SAMPLES = 5; // 0 indexed
constexpr int NUM_SAMPLES = 15;
constexpr int LOWER_TDC = 16000;
constexpr int UPPER_TDC = 64000;

//--------------------------------------------------------------------
//
//
class ArtBinaryPacketsFromDigis : public art::EDProducer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<size_t> generateTimestampTable{Name("generateTimestampTable"),
                                               Comment("generate Timestamp table")};
    fhicl::Atom<std::string> tableFile{Name("tableFile"), Comment("table file name")};
    fhicl::Atom<size_t> timestampOffset{Name("timestampOffset"), Comment("timestamp Offset")};
    fhicl::Atom<int> includeTracker{Name("includeTracker"), Comment("include Tracker digis")};
    fhicl::Atom<int> includeCalorimeter{Name("includeCalorimeter"),
                                        Comment("include Calorimeter digis")};
    fhicl::Atom<int> includeCrv{Name("includeCrv"), Comment("include Crv digis")};
    fhicl::Atom<int> includeDMAHeaders{Name("includeDMAHeaders"), Comment("include DMA Headers")};
    fhicl::Atom<int> generateBinaryFile{Name("generateBinaryFile"), Comment("generate BinaryFile")};
    fhicl::Atom<std::string> outputFile{Name("outputFile"), Comment("output File name")};
    fhicl::Atom<int> generateTextFile{Name("generateTextFile"), Comment("generate Text File")};
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagLevel")};
    fhicl::Atom<int> maxFullPrint{Name("maxFullPrint"), Comment("maxFullPrint")};
    fhicl::Atom<art::InputTag> sdtoken{Name("strawDigiCollection"),
                                       Comment("Straw digi collection name")};
    fhicl::Atom<art::InputTag> cdtoken{Name("caloDigiCollection"),
                                       Comment("Calo digi collection name")};
    fhicl::Atom<art::InputTag> crvtoken{Name("crvDigiCollection"),
                                        Comment("Crv digi collection name")};
  };

  explicit ArtBinaryPacketsFromDigis(const art::EDProducer::Table<Config>& config);

  virtual void beginJob() override;
  virtual void beginRun(art::Run&) override;

  virtual void endJob() override;

  void produce(art::Event&) override;

private:
  size_t _generateTimestampTable;

  // Table used for mapping between DTC timestamp and art EventID
  std::string _tableFile;
  std::vector<std::pair<timestamp, timestamp>> tsTable;

  size_t _timestampOffset;
  int _includeTracker;
  int _includeCalorimeter;
  int _includeCrv;
  int _includeDMAHeaders;

  // -- include proditions handling
  ProditionsHandle<CaloDAQMap> _calodaqconds_h;
  ProditionsHandle<CRVOrdinal> _crvChannelMap_h;
  // Set to 1 to save packet data to a binary file
  int _generateBinaryFile;

  std::string _outputFile;
  std::ofstream outputStream;

  //--------------------------------------------------------------------------------
  // TRACKER ROC/DTC INFO
  //--------------------------------------------------------------------------------
  // 96 straws per panel
  // 1 ROC per panel
  // 216 panels
  //
  // 6 ROCs per DTC
  // 36 DTCs
  //    const size_t number_of_rocs = 216;
  const size_t number_of_rocs = 240;
  const size_t number_of_straws_per_roc = 96; // Each panel in the tracker has 96 straws
  const size_t number_of_rocs_per_dtc = 6;
  const size_t numADCSamples = 15;

  //--------------------------------------------------------------------------------
  // CALORIEMTER ROC/DTC INFO
  //--------------------------------------------------------------------------------
  // 6 rocs per DTC => 23 DTCs
  // 136 rocs * 10 crystals per roc => 1348
  const size_t number_of_calo_rocs = 136;
  const size_t number_of_crystals_per_roc = 10;
  const size_t number_of_calo_rocs_per_dtc = 6;

  //--------------------------------------------------------------------------------
  // CRV ROC INFO
  //--------------------------------------------------------------------------------

  const size_t number_of_crv_rocs = 17;

  //--------------------------------------------------------------------------------

  int _generateTextFile;

  // Label of the module that made the hits.
  art::ProductToken<StrawDigiCollection> const _sdtoken;
  art::ProductToken<StrawDigiADCWaveformCollection> const _sdadctoken;
  art::ProductToken<CaloDigiCollection> const _cdtoken;
  art::ProductToken<CrvDigiCollection> const _crvtoken;

  size_t _numWordsWritten;
  size_t _numEventsProcessed;

  const Calorimeter* _calorimeter; // cached pointer to the calorimeter geometry

  void fillEmptyHeaderDataPacket(DataBlockHeader& HeaderData, uint64_t& EventNum, uint8_t& ROCId,
                                 uint8_t& DTCId, uint8_t Subsys);
  void printHeader(DataBlockHeader const& headerDataBlock);

  void putBlockInEvent(DTCLib::DTC_Event& currentEvent, uint8_t dtcID, DTCLib::DTC_Subsystem subsys,
                       DTCLib::DTC_DataBlock thisBlock) {
    auto subEvt = currentEvent.GetSubEventByDTCID(dtcID, subsys);
    if (subEvt == nullptr) {
      DTCLib::DTC_SubEvent newSubEvt;
      newSubEvt.SetEventWindowTag(currentEvent.GetEventWindowTag());
      newSubEvt.SetSourceDTC(dtcID, subsys);
      newSubEvt.AddDataBlock(thisBlock);
      currentEvent.AddSubEvent(newSubEvt);

    } else {

      subEvt->AddDataBlock(thisBlock);
    }
  }

  //--------------------------------------------------------------------------------
  //  methods used to process the tracker data
  //--------------------------------------------------------------------------------
  void fillTrackerDataPacket(const StrawDigi& SD, const StrawDigiADCWaveform& SDADC,
                             TrackerFullHitFormat& TrkData, DataBlockHeader& headerData);

  void processTrackerData(art::Event& evt, uint64_t& eventNum,
                          tracker_data_block_list_t& trackerData);

  void fillTrackerDMABlocks(DTCLib::DTC_Event& currentEvent,
                            tracker_data_block_list_t const& trackerData);

  void fillTrackerDataStream(DTCLib::DTC_Event& currentEvent, tracker_data_block_t const& trkData);

  void printTrackerData(std::vector<TrackerFullHitFormat> const& curDataBlock);

  //--------------------------------------------------------------------------------
  //  methods used to handle the calorimeter data
  //--------------------------------------------------------------------------------
  void fillCalorimeterDataPacket(CaloDAQMap const& calodaqconds, const CaloDigi& SD,
                                 CaloDataPacket& caloData);

  void addCaloHitToCaloPacket(calo_data_block_t& dataBlock, CaloDataPacket& caloData);

  void fillCalorimeterHeaderDataPacket(CaloDAQMap const& calodaqconds, const CaloDigi& SD,
                                       DataBlockHeader& HeaderData, uint64_t& EventNum);

  void fillHeaderByteAndPacketCounts(calo_data_block_t& caloData);

  void processCalorimeterData(art::Event& evt, uint64_t& eventNum,
                              calo_data_block_list_t& caloDataBlocks);

  void fillCalorimeterDMABlocks(DTCLib::DTC_Event& currentEvent, calo_data_block_list_t& caloData);

  void fillCalorimeterDataStream(DTCLib::DTC_Event& currentEvent, calo_data_block_t& caloData);

  void printCalorimeterData(CaloDataPacket const& curDataBlock);

  const size_t waveformMaximumIndex(std::vector<adc_t> const& waveform);

  //--------------------------------------------------------------------------------
  //  methods used to handle the crv data
  //--------------------------------------------------------------------------------

  void processCrvData(art::Event& evt, uint64_t& eventNum, crv_data_block_list_t& crvDataBlocks);
  //  uint8_t compressCrvDigi(int adc);
  int16_t compressCrvDigi(int16_t adc);
  void fillCrvDataPacket(const CRVOrdinal& crvChannelMap, const CrvDigi& digi, CRVHit& hit,
                         int& rocID);
  void fillCrvHeaderPacket(const CRVOrdinal& crvChannelMap, CrvDataPacket& crvData, uint8_t rocID,
                           uint64_t eventNum);
  void fillCrvDMABlocks(DTCLib::DTC_Event& currentEvent, const crv_data_block_list_t& crvData);
  void fillCrvDataStream(DTCLib::DTC_Event& currentEvent, const CrvDataPacket& crvData);
  void printCrvData(const CrvDataPacket& curDataBlock);

  //--------------------------------------------------------------------------------
};

//--------------------------------------------------------------------------------
// temporary function used to find the location of the waveform peak in the
// calorimeter digitized waveform
//--------------------------------------------------------------------------------
const size_t ArtBinaryPacketsFromDigis::waveformMaximumIndex(std::vector<adc_t> const& waveform) {
  size_t indexMax(0), content(0);
  for (size_t i = 0; i < waveform.size(); ++i) {
    if (waveform[i] > content) {
      content = waveform[i];
      indexMax = i;
    }
  }

  return indexMax;
}

void ArtBinaryPacketsFromDigis::printHeader(DataBlockHeader const& headerDataBlock) {
  TLOG(TLVL_DEBUG + 12) << "START header print";
  TLOG(TLVL_DEBUG + 12) << "ByteCount     : " << headerDataBlock.TransferByteCount;
  TLOG(TLVL_DEBUG + 12) << "Resv1         : " << headerDataBlock.Resv1;
  TLOG(TLVL_DEBUG + 12) << "PacketType    : " << headerDataBlock.PacketType;
  TLOG(TLVL_DEBUG + 12) << "ROCID         : " << headerDataBlock.LinkID;
  TLOG(TLVL_DEBUG + 12) << "SubsystemID   : " << headerDataBlock.SubsystemID;
  TLOG(TLVL_DEBUG + 12) << "Valid         : " << headerDataBlock.Valid;
  TLOG(TLVL_DEBUG + 12) << "PacketCount   : " << headerDataBlock.PacketCount;
  TLOG(TLVL_DEBUG + 12) << "Resv2         : " << headerDataBlock.Resv2;
  TLOG(TLVL_DEBUG + 12) << "ts10          : " << headerDataBlock.ts10;
  TLOG(TLVL_DEBUG + 12) << "ts32          : " << headerDataBlock.ts32;
  TLOG(TLVL_DEBUG + 12) << "ts54          : " << headerDataBlock.ts54;
  TLOG(TLVL_DEBUG + 12) << "Status        : " << headerDataBlock.Status;
  TLOG(TLVL_DEBUG + 12) << "FormatVersion : " << headerDataBlock.Version;
  TLOG(TLVL_DEBUG + 12) << "DTCID         : " << headerDataBlock.DTCID;
  TLOG(TLVL_DEBUG + 12) << "EVBMode       : " << headerDataBlock.EventWindowMode;
}

void ArtBinaryPacketsFromDigis::printTrackerData(std::vector<TrackerFullHitFormat> const& trkData) {
  TLOG(TLVL_DEBUG + 13) << "START tracker-data print";
  for (size_t i = 0; i < trkData.size(); i++) {
    TLOG(TLVL_DEBUG + 13) << "StrawIndex : " << (int)trkData[i].mainPacket.StrawIndex;
    TLOG(TLVL_DEBUG + 13) << "TDC0       : " << (int)trkData[i].mainPacket.TDC0();
    TLOG(TLVL_DEBUG + 13) << "TDC1       : " << (int)trkData[i].mainPacket.TDC1();
    TLOG(TLVL_DEBUG + 13) << "TOT0       : " << (int)trkData[i].mainPacket.TOT0;
    TLOG(TLVL_DEBUG + 13) << "TOT1       : " << (int)trkData[i].mainPacket.TOT1;
    TLOG(TLVL_DEBUG + 13) << "PMP        : " << (int)trkData[i].mainPacket.PMP;
    TLOG(TLVL_DEBUG + 13) << "ADC00      : " << (int)trkData[i].mainPacket.ADC00;
    TLOG(TLVL_DEBUG + 13) << "ADC01      : " << (int)trkData[i].mainPacket.ADC01();
    TLOG(TLVL_DEBUG + 13) << "ADC02      : " << (int)trkData[i].mainPacket.ADC02;
    TLOG(TLVL_DEBUG + 13) << "ErrorFlags : " << (int)trkData[i].mainPacket.ErrorFlags;
  }
}

void ArtBinaryPacketsFromDigis::printCalorimeterData(CaloDataPacket const& caloData) {
  // CalorimeterHitDataPacket packet = caloData.dataPacket;
  CalorimeterFooterPacket footerpacket = caloData.dataFooterPacket;
  size_t nHits = caloData.hitPacketVec.size();
  TLOG(TLVL_DEBUG + 14) << "START calorimeter-data print";
  TLOG(TLVL_DEBUG + 14) << "NumberofHits        : " << (int)nHits;
  TLOG(TLVL_DEBUG + 14) << "BoardID             : " << (int)footerpacket.BoardID;
  TLOG(TLVL_DEBUG + 14) << "ChannelStatusFlagsA : " << (int)footerpacket.ChannelStatusFlagA;
  TLOG(TLVL_DEBUG + 14) << "ChannelStatusFlagsB : " << (int)footerpacket.ChannelStatusFlagC;
  TLOG(TLVL_DEBUG + 14) << "unused              : " << (int)footerpacket.unused;

  for (size_t i = 0; i < nHits; ++i) {
    CalorimeterHitDataPacket const& hit = caloData.hitPacketVec[i];
    TLOG(TLVL_DEBUG + 14) << "hit           : " << (int)i;
    TLOG(TLVL_DEBUG + 14) << "ChannelNumber : " << (int)hit.ChannelNumber;
    TLOG(TLVL_DEBUG + 14) << "DIRACA        : " << (int)hit.DIRACA;
    TLOG(TLVL_DEBUG + 14) << "DIRACB        : " << (int)hit.DIRACB;
    TLOG(TLVL_DEBUG + 14) << "ErrorFlags    : " << (int)hit.ErrorFlags;
    TLOG(TLVL_DEBUG + 14) << "Time          : " << (int)hit.Time;
    // TLOG(TLVL_DEBUG + 14) << "NumberOfSamples : " << (int)hit.NumberOfSamples;// TODO
    TLOG(TLVL_DEBUG + 14) << "IndexOfMaxDigitizerSample : " << (int)hit.IndexOfMaxDigitizerSample;
  }
}

void ArtBinaryPacketsFromDigis::printCrvData(CrvDataPacket const& crvData) {
  size_t nHits = crvData.hits.size();
  TLOG(TLVL_DEBUG + 5) << "START crv-data print";
  TLOG(TLVL_DEBUG + 5) << "ROC controller ID   : " << (int)crvData.rocStatus.ControllerID;
  TLOG(TLVL_DEBUG + 5) << "NHits               : " << (int)nHits;

  for (size_t i = 0; i < nHits; ++i) {
    TLOG(TLVL_DEBUG + 5) << "hit           : " << (int)i;
    TLOG(TLVL_DEBUG + 5) << "Channel       : " << (int)crvData.hits[i].first.febChannel;
    TLOG(TLVL_DEBUG + 5) << "FEB           : " << (int)crvData.hits[i].first.portNumber;
    TLOG(TLVL_DEBUG + 5) << "Time          : " << (int)crvData.hits[i].first.HitTime;
    TLOG(TLVL_DEBUG + 5) << "NumOfSamples  : " << (int)crvData.hits[i].first.NumSamples;
  }
}

void ArtBinaryPacketsFromDigis::fillTrackerDataStream(DTCLib::DTC_Event& currentEvent,
                                                      tracker_data_block_t const& trackerData) {

  auto sz = sizeof(DataBlockHeader);
  // check that the trkDataBlock is not empty

  if (trackerData.first.PacketCount > 0) {
    sz += sizeof(TrackerDataPacket) * trackerData.first.PacketCount;
  }

  uint8_t dtcID = trackerData.first.DTCID;
  DTCLib::DTC_DataBlock thisBlock(sz);

  if (thisBlock.blockPointer == nullptr) {
    throw cet::exception("MemoryAllocationError")
        << "Unable to allocate memory for Tracker block! sz=" << sz;
  }

  auto pos = 0;
  memcpy(thisBlock.allocBytes->data(), &trackerData.first, sizeof(DataBlockHeader));
  pos += sizeof(DataBlockHeader);

  if (trackerData.first.PacketCount > 0) {
    if (sizeof(TrackerDataPacket) % 16 !=
        0) { // Make sure that TrackerDataPacket is an even number of DataPackets!
      throw cet::exception("Online-RECO") << "ArtBinaryPacketsFromDigis::fillTrackerDataStream : "
                                             "sizeof(TrackerDataPacket) % 16 == 0";
    }
    for (size_t ipkt = 0; ipkt < trackerData.second.size(); ipkt++) {
      auto ptr = &(trackerData.second[ipkt]);
      size_t num_packets = ptr->mainPacket.NumADCPackets;
      memcpy(thisBlock.allocBytes->data() + pos, &(ptr->mainPacket), sizeof(TrackerDataPacket));
      pos += sizeof(TrackerDataPacket);
      if (trackerData.second[ipkt].adcPacketVec.size() > 0) {
        memcpy(thisBlock.allocBytes->data() + pos, &(ptr->adcPacketVec[0]),
               sizeof(TrackerADCPacket) * num_packets);
        pos += sizeof(TrackerADCPacket) * num_packets;
      }
    }
  }
  putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_Tracker, thisBlock);
}

void ArtBinaryPacketsFromDigis::fillTrackerDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                     tracker_data_block_list_t const& trkData) {

  auto curDTCID = trkData.front().first.DTCID;
  bool first = true;

  TLOG(TLVL_DEBUG + 1) << "trkData.size() = " << trkData.size();

  for (auto& dataBlock : trkData) {

    fillTrackerDataStream(currentEvent, dataBlock);

    if (dataBlock.first.DTCID != curDTCID || first) {
      TLOG(TLVL_DEBUG + 1) << "\t\tDTCID: " << (int)dataBlock.first.DTCID;
      TLOG(TLVL_DEBUG + 1) << "\t\tSYSID: " << (int)dataBlock.first.SubsystemID;
      curDTCID = dataBlock.first.DTCID;
      first = false;
    }
    if (dataBlock.first.PacketCount > 0) {
      printHeader(dataBlock.first);
      printTrackerData(dataBlock.second);
    }

  } // End loop over DataBlocks
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillEmptyHeaderDataPacket(DataBlockHeader& headerData,
                                                          uint64_t& EventNum, uint8_t& ROCId,
                                                          uint8_t& DTCId, uint8_t Subsys) {

  bzero(&headerData, sizeof(DataBlockHeader));
  // Fill in the byte count field of the header packet
  // Word 0
  headerData.TransferByteCount = sizeof(DataBlockHeader);
  // Word 1
  headerData.Resv1 = 0;      // ask Eric!!!//FIX ME!
  headerData.PacketType = 5; // PacketType::Dataheader;

  headerData.LinkID = ROCId;

  headerData.SubsystemID = Subsys; // DTCLib::DTC_Subsystem_Tracker; //: 3;

  headerData.Valid = 1;
  // Word 2
  headerData.PacketCount = 0;
  headerData.Resv2 = 0; // : 5;
  // Word 3
  uint64_t timestamp = EventNum;
  headerData.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  headerData.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  headerData.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  headerData.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  headerData.Version = format_version;
  // Word 7
  headerData.DTCID = DTCId;
  uint8_t evbMode = 0; // maybe off-spill vs on-spill?
  headerData.EventWindowMode = evbMode;
}

void ArtBinaryPacketsFromDigis::fillTrackerDataPacket(const StrawDigi& SD,
                                                      const StrawDigiADCWaveform& SDADC,
                                                      TrackerFullHitFormat& TrkData,
                                                      DataBlockHeader& headerData) {

  TrkData.mainPacket.StrawIndex = SD.strawId().asUint16();
  TrkData.mainPacket.SetTDC0(SD.TDC(StrawEnd::cal));
  TrkData.mainPacket.SetTDC1(SD.TDC(StrawEnd::hv));
  TrkData.mainPacket.TOT0 = SD.TOT(StrawEnd::cal);
  TrkData.mainPacket.TOT1 = SD.TOT(StrawEnd::hv);
  TrkData.mainPacket.EWMCounter = headerData.ts10 & 0xF;
  TrkData.mainPacket.PMP = SD.PMP();
  TrkData.mainPacket.ErrorFlags = 0; // FIXME
  TrkData.mainPacket.unused1 = 0;

  headerData.TransferByteCount += sizeof(TrackerDataPacket);
  headerData.PacketCount++;

  TrkTypes::ADCWaveform const& theWaveform = SDADC.samples();
  size_t numADCPackets = static_cast<size_t>((theWaveform.size() - 3) / 12);
  TrkData.mainPacket.NumADCPackets = numADCPackets;
  for (size_t i = 0; i < 3; i++) {
    TrkData.mainPacket.SetWaveform(i, theWaveform[i]);
  }
  for (size_t i = 0; i < numADCPackets; i++) {
    TrackerADCPacket adcPacket;
    for (size_t j = 0; j < 12; j++) {
      adcPacket.SetWaveform(j, theWaveform[3 + i * 12 + j]);
    }
    TrkData.adcPacketVec.push_back(adcPacket);
    headerData.TransferByteCount += sizeof(TrackerADCPacket);
    headerData.PacketCount++;
  }
}

ArtBinaryPacketsFromDigis::ArtBinaryPacketsFromDigis(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _generateTimestampTable(config().generateTimestampTable()),
    _tableFile(config().tableFile()), _timestampOffset(config().timestampOffset()),
    _includeTracker(config().includeTracker()), _includeCalorimeter(config().includeCalorimeter()),
    _includeCrv(config().includeCrv()), _includeDMAHeaders(config().includeDMAHeaders()),
    _generateBinaryFile(config().generateBinaryFile()), _outputFile(config().outputFile()),
    _generateTextFile(config().generateTextFile()),
    _sdtoken{consumes<mu2e::StrawDigiCollection>(config().sdtoken())},
    _sdadctoken{consumes<mu2e::StrawDigiADCWaveformCollection>(config().sdtoken())},
    _cdtoken{consumes<mu2e::CaloDigiCollection>(config().cdtoken())},
    _crvtoken{consumes<mu2e::CrvDigiCollection>(config().crvtoken())}, _numWordsWritten(0),
    _numEventsProcessed(0) {

  produces<timestamp>();

  if (_generateBinaryFile == 1) {
    outputStream.open(_outputFile, std::ios::out | std::ios::binary);
  }
}

void ArtBinaryPacketsFromDigis::beginJob() {}

void ArtBinaryPacketsFromDigis::beginRun(art::Run&) {
  mu2e::GeomHandle<mu2e::Calorimeter> ch;
  _calorimeter = ch.get();
}

void ArtBinaryPacketsFromDigis::endJob() {
  if (_generateBinaryFile == 1) {
    outputStream.close();
  }

  if (_generateTimestampTable) {
    std::ofstream tsTableStream;
    tsTableStream.open(_tableFile, std::ios::out | std::ios::binary);
    for (size_t idx = 0; idx < tsTable.size(); idx++) {
      tsTableStream.write(reinterpret_cast<const char*>(&(tsTable[idx].first)), sizeof(timestamp));
      tsTableStream.write(reinterpret_cast<const char*>(&(tsTable[idx].second)), sizeof(timestamp));

      TLOG(TLVL_DEBUG + 3) << "TIMESTAMP_MAPPING: timestamp: " << tsTable[idx].first
                           << " uniqueid: " << tsTable[idx].second;
    }
    tsTableStream << std::flush;
    tsTableStream.close();
  }

  TLOG(TLVL_DEBUG + 0) << "BinaryPacketsFromDataBlocks: "
                       << "Finished writing " << _numWordsWritten << " words from "
                       << _numEventsProcessed << " events to " << _outputFile;
}

void ArtBinaryPacketsFromDigis::produce(art::Event& evt) {

  // unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

  uint64_t eventNum = evt.id().event(); // is not unique! internal counter???//FIXME!
  uint64_t ts = _numEventsProcessed + _timestampOffset;

  TLOG(TLVL_DEBUG + 2) << "ArtBinaryPacketsFromDigis: eventNum: " << eventNum;

  tracker_data_block_list_t trackerData;
  calo_data_block_list_t caloData;
  crv_data_block_list_t crvData;

  if (_includeTracker > 0) {
    processTrackerData(evt, ts, trackerData);
  }

  if (_includeCalorimeter > 0) {
    processCalorimeterData(evt, ts, caloData);
  }

  if (_includeCrv > 0) {
    processCrvData(evt, ts, crvData);
  }

  DTCLib::DTC_Event thisEvent;
  thisEvent.SetEventWindowTag(DTCLib::DTC_EventWindowTag(ts));

  if (_includeTracker > 0) {

    fillTrackerDMABlocks(thisEvent, trackerData);
  }

  if (_includeCalorimeter > 0) {

    fillCalorimeterDMABlocks(thisEvent, caloData);
  }

  if (_includeCrv > 0) {
    fillCrvDMABlocks(thisEvent, crvData);
  }

  // Write all values, including superblock header and DMA header values, to output buffer
  if (_generateBinaryFile == 1) {
    thisEvent.WriteEvent(outputStream);
    outputStream.flush();
  }

  _numEventsProcessed += 1;

  // Store the timestamp and DataBlockCollection in the event
  evt.put(std::unique_ptr<timestamp>(new timestamp(ts)));
  //                evt.put(std::make_unique< raw_data_list_t >(dataStream));

} // end of ::produce

// method....
void ArtBinaryPacketsFromDigis::processCalorimeterData(art::Event& evt, uint64_t& eventNum,
                                                       calo_data_block_list_t& caloDataBlocks) {
  auto const& cdH = evt.getValidHandle(_cdtoken);
  const CaloDigiCollection& hits_CD(*cdH);
  CaloDAQMap const& calodaqconds = _calodaqconds_h.get(evt.id()); // Get calo daq cond

  calo_data_block_list_t tmpCaloDataBlockList;

  for (size_t i = 0; i < hits_CD.size(); ++i) {
    CaloDigi const& CD = hits_CD.at(i);
    // Fill struct with info for current hit
    DataBlockHeader headerData;
    fillCalorimeterHeaderDataPacket(calodaqconds, CD, headerData, eventNum);
    CaloDataPacket caloData;
    fillCalorimeterDataPacket(calodaqconds, CD, caloData);

    printCalorimeterData(caloData);

    tmpCaloDataBlockList.push_back(
        std::pair<DataBlockHeader, CaloDataPacket>(headerData, caloData));
  }

  TLOG(TLVL_DEBUG + 1)
      << "[ArtBinaryPacketsFromDigis::processCalorimeterData ] Total number of calorimeter "
         "non-empty DataBlocks = "
      << tmpCaloDataBlockList.size();

  uint8_t max_dtc_id = number_of_calo_rocs / number_of_calo_rocs_per_dtc;
  if (number_of_calo_rocs % number_of_calo_rocs_per_dtc > 0) {
    max_dtc_id += 1;
  }

  // Loop over the DTC/ROC pairs and generate datablocks for each ROC
  for (uint8_t dtcID = 0; dtcID < max_dtc_id; dtcID++) {

    for (uint8_t rocID = 0; rocID < number_of_calo_rocs_per_dtc; ++rocID) {
      // Find all hits for this event coming from the specified DTC/ROC combination
      bool is_first(true);
      for (size_t curHitIdx = 0; curHitIdx < tmpCaloDataBlockList.size(); curHitIdx++) {
        if (tmpCaloDataBlockList[curHitIdx].first.DTCID == dtcID &&
            tmpCaloDataBlockList[curHitIdx].first.LinkID == rocID) {

          TLOG(TLVL_DEBUG + 1)
              << "[ArtBinaryPacketsFromDigis::processCalorimeterData ] filling Hit from DTCID = "
              << (int)dtcID << " ROCID = " << (int)rocID;

          if (is_first) {
            is_first = false;
            caloDataBlocks.push_back(tmpCaloDataBlockList[curHitIdx]);
          } else {
            addCaloHitToCaloPacket(caloDataBlocks.back(), tmpCaloDataBlockList[curHitIdx].second);
          }
        }
      }

      if (is_first) {
        // No hits, so just fill a header packet and no data packets
        DataBlockHeader headerData;
        CaloDataPacket caloData;

        fillEmptyHeaderDataPacket(headerData, eventNum, rocID, dtcID,
                                  DTCLib::DTC_Subsystem_Calorimeter);

        caloDataBlocks.push_back(std::pair<DataBlockHeader, CaloDataPacket>(headerData, caloData));
      } else {
        fillHeaderByteAndPacketCounts(caloDataBlocks.back());
      }

    } // Done looping over the ROCs in a given DTC
  }
}

//--------------------------------------------------------------------------------
// Fix header ByteCount and PacketCount fields
//--------------------------------------------------------------------------------

void ArtBinaryPacketsFromDigis::fillHeaderByteAndPacketCounts(calo_data_block_t& caloData) {
  caloData.first.TransferByteCount =
      16 /*header packet*/ + sizeof(uint16_t) /* num hits */ +
      (sizeof(uint16_t) + sizeof(CalorimeterHitDataPacket)) * caloData.second.hitPacketVec.size();

  auto idxPos = sizeof(uint16_t) + sizeof(uint16_t) * caloData.second.hitPacketVec.size();
  for (auto& vec : caloData.second.waveformVec) {
    caloData.first.TransferByteCount += sizeof(adc_t) * vec.size();
    caloData.second.hitIndex.push_back(idxPos);

    idxPos += sizeof(CalorimeterHitDataPacket) + sizeof(adc_t) * vec.size();
  }

  while (caloData.first.TransferByteCount % 16 != 0)
    caloData.first.TransferByteCount++;

  caloData.first.PacketCount = (caloData.first.TransferByteCount - 16) / 16;
}

//--------------------------------------------------------------------------------
// crate a caloPacket from the digi
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterDataPacket(CaloDAQMap const& calodaqconds,
                                                          const CaloDigi& CD,
                                                          CaloDataPacket& CaloData) {
  // CaloData.dataPacket.NumberOfHits = 1;

  // CalorimeterBoardID ccBoardID;

  CaloSiPMId offId(CD.SiPMID());
  //  uint16_t roId      = CD.SiPMID();
  uint16_t crystalId = offId.crystal().id();
  TLOG(TLVL_DEBUG + 1) << "...FromDigis: cryId " << crystalId << " roId " << offId.id();

  CaloRawSiPMId rawId = calodaqconds.rawId(offId);
  uint16_t globalROCID = rawId.dirac();
  uint16_t DiracChannel = rawId.ROCchannel();
  uint16_t DetType = offId.detType();
  uint16_t packetId = globalROCID | (DiracChannel << 8) | (DetType << 13);

  TLOG(TLVL_DEBUG + 1) << "..FromDigis: DTYPE " << DetType << " ROCID  " << globalROCID << " CHAN "
                       << DiracChannel << (DetType == 1 ? " Caphri" : "");

  CaloData.dataFooterPacket.BoardID = globalROCID % number_of_calo_rocs_per_dtc;
  CaloData.dataFooterPacket.ChannelStatusFlagA = 0;
  CaloData.dataFooterPacket.ChannelStatusFlagC = 0;
  CaloData.dataFooterPacket.unused = 0;

  CaloData.dataPacket.ChannelNumber = DiracChannel; // modified as it should be in the packet
  CaloData.dataPacket.DIRACA = packetId;            // Change-5
  CaloData.dataPacket.DIRACB =
      (((CD.SiPMID() % 2) << 12) | (crystalId)); // this is useless for the moment .. can be a test
  CaloData.dataPacket.ErrorFlags = 0;
  CaloData.dataPacket.Time = CD.t0();
  std::vector<adc_t> theWaveform;
  for (size_t i = 0; i < CD.waveform().size(); ++i) {
    theWaveform.push_back((adc_t)CD.waveform().at(i));
  }
  CaloData.dataPacket.NumberOfSamples = theWaveform.size();
  CaloData.dataPacket.IndexOfMaxDigitizerSample = waveformMaximumIndex(theWaveform);
  CaloData.hitPacketVec.push_back(CaloData.dataPacket); // TODO - where from????

  CaloData.waveformVec.push_back(theWaveform);
}

//--------------------------------------------------------------------------------
// add a caloHit to a caloPacketVector
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::addCaloHitToCaloPacket(calo_data_block_t& caloDataBlock,
                                                       CaloDataPacket& caloHit) {
  // caloDataBlock.second.dataPacket.NumberOfHits += 1;

  caloDataBlock.second.hitPacketVec.push_back(caloHit.hitPacketVec[0]); // hitPacket);
  caloDataBlock.second.waveformVec.push_back(caloHit.waveformVec[0]);

  // increase the size of the block in the header
  caloDataBlock.first.TransferByteCount +=
      sizeof(uint16_t) * (caloHit.hitPacketVec[0].NumberOfSamples + 1) +
      sizeof(CalorimeterHitDataPacket);
  caloDataBlock.first.PacketCount = std::ceil((caloDataBlock.first.TransferByteCount - 16) / 16);
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                         calo_data_block_list_t& caloData) {

  bool first = true;
  auto curDTCID = caloData.front().first.DTCID;
  for (size_t dataBlockIdx = 0; dataBlockIdx < caloData.size(); dataBlockIdx++) {

    // Add the current DataBlock to the current SuperBlock
    // curDataBlock.setTimestamp(ts); // Overwrite the timestamp

    fillCalorimeterDataStream(currentEvent, caloData[dataBlockIdx]);

    if (first || curDTCID != caloData[dataBlockIdx].first.DTCID) {
      TLOG(TLVL_DEBUG + 1) << "\t\tDTCID: " << (int)caloData[dataBlockIdx].first.DTCID;
      TLOG(TLVL_DEBUG + 1) << "\t\tSYSID: " << (int)caloData[dataBlockIdx].first.SubsystemID;
      first = false;
      curDTCID = caloData[dataBlockIdx].first.DTCID;
    }
    if (caloData[dataBlockIdx].first.PacketCount > 0) {
      printHeader(caloData[dataBlockIdx].first);

      printCalorimeterData(caloData[dataBlockIdx].second);
    }

  } // End loop over DataBlocks
}

//--------------------------------------------------------------------------------
//  method to fill the datastream with the calorimeter packets
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterDataStream(DTCLib::DTC_Event& currentEvent,
                                                          calo_data_block_t& caloData) {

  size_t sz = sizeof(DataBlockHeader);
  // check that the trkDataBlock is not empty
  if (caloData.second.hitPacketVec.size() != 0) {
    sz += sizeof(CalorimeterHitDataPacket) + caloData.second.hitPacketVec.size() * sizeof(uint16_t);
    for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {
      auto nSamples = caloData.second.hitPacketVec[i].NumberOfSamples;
      sz += sizeof(uint16_t) * nSamples + sizeof(CalorimeterHitDataPacket);
    }
  }
  while (sz % 16 != 0)
    sz++;

  if (sz >= 0x10000) { // Maximum transfer size from driver
    throw cet::exception("Online-RECO")
        << "ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz < sizeof(mu2e_databuff_t)";
  }
  if (sz != caloData.first.TransferByteCount) {
    throw cet::exception("Online-RECO")
        << "ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz == caloData.first.ByteCount";
  }

  uint8_t dtcID = caloData.first.DTCID;
  DTCLib::DTC_DataBlock thisBlock(sz);

  if (thisBlock.blockPointer == nullptr) {
    throw cet::exception("MemoryAllocationError")
        << "Unable to allocate memory for Tracker block! sz=" << sz;
  }

  auto pos = 0;
  memcpy(thisBlock.allocBytes->data(), &caloData.first, sizeof(DataBlockHeader));
  // Copies the values of num bytes from the location pointed to by source directly to the memory
  // block pointed to by destination (destination, source, bytes)

  pos += sizeof(DataBlockHeader);

  if (caloData.second.hitPacketVec.size() != 0) {

    uint16_t hitCount = caloData.second.hitPacketVec.size();
    memcpy(thisBlock.allocBytes->data() + pos, &hitCount, sizeof(uint16_t));
    pos += sizeof(uint16_t);

    memcpy(thisBlock.allocBytes->data() + pos, &caloData.second.hitIndex[0],
           sizeof(uint16_t) * caloData.second.hitIndex.size());
    pos += sizeof(uint16_t) * caloData.second.hitIndex.size();

    for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {

      memcpy(thisBlock.allocBytes->data() + pos, &(caloData.second.hitPacketVec[i]),
             sizeof(CalorimeterHitDataPacket));
      pos += sizeof(CalorimeterHitDataPacket);

      auto waveform_size = sizeof(uint16_t) * (caloData.second.waveformVec[i].size());
      memcpy(thisBlock.allocBytes->data() + pos, &(caloData.second.waveformVec[i][0]),
             waveform_size);
      pos += waveform_size;
    } // end loop over the CalorimeterHitDataPacketVector
  }
    putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_Calorimeter, thisBlock);
}

//--------------------------------------------------------------------------------
// create the header for the caloPacket
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterHeaderDataPacket(CaloDAQMap const& calodaqconds,
                                                                const CaloDigi& CD,
                                                                DataBlockHeader& HeaderData,
                                                                uint64_t& EventNum) {
  bzero(&HeaderData, sizeof(DataBlockHeader));
  // Word 0
  adc_t nBytes =
      sizeof(DataBlockHeader) +
      sizeof(
          CalorimeterHitDataPacket); // this needs to be increased every time a new hit is addeded!
  HeaderData.TransferByteCount = nBytes;
  // Word 1
  HeaderData.PacketType = 5; // PacketType::Dataheader;

  // get only Dirac# and DetType from roID and DMAP ....
  // ---------------------------------------------------------------
  CaloSiPMId offId = CaloSiPMId(CD.SiPMID());
  CaloRawSiPMId rawId = calodaqconds.rawId(offId);
  uint16_t globalROCID = rawId.dirac();
  uint16_t DetType = offId.detType();
  // ----------------------------------------------------------------
  if (DetType == 1)
    TLOG(TLVL_DEBUG + 2) << " CAPHRI !!!";

  HeaderData.LinkID = globalROCID % number_of_calo_rocs_per_dtc; // from ROCID call it now LinkID
  HeaderData.SubsystemID = DTCLib::DTC_Subsystem_Calorimeter;
  HeaderData.Valid = 1;
  // Word 2
  HeaderData.PacketCount = 1; // NEEDS TO BE INCREASED EVERY TIME A NEW HIT IS ADDED!
  // Word 3
  uint64_t timestamp = EventNum;
  HeaderData.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  HeaderData.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  HeaderData.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  HeaderData.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  HeaderData.Version = format_version;
  // Word 7
  HeaderData.DTCID = static_cast<uint8_t>(globalROCID / number_of_calo_rocs_per_dtc);
  uint8_t evbMode = 0; // ask Eric
  HeaderData.EventWindowMode = evbMode;
  TLOG(TLVL_DEBUG + 2) << " >>FromDigi-Header: Dtyp " << DetType << " Dirac# " << globalROCID
                       << " Link-DTC " << HeaderData.LinkID << " DTC " << HeaderData.DTCID;
}

//--------------------------------------------------------------------------------
//  method that process the tracker data
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::processTrackerData(art::Event& evt, uint64_t& eventNum,
                                                   tracker_data_block_list_t& trackerData) {
  auto const& sdH = evt.getValidHandle(_sdtoken);
  const StrawDigiCollection& hits_SD(*sdH);
  auto const& sdadcH = evt.getValidHandle(_sdadctoken);
  const StrawDigiADCWaveformCollection& hits_SDADC(*sdadcH);

  tracker_data_block_list_t tmpTrackerData;

  uint8_t max_dtc_id = number_of_rocs / number_of_rocs_per_dtc;
  if (number_of_rocs % number_of_rocs_per_dtc > 0) {
    max_dtc_id += 1;
  }

  // Loop over the DTC/ROC pairs and generate datablocks for each ROC
  for (uint8_t dtcID = 0; dtcID < max_dtc_id; dtcID++) {

    for (uint8_t rocID = 0; rocID < number_of_rocs_per_dtc; ++rocID) {
      DataBlockHeader headerData;
      fillEmptyHeaderDataPacket(headerData, eventNum, rocID, dtcID, DTCLib::DTC_Subsystem_Tracker);
      std::vector<TrackerFullHitFormat> rocData;
      trackerData.push_back(
          std::pair<DataBlockHeader, std::vector<TrackerFullHitFormat>>(headerData, rocData));

      // Find all hits for this event coming from the specified DTC/ROC combination
      for (size_t curHitIdx = 0; curHitIdx < hits_SD.size(); curHitIdx++) {
        StrawDigi const& SD = hits_SD.at(curHitIdx);
        StrawDigiADCWaveform const& SDADC = hits_SDADC.at(curHitIdx);

        int panel = SD.strawId().getPanel();
        int plane = SD.strawId().getPlane();

        // ROC ID, counting from 0 across all DTCs (for the tracker)
        //    uint8_t localROCID = panel;
        uint8_t globalROCID =
            (plane * 6) + panel; // strawId().uniquePanel() would provide the ROCID
        uint8_t thisDTCID = static_cast<uint8_t>(globalROCID / number_of_rocs_per_dtc);
        if (panel == rocID && thisDTCID == dtcID) {
          trackerData.back().second.emplace_back();
          fillTrackerDataPacket(SD, SDADC, trackerData.back().second.back(),
                                trackerData.back().first);
        }
      }
    } // Done looping over the ROCs in a given DTC
  }
}

//------------------------------------
// Crv Methods
//------------------------------------
void ArtBinaryPacketsFromDigis::processCrvData(art::Event& evt, uint64_t& eventNum,
                                               crv_data_block_list_t& crvDataBlocks) {
  auto const& crvdH = evt.getValidHandle(_crvtoken);
  const CrvDigiCollection& digis(*crvdH);

  auto const& crvChannelMap = _crvChannelMap_h.get(evt.id());

  for (size_t i = 0; i < digis.size(); ++i) {
    CrvDigi const& digi = digis.at(i);

    // Fill struct with info for current hit
    CRVHit hit;
    int rocID;
    fillCrvDataPacket(crvChannelMap, digi, hit, rocID);
    crvDataBlocks[rocID].hits.push_back(hit);
  }

  TLOG(TLVL_DEBUG + 1) << "Total number of CRV digis = " << digis.size();

  // Loop over all ROCs, fill headers for each ROC - even for ROCs without hits
  for (uint8_t rocID = 1; rocID <= number_of_crv_rocs; ++rocID) {
    fillCrvHeaderPacket(crvChannelMap, crvDataBlocks[rocID], rocID,
                        eventNum); // this will create a new entry for ROCs without hits
  }
}

//--------------------------------------------------------------------------------
// crate a crvPacket from the digi
//--------------------------------------------------------------------------------
int16_t ArtBinaryPacketsFromDigis::compressCrvDigi(int16_t adc) {
  // TODO: Temporary implementation until we have the real compression used at the FEBs
  // FEBs use only 12 bits out of the 16 bits
  if (adc > 2047)
    adc = 2047;
  if (adc < -2048)
    adc = -2048;
  return adc;
}

void ArtBinaryPacketsFromDigis::fillCrvDataPacket(const CRVOrdinal& crvChannelMap,
                                                  const CrvDigi& digi, CRVHit& hit, int& rocID) {
  int crvSiPMNumber = digi.GetSiPMNumber();
  uint16_t crvBarIndex = digi.GetScintillatorBarIndex().asUint();
  uint16_t offlineChannel = crvBarIndex * 4 + crvSiPMNumber;

  CRVROC onlineChannel = crvChannelMap.online(offlineChannel);
  rocID = onlineChannel.ROC();
  uint16_t rocPort = onlineChannel.FEB();
  uint16_t febChannel = onlineChannel.FEBchannel();

  hit.first.febChannel = febChannel;
  hit.first.portNumber = rocPort;
  hit.first.controllerNumber = rocID;
  hit.first.HitTime = digi.GetStartTDC();
  hit.first.NumSamples = CrvDigi::NSamples;
  hit.second.resize(CrvDigi::NSamples);
  for (size_t i = 0; i < CrvDigi::NSamples; ++i)
    hit.second.at(i).ADC = compressCrvDigi(digi.GetADCs().at(i));
}

//--------------------------------------------------------------------------------
// create the header for the crvPacket
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvHeaderPacket(const CRVOrdinal& crvChannelMap,
                                                    CrvDataPacket& crvData, uint8_t rocID,
                                                    uint64_t eventNum) {
  size_t nHits = crvData.hits.size();

  //----------------------------------------------
  // DataBlockHeader //TODO: This may have changed
  //----------------------------------------------
  // Word 0
  adc_t nBytes = sizeof(DataBlockHeader) + sizeof(CRVROCStatusPacket) +
                 (sizeof(CRVHitInfo) + sizeof(CRVHitWaveformSample) * CrvDigi::NSamples) * nHits;
  while (nBytes % 16 != 0)
    nBytes++;
  crvData.header.TransferByteCount = nBytes;
  // Word 1
  crvData.header.PacketType = DTCLib::DTC_PacketType_DataHeader;

  crvData.header.LinkID = rocID;
  crvData.header.SubsystemID = DTCLib::DTC_Subsystem_CRV;
  crvData.header.Valid = 1;
  // Word 2
  // That's how pcie_linux_kernel_module/dtcInterfaceLib/DTC.cpp
  // interpretes it, but it seems redundant
  crvData.header.PacketCount = (crvData.header.TransferByteCount - 16) / 16;
  // Word 3
  uint64_t timestamp =
      eventNum; // TODO: seems to be identical to the microbunch number and EventWindowTag
  crvData.header.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  crvData.header.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  crvData.header.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  crvData.header.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  crvData.header.Version = format_version;
  // Word 7
  crvData.header.DTCID = (rocID - 1) / 9; // DTC0: ROCs 1...9, DTC1: ROCs 10...17
  uint8_t evbMode = 0;                    // ask Eric
  crvData.header.EventWindowMode = evbMode;

  //------------------
  // CRVROCStatusPacket
  //------------------
  // Word 0
  crvData.rocStatus.unused1 = 0;
  crvData.rocStatus.PacketType = 0x06;
  crvData.rocStatus.ControllerID = rocID;
  // Word 1
  crvData.rocStatus.ControllerEventWordCount =
      (sizeof(CRVROCStatusPacket) +
       (sizeof(CRVHitInfo) + sizeof(CRVHitWaveformSample) * CrvDigi::NSamples) * nHits) /
      2;
  // Word 2
  crvData.rocStatus.ActiveFEBFlags2 = 0xFF;
  crvData.rocStatus.unused2 = 0;
  // Word 3
  crvData.rocStatus.ActiveFEBFlags0 = 0xFF;
  crvData.rocStatus.ActiveFEBFlags1 = 0xFF;
  // Word 4
  static uint16_t triggerCount = 0;
  crvData.rocStatus.TriggerCount = ++triggerCount; // TODO: This seems to be a running number
  // Word 5
  crvData.rocStatus.MicroBunchStatus = 0x0FFF;
  // Word 6
  crvData.rocStatus.EventWindowTag1 = (eventNum >> 16);
  // Word 7
  crvData.rocStatus.EventWindowTag0 = eventNum;
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                 const crv_data_block_list_t& crvDataBlocks) {
  // Loop over all ROCs
  uint8_t currentDTCID = 0;
  for (uint8_t rocID = 1; rocID <= number_of_crv_rocs; ++rocID) {
    // Add the current DataBlock to the current SuperBlock
    // curDataBlock.setTimestamp(ts); // Overwrite the timestamp
    const CrvDataPacket& crvData = crvDataBlocks.at(rocID);
    fillCrvDataStream(currentEvent, crvData);

    if (rocID == 1 || currentDTCID != crvData.header.DTCID) {
      TLOG(TLVL_DEBUG + 1) << "\t\tDTCID: " << (int)crvData.header.DTCID;
      TLOG(TLVL_DEBUG + 1) << "\t\tSYSID: " << (int)crvData.header.SubsystemID;
      currentDTCID = crvData.header.DTCID;
    }
    if (crvData.header.PacketCount > 0) {
      printHeader(crvData.header);
      printCrvData(crvData);
    }

  } // End loop over DataBlocks
}

//--------------------------------------------------------------------------------
//  method to fill the datastream with the crv packets
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvDataStream(DTCLib::DTC_Event& currentEvent,
                                                  const CrvDataPacket& crvData) {
  size_t sz =
      crvData.header.TransferByteCount; // byte count was increased to get full chunks of 16 bytes

  uint8_t dtcID = crvData.header.DTCID;
  DTCLib::DTC_DataBlock thisBlock(sz);

  if (thisBlock.blockPointer == nullptr) {
    throw cet::exception("MemoryAllocationError")
        << "Unable to allocate memory for CRV block! sz=" << sz;
  }

  auto pos = 0;
  memcpy(thisBlock.allocBytes->data(), &crvData.header, sizeof(DataBlockHeader));
  pos += sizeof(DataBlockHeader);
  memcpy(thisBlock.allocBytes->data() + pos, &crvData.rocStatus, sizeof(CRVROCStatusPacket));
  pos += sizeof(CRVROCStatusPacket);

  uint16_t hitCount = crvData.hits.size();

  for (size_t i = 0; i < hitCount; i++) {
    memcpy(thisBlock.allocBytes->data() + pos, &crvData.hits[i].first, sizeof(CRVHitInfo));
    pos += sizeof(CRVHitInfo);
    memcpy(thisBlock.allocBytes->data() + pos, &crvData.hits[i].second[0],
           sizeof(CRVHitWaveformSample) * CrvDigi::NSamples);
    pos += sizeof(CRVHitWaveformSample) * CrvDigi::NSamples;
  }

    putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_CRV, thisBlock);
  
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ArtBinaryPacketsFromDigis)
