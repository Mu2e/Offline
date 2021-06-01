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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// geometry
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"

// mu2e-artdaq-core includes
#include "mu2e-artdaq-core/Overlays/CRVFragment.hh"
#include "mu2e-artdaq-core/Overlays/CalorimeterFragment.hh"
#include "mu2e-artdaq-core/Overlays/TrackerFragment.hh"

// pci_linux_kernel_module includes
#include "dtcInterfaceLib/DTC_Packets.h"

// Mu2e includes.
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "DAQDataProducts/inc/DataBlockCollection.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "CaloConditions/inc/CaloDAQMap.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>

// Typedefs needed for compatibility with HLS codeblock
using timestamp = uint64_t;

using DataBlockHeader = DataHeaderPacket;
using TrackerDataPacket = mu2e::TrackerFragment::TrackerDataPacket;
using TrackerADCPacket = mu2e::TrackerFragment::TrackerADCPacket;
using adc_t = uint16_t;
using CalorimeterDataPacket = mu2e::CalorimeterFragment::CalorimeterDataPacket;
using CalorimeterBoardID = mu2e::CalorimeterFragment::CalorimeterBoardID;
using CalorimeterHitReadoutPacket = mu2e::CalorimeterFragment::CalorimeterHitReadoutPacket;
using CRVROCStatusPacket = mu2e::CRVFragment::CRVROCStatusPacket;
using CRVHitReadoutPacket = mu2e::CRVFragment::CRVHitReadoutPacket;

// data struct for the calorimeter
struct CaloDataPacket {
  CalorimeterDataPacket dataPacket;
  CalorimeterBoardID boardID;
  std::vector<CalorimeterHitReadoutPacket> hitPacketVec;
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
  std::vector<CRVHitReadoutPacket> hits;

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
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostic Level")};
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
  int    _includeTracker;
  int    _includeCalorimeter;
  int    _includeCrv;
  int    _includeDMAHeaders;

  // -- include proditions handling 
  ProditionsHandle<CaloDAQMap> _calodaqconds_h;
  // Set to 1 to save packet data to a binary file
  int    _generateBinaryFile;

  std::string   _outputFile;
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
  const size_t number_of_calo_rocs         = 136;
  const size_t number_of_crystals_per_roc  = 10;
  const size_t number_of_calo_rocs_per_dtc = 6;

  //--------------------------------------------------------------------------------
  // CRV ROC/DTC INFO
  //--------------------------------------------------------------------------------

  const size_t number_of_crv_rocs = 16;
  const size_t number_of_crv_rocs_per_dtc = 8;

  //--------------------------------------------------------------------------------

  int _generateTextFile;

  // Diagnostics level.
  int _diagLevel;

  // Limit on number of events for which there will be full printout.
  int _maxFullPrint;

  // Label of the module that made the hits.
  art::ProductToken<StrawDigiCollection> const _sdtoken;
  art::ProductToken<StrawDigiADCWaveformCollection> const _sdadctoken;
  art::ProductToken<CaloDigiCollection> const _cdtoken;
  art::ProductToken<CrvDigiCollection> const _crvtoken;

  size_t _numWordsWritten;
  size_t _numEventsProcessed;

  const Calorimeter* _calorimeter; // cached pointer to the calorimeter geometry
  const CosmicRayShield* _crv;     // cached pointer to the crv geometry

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
  void   fillCalorimeterDataPacket(CaloDAQMap const& calodaqconds,
				   const CaloDigi& SD, CaloDataPacket& caloData);
  
  void addCaloHitToCaloPacket(calo_data_block_t& dataBlock, CaloDataPacket& caloData);

  void   fillCalorimeterHeaderDataPacket(CaloDAQMap const& calodaqconds,
					 const CaloDigi& SD, DataBlockHeader& HeaderData, uint64_t& EventNum);
 
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
  uint8_t compressCrvDigi(int adc);
  void fillCrvDataPacket(const CrvDigi& digi, CRVHitReadoutPacket& hit, int& globalRocID);
  void fillCrvHeaderPacket(CrvDataPacket& crvData, uint8_t globalRocID, uint64_t eventNum);
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
  printf("[ArtBinaryPacketsFromDigis::printHeader] START header print  \n");
  printf("[ArtBinaryPacketsFromDigis::printHeader] ByteCount      : %i \n",
         headerDataBlock.s.TransferByteCount);
  printf("[ArtBinaryPacketsFromDigis::printHeader] Resv1       : %i \n", headerDataBlock.s.Resv1);
  printf("[ArtBinaryPacketsFromDigis::printHeader] PacketType     : %i \n",
         headerDataBlock.s.PacketType);
  printf("[ArtBinaryPacketsFromDigis::printHeader] ROCID 	   : %i \n", headerDataBlock.s.LinkID);
  printf("[ArtBinaryPacketsFromDigis::printHeader] SubsystemID    : %i \n",
         headerDataBlock.s.SubsystemID);
  printf("[ArtBinaryPacketsFromDigis::printHeader] Valid 	   : %i \n", headerDataBlock.s.Valid);
  printf("[ArtBinaryPacketsFromDigis::printHeader] PacketCount    : %i \n",
         headerDataBlock.s.PacketCount);
  printf("[ArtBinaryPacketsFromDigis::printHeader] Resv2        : %i \n", headerDataBlock.s.Resv2);
  printf("[ArtBinaryPacketsFromDigis::printHeader] ts10   : %i \n", headerDataBlock.s.ts10);
  printf("[ArtBinaryPacketsFromDigis::printHeader] ts32   : %i \n", headerDataBlock.s.ts32);
  printf("[ArtBinaryPacketsFromDigis::printHeader] ts54  : %i \n", headerDataBlock.s.ts54);
  printf("[ArtBinaryPacketsFromDigis::printHeader] Status	   : %i \n", headerDataBlock.s.Status);
  printf("[ArtBinaryPacketsFromDigis::printHeader] FormatVersion  : %i \n",
         headerDataBlock.s.Version);
  printf("[ArtBinaryPacketsFromDigis::printHeader] DTCID	   : %i \n", headerDataBlock.s.DTCID);
  printf("[ArtBinaryPacketsFromDigis::printHeader] EVBMode        : %i \n",
         headerDataBlock.s.EventWindowMode);
}

void ArtBinaryPacketsFromDigis::printTrackerData(std::vector<TrackerFullHitFormat> const& trkData) {
  printf("[ArtBinaryPacketsFromDigis::printTrackerData] START tracker-data print \n");
  for (size_t i = 0; i < trkData.size(); i++) {
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] StrawIndex    : %i \n",
           (int)trkData[i].mainPacket.StrawIndex);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC0		: %i \n",
           (int)trkData[i].mainPacket.TDC0());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC1		: %i \n",
           (int)trkData[i].mainPacket.TDC1());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT0		: %i \n",
           (int)trkData[i].mainPacket.TOT0);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT1		: %i \n",
           (int)trkData[i].mainPacket.TOT1);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] PMP             : %i \n",
           (int)trkData[i].mainPacket.PMP);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC00         : %i \n",
           (int)trkData[i].mainPacket.ADC00);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC01  	: %i \n",
           (int)trkData[i].mainPacket.ADC01());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC02  	: %i \n",
           (int)trkData[i].mainPacket.ADC02);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ErrorFlags : %i \n",
           (int)trkData[i].mainPacket.ErrorFlags);
  }
}

void ArtBinaryPacketsFromDigis::printCalorimeterData(CaloDataPacket const& caloData) {
  CalorimeterDataPacket packet = caloData.dataPacket;
  CalorimeterBoardID boardId = caloData.boardID;
  size_t nHits = caloData.hitPacketVec.size();
  printf("[ArtBinaryPacketsFromDigis::printCaloData] START calorimeter-data print \n");
  printf("[ArtBinaryPacketsFromDigis::printCaloData] NumberofHits        : %i \n",
         (int)packet.NumberOfHits);
  printf("[ArtBinaryPacketsFromDigis::printCaloData] BoardID             : %i \n",
         (int)boardId.BoardID);
  printf("[ArtBinaryPacketsFromDigis::printCaloData] ChannelStatusFlagsA : %i \n",
         (int)boardId.ChannelStatusFlagsA);
  printf("[ArtBinaryPacketsFromDigis::printCaloData] ChannelStatusFlagsB : %i \n",
         (int)boardId.ChannelStatusFlagsB);
  printf("[ArtBinaryPacketsFromDigis::printCaloData] unused              : %i \n",
         (int)boardId.unused);
  printf("[ArtBinaryPacketsFromDigis::printCaloData] NHits               : %i \n", (int)nHits);

  for (size_t i = 0; i < nHits; ++i) {
    CalorimeterHitReadoutPacket const& hit = caloData.hitPacketVec[i];
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t hit : %i \n", (int)i);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t ChannelNumber : %i \n",
           (int)hit.ChannelNumber);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t DIRACA        : %i \n", (int)hit.DIRACA);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t DIRACB        : %i \n", (int)hit.DIRACB);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t ErrorFlags    : %i \n", (int)hit.ErrorFlags);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t Time          : %i \n", (int)hit.Time);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t NumberOfSamples : %i \n",
           (int)hit.NumberOfSamples);
    printf("[ArtBinaryPacketsFromDigis::printCaloData]\t IndexOfMaxDigitizerSample : %i \n",
           (int)hit.IndexOfMaxDigitizerSample);
  }
}

void ArtBinaryPacketsFromDigis::printCrvData(CrvDataPacket const& crvData) {
  size_t nHits = crvData.hits.size();
  printf("[ArtBinaryPacketsFromDigis::printCrvData] START crv-data print \n");
  printf("[ArtBinaryPacketsFromDigis::printCrvData] ROC controller ID   : %i \n",
         (int)crvData.rocStatus.ControllerID);
  printf("[ArtBinaryPacketsFromDigis::printCrvData] Errors              : %i \n",
         (int)crvData.rocStatus.Errors);
  printf("[ArtBinaryPacketsFromDigis::printCrvData] NHits               : %i \n", (int)nHits);

  for (size_t i = 0; i < nHits; ++i) {
    printf("[ArtBinaryPacketsFromDigis::printCrvData] hit : %i \n", (int)i);
    printf("[ArtBinaryPacketsFromDigis::printCrvData] Channel       : %i \n",
           (int)(crvData.hits[i].SiPMID & 0x7F));
    printf("[ArtBinaryPacketsFromDigis::printCrvData] FEB           : %i \n",
           (int)(crvData.hits[i].SiPMID >> 7));
    printf("[ArtBinaryPacketsFromDigis::printCrvData] Time          : %i \n",
           (int)crvData.hits[i].HitTime);
    printf("[ArtBinaryPacketsFromDigis::printCrvData] NumOfSamples  : %i \n",
           (int)crvData.hits[i].NumSamples);
  }
}

void ArtBinaryPacketsFromDigis::fillTrackerDataStream(DTCLib::DTC_Event& currentEvent,
                                                      tracker_data_block_t const& trackerData) {

  auto sz = sizeof(DataBlockHeader);
  // check that the trkDataBlock is not empty

  if (trackerData.first.s.PacketCount > 0) {
    sz += sizeof(TrackerDataPacket) * trackerData.first.s.PacketCount;
  }

  uint8_t dtcID = trackerData.first.s.DTCID;
  DTCLib::DTC_DataBlock thisBlock(sz);

  if (thisBlock.blockPointer == nullptr) {
    throw cet::exception("MemoryAllocationError")
        << "Unable to allocate memory for Tracker block! sz=" << sz;
  }

  auto pos = 0;
  memcpy(thisBlock.allocBytes->data(), &trackerData.first, sizeof(DataBlockHeader));
  pos += sizeof(DataBlockHeader);

  if (trackerData.first.s.PacketCount > 0) {
    if (sizeof(TrackerDataPacket) % 16 !=
        0) { // Make sure that TrackerDataPacket is an even number of DataPackets!
      throw cet::exception("Online-RECO") << "ArtBinaryPacketsFromDigis::fillTrackerDataStream : "
                                             "sizeof(TrackerDataPacket) % 16 == 0"
                                          << std::endl;
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
    putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_Tracker, thisBlock);
  }
}

void ArtBinaryPacketsFromDigis::fillTrackerDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                     tracker_data_block_list_t const& trkData) {

  auto curDTCID = trkData.front().first.s.DTCID;
  bool first = true;
  if (_diagLevel > 1) {
    std::cout << "[ArtBinaryPacketsFromDigis::fillTrackerDMABlocks] trkData.size() = "
              << trkData.size() << std::endl;
  }
  for (auto& dataBlock : trkData) {

    fillTrackerDataStream(currentEvent, dataBlock);

    if (_diagLevel > 1) {
      if (dataBlock.first.s.DTCID != curDTCID || first) {
        std::cout << "================================================" << std::endl;
        // std::cout << "\t\tTimestamp: " << ts << std::endl;
        std::cout << "\t\tDTCID: " << (int)dataBlock.first.s.DTCID << std::endl;
        std::cout << "\t\tSYSID: " << (int)dataBlock.first.s.SubsystemID << std::endl;
        curDTCID = dataBlock.first.s.DTCID;
        first = false;
      }
      if (dataBlock.first.s.PacketCount > 0) {
        printHeader(dataBlock.first);
        if (_diagLevel > 2) {
          printTrackerData(dataBlock.second);
        }
      }
    }

  } // End loop over DataBlocks
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillEmptyHeaderDataPacket(DataBlockHeader& headerData,
                                                          uint64_t& EventNum, uint8_t& ROCId,
                                                          uint8_t& DTCId, uint8_t Subsys) {

  bzero(&headerData.s, sizeof(DataBlockHeader));
  // Fill in the byte count field of the header packet
  // Word 0
  headerData.s.TransferByteCount = sizeof(DataBlockHeader);
  // Word 1
  headerData.s.Resv1 = 0;      // ask Eric!!!//FIX ME!
  headerData.s.PacketType = 5; // PacketType::Dataheader;

  headerData.s.LinkID = ROCId;

  headerData.s.SubsystemID = Subsys; // DTCLib::DTC_Subsystem_Tracker; //: 3;

  headerData.s.Valid = 1;
  // Word 2
  headerData.s.PacketCount = 0;
  headerData.s.Resv2 = 0; // : 5;
  // Word 3
  uint64_t timestamp = EventNum;
  headerData.s.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  headerData.s.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  headerData.s.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  headerData.s.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  headerData.s.Version = format_version;
  // Word 7
  headerData.s.DTCID = DTCId;
  uint8_t evbMode = 0; // maybe off-spill vs on-spill?
  headerData.s.EventWindowMode = evbMode;
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
  TrkData.mainPacket.EWMCounter = headerData.s.ts10 & 0xF;
  TrkData.mainPacket.PMP = SD.PMP();
  TrkData.mainPacket.ErrorFlags = 0; // FIXME
  TrkData.mainPacket.unused1 = 0;

  headerData.s.TransferByteCount += sizeof(TrackerDataPacket);
  headerData.s.PacketCount++;

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
    headerData.s.TransferByteCount += sizeof(TrackerADCPacket);
    headerData.s.PacketCount++;
  }
}

ArtBinaryPacketsFromDigis::ArtBinaryPacketsFromDigis(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _generateTimestampTable(config().generateTimestampTable()),
    _tableFile(config().tableFile()), _timestampOffset(config().timestampOffset()),
    _includeTracker(config().includeTracker()), _includeCalorimeter(config().includeCalorimeter()),
    _includeCrv(config().includeCrv()), _includeDMAHeaders(config().includeDMAHeaders()),
    _generateBinaryFile(config().generateBinaryFile()), _outputFile(config().outputFile()),
    _generateTextFile(config().generateTextFile()), _diagLevel(config().diagLevel()),
    _maxFullPrint(config().maxFullPrint()), _sdtoken{consumes<mu2e::StrawDigiCollection>(
                                                config().sdtoken())},
    _sdadctoken{consumes<mu2e::StrawDigiADCWaveformCollection>(config().sdtoken())},
    _cdtoken{consumes<mu2e::CaloDigiCollection>(config().cdtoken())},
    _crvtoken{consumes<mu2e::CrvDigiCollection>(config().crvtoken())}, _numWordsWritten(0),
    _numEventsProcessed(0) {

  produces<timestamp>();

  if (_generateBinaryFile == 1) {
    outputStream.open(_outputFile, std::ios::out | std::ios::binary);
  }
}

void ArtBinaryPacketsFromDigis::beginJob() {

  if (_diagLevel > 0) {
    std::cout << "ArtBinaryPacketsFromDigis Diaglevel: " << _diagLevel << " " << _maxFullPrint
              << std::endl;
  }
}

void ArtBinaryPacketsFromDigis::beginRun(art::Run&) {
  mu2e::GeomHandle<mu2e::Calorimeter> ch;
  _calorimeter = ch.get();

  mu2e::GeomHandle<mu2e::CosmicRayShield> crvHandle;
  _crv = crvHandle.get();
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

      if (_diagLevel > 3) {
        std::cout << "TIMESTAMP_MAPPING: timestamp: " << tsTable[idx].first
                  << " uniqueid: " << tsTable[idx].second << std::endl;
      }
    }
    tsTableStream << std::flush;
    tsTableStream.close();
  }

  if (_diagLevel > 0) {
    std::cout << "BinaryPacketsFromDataBlocks: "
              << "Finished writing " << _numWordsWritten << " words from " << _numEventsProcessed
              << " events to " << _outputFile << std::endl;
  }
}

void ArtBinaryPacketsFromDigis::produce(art::Event& evt) {

  // unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

  uint64_t eventNum = evt.id().event(); // is not unique! internal counter???//FIXME!
  uint64_t ts = _numEventsProcessed + _timestampOffset;

  if (_diagLevel > 2) {
    std::cout << "ArtBinaryPacketsFromDigis: eventNum: " << eventNum << std::endl;
  }

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
  //		evt.put(std::make_unique< raw_data_list_t >(dataStream));

} // end of ::produce

// method....
void ArtBinaryPacketsFromDigis::processCalorimeterData(art::Event& evt, uint64_t& eventNum,
                                                       calo_data_block_list_t& caloDataBlocks) {
  auto const& cdH = evt.getValidHandle(_cdtoken);
  const CaloDigiCollection& hits_CD(*cdH);
  CaloDAQMap const& calodaqconds = _calodaqconds_h.get(evt.id()); // Get calo daq cond

  calo_data_block_list_t tmpCaloDataBlockList;

  for (size_t i = 0; i < hits_CD.size(); ++i) {  
    CaloDigi   const& CD = hits_CD.at(i);
    // Fill struct with info for current hit
    DataBlockHeader headerData;
    fillCalorimeterHeaderDataPacket(calodaqconds, CD, headerData, eventNum);
    CaloDataPacket caloData;
    fillCalorimeterDataPacket(calodaqconds, CD, caloData);
    if (_diagLevel > 1) {
      printCalorimeterData(caloData);
    }
    tmpCaloDataBlockList.push_back(
        std::pair<DataBlockHeader, CaloDataPacket>(headerData, caloData));
  }

  if (_diagLevel > 1) {
    std::cout << "[ArtBinaryPacketsFromDigis::processCalorimeterData ] Total number of calorimeter "
                 "non-empty DataBlocks = "
              << tmpCaloDataBlockList.size() << std::endl;
  }

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
        if (tmpCaloDataBlockList[curHitIdx].first.s.DTCID == dtcID &&
            tmpCaloDataBlockList[curHitIdx].first.s.LinkID == rocID) {

	  if (_diagLevel > 1) {
	    std::cout << "[ArtBinaryPacketsFromDigis::processCalorimeterData ] filling Hit from DTCID = "<< (int)dtcID 
		      << " ROCID = " << (int)rocID
		      << std::endl;
	  }
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
  caloData.first.s.TransferByteCount = 16 /*header packet*/ + sizeof(uint16_t) /* num hits */ +
                                       sizeof(CalorimeterBoardID) +
                                       (sizeof(uint16_t) + sizeof(CalorimeterHitReadoutPacket)) *
                                           caloData.second.hitPacketVec.size();

  auto idxPos = sizeof(uint16_t) + sizeof(CalorimeterBoardID) +
                sizeof(uint16_t) * caloData.second.hitPacketVec.size();
  for (auto& vec : caloData.second.waveformVec) {
    caloData.first.s.TransferByteCount += sizeof(adc_t) * vec.size();
    caloData.second.hitIndex.push_back(idxPos);

    idxPos += sizeof(CalorimeterHitReadoutPacket) + sizeof(adc_t) * vec.size();
  }

  while (caloData.first.s.TransferByteCount % 16 != 0)
    caloData.first.s.TransferByteCount++;

  caloData.first.s.PacketCount = (caloData.first.s.TransferByteCount - 16) / 16;
}

//--------------------------------------------------------------------------------
// crate a caloPacket from the digi
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterDataPacket(CaloDAQMap const& calodaqconds,
							  const CaloDigi& CD,
                                                          CaloDataPacket& CaloData) {
  CaloData.dataPacket.NumberOfHits = 1;

  CalorimeterBoardID ccBoardID;
  // Change # 1: get roid and cryID from Digi
  //=========================================
  
  uint16_t roId      = CD.SiPMID();
  uint16_t crystalId = _calorimeter->caloIDMapper().crystalIDFromSiPMID(roId);
  if( _diagLevel==1) printf( "...FromDigis: cryId %d roId %d \n",crystalId,roId);

  //=========================================================================================
  // Change # 2: get packetId from DMAP and roId to extract: Dirac#, Chan# and Dettype 
  // For the moment (wait OTSDAQ) define BoardId as consecutive with Dirac# and 6 Diracs=1DTC
  //=========================================================================================
  uint16_t packetId     = calodaqconds.caloRoIdToPacketId(roId);
  uint16_t globalROCID  = (packetId & (0x00FF));
  uint16_t DiracChannel = (packetId & (0x1F00)) >> 8;
  uint16_t DetType      = (packetId & (0xE000)) >> 13;
  
  if( _diagLevel==1) printf( "..FromDigis: DTYPE %d ROCID %d CHAN %d \n",DetType,globalROCID,DiracChannel);
  if( _diagLevel==1 && DetType==1) printf("Caphri \n");
    
  ccBoardID.BoardID             = globalROCID % number_of_calo_rocs_per_dtc;
  ccBoardID.ChannelStatusFlagsA = 0;
  ccBoardID.ChannelStatusFlagsB = 0;
  ccBoardID.unused              = 0;

  CaloData.boardID = ccBoardID;

  CalorimeterHitReadoutPacket hitPacket;
  hitPacket.ChannelNumber = DiracChannel;  // modified as it should be in the packet
  hitPacket.DIRACA        = packetId; //Change-5
  hitPacket.DIRACB        = (((CD.SiPMID() % 2) << 12) | (crystalId)); // this is useless for the moment .. can be a test
  hitPacket.ErrorFlags    = 0;
  hitPacket.Time          = CD.t0();
  std::vector<adc_t> theWaveform;
  for (size_t i = 0; i < CD.waveform().size(); ++i) {
    theWaveform.push_back((adc_t)CD.waveform().at(i));
  }
  hitPacket.NumberOfSamples           = theWaveform.size();
  hitPacket.IndexOfMaxDigitizerSample = waveformMaximumIndex(theWaveform);
  CaloData.hitPacketVec.push_back(hitPacket);

  CaloData.waveformVec.push_back(theWaveform);
}

//--------------------------------------------------------------------------------
// add a caloHit to a caloPacketVector
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::addCaloHitToCaloPacket(calo_data_block_t& caloDataBlock,
                                                       CaloDataPacket& caloHit) {
  caloDataBlock.second.dataPacket.NumberOfHits += 1;

  caloDataBlock.second.hitPacketVec.push_back(caloHit.hitPacketVec[0]); // hitPacket);
  caloDataBlock.second.waveformVec.push_back(caloHit.waveformVec[0]);

  // increase the size of the block in the header
  caloDataBlock.first.s.TransferByteCount +=
      sizeof(uint16_t) * (caloHit.hitPacketVec[0].NumberOfSamples + 1) +
      sizeof(CalorimeterHitReadoutPacket);
  caloDataBlock.first.s.PacketCount =
      std::ceil((caloDataBlock.first.s.TransferByteCount - 16) / 16);
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                         calo_data_block_list_t& caloData) {

  bool first = true;
  auto curDTCID = caloData.front().first.s.DTCID;
  for (size_t dataBlockIdx = 0; dataBlockIdx < caloData.size(); dataBlockIdx++) {

    // Add the current DataBlock to the current SuperBlock
    // curDataBlock.setTimestamp(ts); // Overwrite the timestamp

    fillCalorimeterDataStream(currentEvent, caloData[dataBlockIdx]);

    if (_diagLevel > 1) {
      if (first || curDTCID != caloData[dataBlockIdx].first.s.DTCID) {
        std::cout << "================================================" << std::endl;
        // std::cout << "\t\tTimestamp: " << ts << std::endl;
        std::cout << "\t\tDTCID: " << (int)caloData[dataBlockIdx].first.s.DTCID << std::endl;
        std::cout << "\t\tSYSID: " << (int)caloData[dataBlockIdx].first.s.SubsystemID << std::endl;
        first = false;
        curDTCID = caloData[dataBlockIdx].first.s.DTCID;
      }
      if (caloData[dataBlockIdx].first.s.PacketCount > 0) {
        printHeader(caloData[dataBlockIdx].first);
        if (_diagLevel > 2) {

          printCalorimeterData(caloData[dataBlockIdx].second);
        }
      }
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
    sz += sizeof(CalorimeterDataPacket) + caloData.second.hitPacketVec.size() * sizeof(uint16_t) +
          sizeof(CalorimeterBoardID);
    for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {
      auto nSamples = caloData.second.hitPacketVec[i].NumberOfSamples;
      sz += sizeof(uint16_t) * nSamples + sizeof(CalorimeterHitReadoutPacket);
    }
  }
  while (sz % 16 != 0)
    sz++;

  if (sz >= sizeof(mu2e_databuff_t)) {
    throw cet::exception("Online-RECO")
        << "ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz < sizeof(mu2e_databuff_t)"
        << std::endl;
  }
  if (sz != caloData.first.s.TransferByteCount) {
    throw cet::exception("Online-RECO")
        << "ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz == caloData.first.ByteCount"
        << std::endl;
  }

  uint8_t dtcID = caloData.first.s.DTCID;
  DTCLib::DTC_DataBlock thisBlock(sz);

  if (thisBlock.blockPointer == nullptr) {
    throw cet::exception("MemoryAllocationError")
        << "Unable to allocate memory for Tracker block! sz=" << sz;
  }

  auto pos = 0;
  memcpy(thisBlock.allocBytes->data(), &caloData.first, sizeof(DataBlockHeader));
  pos += sizeof(DataBlockHeader);

  if (caloData.second.hitPacketVec.size() != 0) {

    uint16_t hitCount = caloData.second.hitPacketVec.size();
    memcpy(thisBlock.allocBytes->data() + pos, &hitCount, sizeof(uint16_t));
    pos += sizeof(uint16_t);

    memcpy(thisBlock.allocBytes->data() + pos, &caloData.second.hitIndex[0],
           sizeof(uint16_t) * caloData.second.hitIndex.size());
    pos += sizeof(uint16_t) * caloData.second.hitIndex.size();

    memcpy(thisBlock.allocBytes->data() + pos, &(caloData.second.boardID),
           sizeof(CalorimeterBoardID));
    pos += sizeof(CalorimeterBoardID);

    for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {

      memcpy(thisBlock.allocBytes->data() + pos, &(caloData.second.hitPacketVec[i]),
             sizeof(CalorimeterHitReadoutPacket));
      pos += sizeof(CalorimeterHitReadoutPacket);

      auto waveform_size = sizeof(uint16_t) * (caloData.second.waveformVec[i].size());
      memcpy(thisBlock.allocBytes->data() + pos, &(caloData.second.waveformVec[i][0]),
             waveform_size);
      pos += waveform_size;
    } // end loop over the calorimeterHitReadoutPacketVector
    putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_Calorimeter, thisBlock);
  }
}

//--------------------------------------------------------------------------------
// create the header for the caloPacket
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCalorimeterHeaderDataPacket(CaloDAQMap const& calodaqconds,
								const CaloDigi& CD,
                                                                DataBlockHeader& HeaderData,
                                                                uint64_t& EventNum) {
  bzero(&HeaderData.s, sizeof(DataBlockHeader));
  // Word 0
  adc_t nBytes =
      sizeof(DataBlockHeader) + sizeof(CalorimeterDataPacket) +
      sizeof(CalorimeterBoardID); // this needs to be increased every time a new hit is addeded!
  HeaderData.s.TransferByteCount = nBytes;
  // Word 1
  HeaderData.s.PacketType = 5; // PacketType::Dataheader;

  // get only Dirac# and DetType from roID and DMAP ....
  // ---------------------------------------------------------------
  size_t roId = CD.SiPMID();
  //    size_t globalROCID = crystalId / number_of_crystals_per_roc;
  uint16_t packetId     = calodaqconds.caloRoIdToPacketId(roId);
  uint16_t globalROCID  = (packetId & (0x00FF));
  uint16_t DetType      = (packetId & (0xE000)) >> 13;
  // ----------------------------------------------------------------
  if( _diagLevel==1 && DetType == 1) printf(" CAPHRI !!! \n");

  HeaderData.s.LinkID      = globalROCID % number_of_calo_rocs_per_dtc;// from ROCID call it now LinkID
  HeaderData.s.SubsystemID = DTCLib::DTC_Subsystem_Calorimeter;
  HeaderData.s.Valid       = 1;
  // Word 2
  HeaderData.s.PacketCount = 1; // NEEDS TO BE INCREASED EVERY TIME A NEW HIT IS ADDED!
  // Word 3
  uint64_t timestamp = EventNum;
  HeaderData.s.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  HeaderData.s.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  HeaderData.s.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  HeaderData.s.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  HeaderData.s.Version = format_version;
  // Word 7
  HeaderData.s.DTCID = static_cast<uint8_t>(globalROCID / number_of_calo_rocs_per_dtc);
  uint8_t evbMode = 0; // ask Eric
  HeaderData.s.EventWindowMode = evbMode;
  if( _diagLevel==1) printf(" >>FromDigi-Header: Dtyp Dirac# Link-DTC DTC %d %d %d %d \n",
			    DetType,globalROCID,HeaderData.s.LinkID,HeaderData.s.DTCID);
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
        uint8_t globalROCID = (plane * 6) + panel;// strawId().uniquePanel() would provide the ROCID
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

  for (size_t i = 0; i < digis.size(); ++i) {
    CrvDigi const& digi = digis.at(i);

    // Fill struct with info for current hit
    CRVHitReadoutPacket hit;
    int globalRocID;
    fillCrvDataPacket(digi, hit, globalRocID);
    crvDataBlocks[globalRocID].hits.push_back(hit);
  }

  if (_diagLevel > 1) {
    std::cout << "[ArtBinaryPacketsFromDigis::processCrvData ] Total number of CRV digis = "
              << digis.size() << std::endl;
  }

  // Loop over all ROCs, fill headers for each ROC - even for ROCs without hits
  for (uint8_t globalRocID = 0; globalRocID < number_of_crv_rocs; globalRocID++) {
    fillCrvHeaderPacket(crvDataBlocks[globalRocID], globalRocID,
                        eventNum); // this will create a new entry for ROCs without hits
  }
}

//--------------------------------------------------------------------------------
// crate a crvPacket from the digi
//--------------------------------------------------------------------------------
uint8_t ArtBinaryPacketsFromDigis::compressCrvDigi(int adc) {
  // TODO: Temporary implementation until we have the real compression used at the FEBs
  adc -= 95;
  if (adc < 0)
    adc = 0;
  uint8_t toReturn = adc;
  if (adc > 50 && adc <= 100)
    toReturn = 50 + (adc - 50) / 2;
  if (adc > 100 && adc <= 200)
    toReturn = 75 + (adc - 100) / 4;
  if (adc > 200 && adc <= 400)
    toReturn = 100 + (adc - 200) / 8;
  if (adc > 400 && adc <= 2480)
    toReturn = 125 + (adc - 400) / 16;
  if (adc > 2480)
    toReturn = 255;
  return toReturn;
}

void ArtBinaryPacketsFromDigis::fillCrvDataPacket(const CrvDigi& digi, CRVHitReadoutPacket& hit,
                                                  int& globalRocID) {
  // TODO: This is a temporary implementation.
  // There will be a major change on the barIndex+SiPMNumber system,
  // which will be replaced by a channel ID system
  int crvSiPMNumber = digi.GetSiPMNumber();
  mu2e::CRSScintillatorBarIndex crvBarIndex = digi.GetScintillatorBarIndex();
  // Only a toy model is used here. The real implementation will follow.
  int channel = (crvBarIndex.asUint() * 4 + crvSiPMNumber) % 64; // channel within an FEB
  int FEB = (crvBarIndex.asUint() * 4 + crvSiPMNumber) / 64;     // globale FEBId
  uint16_t SiPMID = (FEB << 7) | channel;
  globalRocID = FEB / 24; // global ROCId

  hit.SiPMID = SiPMID;
  hit.HitTime = digi.GetStartTDC();
  hit.NumSamples = 8;
  hit.WaveformSample0 = compressCrvDigi(
      digi.GetADCs().at(0)); // TODO: There should be a better way of filling the waveform
  hit.WaveformSample1 = compressCrvDigi(digi.GetADCs().at(1));
  hit.WaveformSample2 = compressCrvDigi(digi.GetADCs().at(2));
  hit.WaveformSample3 = compressCrvDigi(digi.GetADCs().at(3));
  hit.WaveformSample4 = compressCrvDigi(digi.GetADCs().at(4));
  hit.WaveformSample5 = compressCrvDigi(digi.GetADCs().at(5));
  hit.WaveformSample6 = compressCrvDigi(digi.GetADCs().at(6));
  hit.WaveformSample7 = compressCrvDigi(digi.GetADCs().at(7));
}

//--------------------------------------------------------------------------------
// create the header for the crvPacket
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvHeaderPacket(CrvDataPacket& crvData, uint8_t globalRocID,
                                                    uint64_t eventNum) {
  size_t nHits = crvData.hits.size();

  //--------------
  // DataBlocHeader
  //--------------
  // Word 0
  adc_t nBytes =
      sizeof(DataBlockHeader) + sizeof(CRVROCStatusPacket) + sizeof(CRVHitReadoutPacket) * nHits;
  while (nBytes % 16 != 0)
    nBytes++;
  crvData.header.s.TransferByteCount = nBytes;
  // Word 1
  crvData.header.s.PacketType = 5; // PacketType::Dataheader;

  crvData.header.s.LinkID = globalRocID % number_of_crv_rocs_per_dtc; // TODO: Is this correct?
  crvData.header.s.SubsystemID = DTCLib::DTC_Subsystem_CRV;
  crvData.header.s.Valid = 1;
  // Word 2
  crvData.header.s.PacketCount =
      (crvData.header.s.TransferByteCount - 16) /
      16; // TODO: That's how pcie_linux_kernel_module/dtcInterfaceLib/DTC.cpp
          // interpretes it, but it seems redundant

  // Word 3
  uint64_t timestamp = eventNum; // TODO: Is this correct?
  crvData.header.s.ts10 = static_cast<adc_t>(timestamp & 0xFFFF);
  // Word 4
  crvData.header.s.ts32 = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
  // Word 5
  crvData.header.s.ts54 = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
  // Word 6
  crvData.header.s.Status = 0; // 0 corresponds to "TimeStamp had valid data"
  crvData.header.s.Version = format_version;
  // Word 7
  crvData.header.s.DTCID = globalRocID / number_of_crv_rocs_per_dtc;
  uint8_t evbMode = 0; // ask Eric
  crvData.header.s.EventWindowMode = evbMode;

  //------------------
  // CRVROCStatusPacket
  //------------------
  // Word 0
  crvData.rocStatus.unused1 = 0;
  crvData.rocStatus.PacketType = 0x06;
  crvData.rocStatus.ControllerID =
      globalRocID % number_of_crv_rocs_per_dtc; // TODO: Is this correct?
  // Word 1
  crvData.rocStatus.ControllerEventWordCount =
      sizeof(CRVROCStatusPacket) +
      sizeof(CRVHitReadoutPacket) * nHits; // TODO: ArtFragmentReader::GetCRVHitCount() seems to
                                           // interpret this as byte counter and not as word count
  // Word 2
  crvData.rocStatus.ActiveFEBFlags2 = 0xFF;
  crvData.rocStatus.unused2 = 0;
  // Word 3
  crvData.rocStatus.ActiveFEBFlags0 = 0xFF;
  crvData.rocStatus.ActiveFEBFlags1 = 0xFF;
  // Word 3
  crvData.rocStatus.unused3 = 0;
  crvData.rocStatus.unused4 = 0;
  // Word 4
  crvData.rocStatus.TriggerCount =
      nHits; // TODO: Is this is what is meant by TriggerCount? Why isn't this number used in
             // ArtFragmentReader::GetCRVHitCount()?
  // Word 5
  crvData.rocStatus.unused5 = 0;
  crvData.rocStatus.unused6 = 0;
  // Word 6
  crvData.rocStatus.Errors = 0x0;
  crvData.rocStatus.EventType = 0; // TODO: How is this defined?
}

//--------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvDMABlocks(DTCLib::DTC_Event& currentEvent,
                                                 const crv_data_block_list_t& crvDataBlocks) {
  // Loop over all ROCs
  uint8_t currentDTCID = 0;
  for (uint8_t globalRocID = 0; globalRocID < number_of_crv_rocs; globalRocID++) {
    // Add the current DataBlock to the current SuperBlock
    // curDataBlock.setTimestamp(ts); // Overwrite the timestamp
    const CrvDataPacket& crvData = crvDataBlocks.at(globalRocID);
    fillCrvDataStream(currentEvent, crvData);

    if (_diagLevel > 1) {
      if (globalRocID == 0 || currentDTCID != crvData.header.s.DTCID) {
        std::cout << "================================================" << std::endl;
        // std::cout << "\t\tTimestamp: " << ts << std::endl;
        std::cout << "\t\tDTCID: " << (int)crvData.header.s.DTCID << std::endl;
        std::cout << "\t\tSYSID: " << (int)crvData.header.s.SubsystemID << std::endl;
        currentDTCID = crvData.header.s.DTCID;
      }
      if (crvData.header.s.PacketCount > 0) {
        printHeader(crvData.header);
        printCrvData(crvData);
      }
    }

  } // End loop over DataBlocks
}

//--------------------------------------------------------------------------------
//  method to fill the datastream with the crv packets
//--------------------------------------------------------------------------------
void ArtBinaryPacketsFromDigis::fillCrvDataStream(DTCLib::DTC_Event& currentEvent,
                                                  const CrvDataPacket& crvData) {
  size_t sz =
      crvData.header.s.TransferByteCount; // byte count was increased to get full chunks of 16 bytes

  assert(sz < sizeof(mu2e_databuff_t));

  uint8_t dtcID = crvData.header.s.DTCID;
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
    memcpy(thisBlock.allocBytes->data() + pos, &(crvData.hits[i]), sizeof(CRVHitReadoutPacket));
    pos += sizeof(CRVHitReadoutPacket);
  }

  if (hitCount > 0) {
    putBlockInEvent(currentEvent, dtcID, DTCLib::DTC_Subsystem_CRV, thisBlock);
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ArtBinaryPacketsFromDigis);
