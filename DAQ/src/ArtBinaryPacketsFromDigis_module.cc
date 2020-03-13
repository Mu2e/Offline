//
// EDProducer module for converting straw hit digis into DTC formatted
// packets (which are stored in DataBlock data products)
//
//
// Original author Tomonari Miyashita
//
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

#include <math.h>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"
//geometry
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

// mu2e-artdaq-core includes
#include "mu2e-artdaq-core/Overlays/ArtFragment.hh"
#include "mu2e-artdaq-core/Overlays/ArtFragmentReader.hh"

//pci_linux_kernel_module includes
#include "dtcInterfaceLib/DTC.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
//#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>


// Definitions and typedefs needed for compatibility with HLS codeblock
#define NUM_PRESAMPLES 4
#define START_SAMPLES 5 // 0 indexed
#define NUM_SAMPLES 15
#define LOWER_TDC 16000
#define UPPER_TDC 64000

typedef uint16_t tdc_type;
typedef uint8_t  tot_type;
typedef uint16_t adc_type;
typedef uint16_t calib_constant_type;
typedef uint8_t  flag_mask_type;
typedef uint64_t timestamp;

using namespace std;

using  DataBlockHeader = mu2e::ArtFragmentReader::DataBlockHeader;
using  TrackerDataPacket = mu2e::ArtFragmentReader::TrackerDataPacket;
using  adc_t = mu2e::ArtFragment::adc_t;
using  CalorimeterDataPacket = mu2e::ArtFragmentReader::CalorimeterDataPacket;
using  CalorimeterBoardID = mu2e::ArtFragmentReader::CalorimeterBoardID;
using  CalorimeterHitReadoutPacket = mu2e::ArtFragmentReader::CalorimeterHitReadoutPacket;

//data struct for the calorimeter
struct CaloDataPacket {
  CalorimeterDataPacket                    dataPacket;
  CalorimeterBoardID                       boardID;
  std::vector<CalorimeterHitReadoutPacket> hitPacketVec;
  std::vector< std::vector<adc_t> >        waveformVec;
  std::vector<uint16_t> hitIndex;
};

typedef std::deque<std::pair<mu2e_databuff_t, size_t>> raw_data_list_t;
typedef std::pair<DataBlockHeader, TrackerDataPacket> tracker_data_block_t;
typedef std::pair<DataBlockHeader, CaloDataPacket> calo_data_block_t;
// typedef std::pair<DataBlockHeader,CRVDataPacket> crv_data_block_t;

typedef std::deque<tracker_data_block_t> tracker_data_block_list_t;
typedef std::deque<calo_data_block_t> calo_data_block_list_t;

namespace mu2e {

  constexpr int format_version = 6;
  //--------------------------------------------------------------------
  //
  //
  class ArtBinaryPacketsFromDigis : public art::EDProducer {
  public:

    explicit ArtBinaryPacketsFromDigis(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;
    virtual void beginRun(art::Run&);

    virtual void endJob();

    void produce(art::Event&) override;

    flag_mask_type filter(TrackerDataPacket& trkData,
			  calib_constant_type clockstart,
			  calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
			  calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
			  calib_constant_type gain_RSHIFT15,
			  calib_constant_type inverse_ionization_energy_LSHIFT26);

  private:

    size_t _generateTimestampTable;

    // Table used for mapping between DTC timestamp and art EventID
    std::string _tableFile;
    std::vector< std::pair<timestamp, timestamp> > tsTable;

    size_t  _timestampOffset;
    size_t  _numWordsWritten;
    size_t  _numEventsProcessed;
    int     _includeTracker;
    int     _includeCalorimeter;
    int     _includeCosmicRayVeto;

    int _includeDMAHeaders;

    // Set to 1 to save packet data to a binary file
    int _generateBinaryFile;

    string                _outputFile;
    ofstream              outputStream;

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
    // 6 rocs per DTC => 27 DTCs
    // 172 rocs * 8 crystals per roc => 1376
    // Note: the highest crystal ID in the old simulation was 1355
    const size_t number_of_calo_rocs = 172;
    const size_t number_of_crystals_per_roc = 8;
    const size_t number_of_calo_rocs_per_dtc = 6;

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    int _enableFPGAEmulation;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    art::ProductToken<StrawDigiCollection> const  _sdtoken;
    art::ProductToken<CaloDigiCollection>  const  _cdtoken;

    const Calorimeter* _calorimeter; // cached pointer to the calorimeter geometry

    void   fillEmptyHeaderDataPacket(DataBlockHeader& HeaderData, uint64_t& EventNum, uint8_t& ROCId, uint8_t& DTCId, uint8_t Subsys);
    void   printHeader(DataBlockHeader const& headerDataBlock);

    //--------------------------------------------------------------------------------
    //  methods used to process the tracker data
    //--------------------------------------------------------------------------------
    void   fillTrackerDataPacket(const StrawDigi& SD, TrackerDataPacket& TrkData);

    void   fillTrackerHeaderDataPacket(const StrawDigi& SD, DataBlockHeader& HeaderData, uint64_t& EventNum);

    void   fillEmptyTrackerDataPacket(TrackerDataPacket& TrkData);

    void   processTrackerData(art::Event& evt, uint64_t& eventNum,
			      tracker_data_block_list_t &trackerData);


    void   fillTrackerDMABlocks(raw_data_list_t& dataStream, tracker_data_block_list_t const& trackerData);

    void   fillTrackerDataStream(raw_data_list_t& dataStream, tracker_data_block_t const& trkData);

    void   printTrackerData(TrackerDataPacket const& curDataBlock);

    //--------------------------------------------------------------------------------
    //  methods used to handle the calorimeter data
    //--------------------------------------------------------------------------------
    void   fillCalorimeterDataPacket(const CaloDigi& SD, CaloDataPacket& caloData);
    void   addCaloHitToCaloPacket(calo_data_block_t& dataBlock,
				  CaloDataPacket& caloData);

    void   fillCalorimeterHeaderDataPacket(const CaloDigi& SD, DataBlockHeader& HeaderData, uint64_t& EventNum);

    void   fillEmptyCalorimeterDataPacket(CaloDataPacket& caloData);
    void   fillHeaderByteAndPacketCounts(calo_data_block_t& caloData);

    void   processCalorimeterData(art::Event& evt, uint64_t& eventNum,
				  calo_data_block_list_t& caloDataBlocks);


    void   fillCalorimeterDMABlocks(raw_data_list_t& dataStream, calo_data_block_list_t& caloData);

    void   fillCalorimeterDataStream(raw_data_list_t& dataStream, calo_data_block_t& caloData);

    void   printCalorimeterData(CaloDataPacket const& curDataBlock);

    size_t waveformMaximumIndex(std::vector<adc_t>& waveform);

    void flushBuffer(raw_data_list_t& buffers);
    void closeDataBuffer(raw_data_list_t& dataStream, bool openNew = true) {     
      uint64_t sz = dataStream.back().second;
      if (_diagLevel > 3) {
	std::cout<< "[closeDataBuffer] 1) sz = "<<sz<<std::endl;
	std::cout << "[closeDataBuffer] 1a) sz = "<< dataStream.back().second << std::endl;
      }
      if (sz < 48) {
	if (_diagLevel > 3) {
	  std::cout<< "[closeDataBuffer] dataStream.back().first = "<<dataStream.back().first<<std::endl;
	  std::cout<< "[closeDataBuffer] dataStream.back().first[sz] = "<<dataStream.back().first[sz]<<std::endl;
	  std::cout<< "[closeDataBuffer] &dataStream.back().first[sz] = "<<&dataStream.back().first[sz]<<std::endl;
	}
	bzero(&dataStream.back().first[sz], 48 - sz);
	sz = 48;
      }
      if (_diagLevel > 3) {
	std::cout << "[closeDataBuffer] dataStream.size = "<< dataStream.size() << std::endl;
	if (dataStream.size()> 0){
	  std::cout <<  "[closeDataBuffer] dataStream.back().second = " << dataStream.back().second <<std::endl;
	}
	std::cout<< "[closeDataBuffer] 2) sz = "<<sz<<std::endl;
      }
      uint64_t ex_sz = sz - 16;

      memcpy(&dataStream.back().first, &sz, sizeof(uint64_t));
      memcpy(&dataStream.back().first[8], &ex_sz, sizeof(uint64_t));

      if (openNew) {
	dataStream.push_back(std::pair<mu2e_databuff_t, size_t>());
	dataStream.back().second = 16; // DMA Size words
      }
    }

    std::vector<adc_t> generateDMABlockHeader(size_t theCount) const;
    std::vector<adc_t> generateEventByteHeader(size_t theCount) const;
  };

  //--------------------------------------------------------------------------------
  // temporary function used to find the location of the waveform peak in the 
  // calorimeter digitized waveform
  //--------------------------------------------------------------------------------
  size_t ArtBinaryPacketsFromDigis::waveformMaximumIndex(std::vector<adc_t>& waveform) {
    size_t  indexMax(0), content(0);
    for (size_t i = 0; i < waveform.size(); ++i) {
      if (waveform[i] > content) {
	content = waveform[i];
	indexMax = i;
      }
    }

    return indexMax;
  }

  void   ArtBinaryPacketsFromDigis::printHeader(DataBlockHeader const& headerDataBlock) {
    printf("[ArtBinaryPacketsFromDigis::printHeader] START header print  \n");
    printf("[ArtBinaryPacketsFromDigis::printHeader] ByteCount      : %i \n", headerDataBlock.ByteCount);
    printf("[ArtBinaryPacketsFromDigis::printHeader] Hopcount       : %i \n", headerDataBlock.Hopcount);
    printf("[ArtBinaryPacketsFromDigis::printHeader] PacketType     : %i \n", headerDataBlock.PacketType);
    printf("[ArtBinaryPacketsFromDigis::printHeader] ROCID 	   : %i \n", headerDataBlock.ROCID);
    printf("[ArtBinaryPacketsFromDigis::printHeader] unused1	   : %i \n", headerDataBlock.unused1);
    printf("[ArtBinaryPacketsFromDigis::printHeader] SubsystemID    : %i \n", headerDataBlock.SubsystemID);
    printf("[ArtBinaryPacketsFromDigis::printHeader] Valid 	   : %i \n", headerDataBlock.Valid);
    printf("[ArtBinaryPacketsFromDigis::printHeader] PacketCount    : %i \n", headerDataBlock.PacketCount);
    printf("[ArtBinaryPacketsFromDigis::printHeader] unused2        : %i \n", headerDataBlock.unused2);
    printf("[ArtBinaryPacketsFromDigis::printHeader] TimestampLow   : %i \n", headerDataBlock.TimestampLow);
    printf("[ArtBinaryPacketsFromDigis::printHeader] TimestampMed   : %i \n", headerDataBlock.TimestampMed);
    printf("[ArtBinaryPacketsFromDigis::printHeader] TimestampHigh  : %i \n", headerDataBlock.TimestampHigh);
    printf("[ArtBinaryPacketsFromDigis::printHeader] Status	   : %i \n", headerDataBlock.Status);
    printf("[ArtBinaryPacketsFromDigis::printHeader] FormatVersion  : %i \n", headerDataBlock.FormatVersion);
    printf("[ArtBinaryPacketsFromDigis::printHeader] DTCID	   : %i \n", headerDataBlock.DTCID);
    printf("[ArtBinaryPacketsFromDigis::printHeader] EVBMode        : %i \n", headerDataBlock.EVBMode);

  }

  void   ArtBinaryPacketsFromDigis::printTrackerData(TrackerDataPacket const& trkData) {
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] START tracker-data print \n");
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] StrawIndex    : %i \n", (int)trkData.StrawIndex);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC0		: %i \n", (int)trkData.TDC0);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC1		: %i \n", (int)trkData.TDC1);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT0		: %i \n", (int)trkData.TOT0);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT1		: %i \n", (int)trkData.TOT1);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC00         : %i \n", (int)trkData.ADC00);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC01  	: %i \n", (int)trkData.ADC01());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC02  	: %i \n", (int)trkData.ADC02());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC03  	: %i \n", (int)trkData.ADC03);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC04  	: %i \n", (int)trkData.ADC04);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC05 	: %i \n", (int)trkData.ADC05());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC06  	: %i \n", (int)trkData.ADC06());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC07  	: %i \n", (int)trkData.ADC07);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC08  	: %i \n", (int)trkData.ADC08);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC09  	: %i \n", (int)trkData.ADC09());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC10  	: %i \n", (int)trkData.ADC10());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC11  	: %i \n", (int)trkData.ADC11);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC12  	: %i \n", (int)trkData.ADC12);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC13 	: %i \n", (int)trkData.ADC13());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC14 	: %i \n", (int)trkData.ADC14());
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] unused1 	: %i \n", (int)trkData.unused1);
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] PreprocessingFlags : %i \n", (int)trkData.PreprocessingFlags);


  }

  void   ArtBinaryPacketsFromDigis::printCalorimeterData(CaloDataPacket const& caloData) {
    CalorimeterDataPacket packet = caloData.dataPacket;
    CalorimeterBoardID    boardId = caloData.boardID;
    size_t      nHits = caloData.hitPacketVec.size();
    printf("[ArtBinaryPacketsFromDigis::printCaloData] START calorimeter-data print \n");
    printf("[ArtBinaryPacketsFromDigis::printCaloData] NumberofHits        : %i \n", (int)packet.NumberOfHits);
    printf("[ArtBinaryPacketsFromDigis::printCaloData] BoardID             : %i \n", (int)boardId.BoardID);
    printf("[ArtBinaryPacketsFromDigis::printCaloData] ChannelStatusFlagsA : %i \n", (int)boardId.ChannelStatusFlagsA);
    printf("[ArtBinaryPacketsFromDigis::printCaloData] ChannelStatusFlagsB : %i \n", (int)boardId.ChannelStatusFlagsB);
    printf("[ArtBinaryPacketsFromDigis::printCaloData] unused              : %i \n", (int)boardId.unused);
    printf("[ArtBinaryPacketsFromDigis::printCaloData] NHits               : %i \n", (int)nHits);

    for (size_t i = 0; i < nHits; ++i) {
      CalorimeterHitReadoutPacket const& hit = caloData.hitPacketVec[i];
      printf("[ArtBinaryPacketsFromDigis::printCaloData] hit : %i \n", (int)i);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] ChannelNumber : %i \n", (int)hit.ChannelNumber);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] DIRACA        : %i \n", (int)hit.DIRACA);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] DIRACB        : %i \n", (int)hit.DIRACB);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] ErrorFlags    : %i \n", (int)hit.ErrorFlags);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] Time          : %i \n", (int)hit.Time);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] NumberOfSamples : %i \n", (int)hit.NumberOfSamples);
      printf("[ArtBinaryPacketsFromDigis::printCaloData] IndexOfMaxDigitizerSample : %i \n", (int)hit.IndexOfMaxDigitizerSample);
    }

  }

  void   ArtBinaryPacketsFromDigis::fillTrackerDataStream(raw_data_list_t& dataStream,
							  tracker_data_block_t const& trackerData) {

    auto sz = sizeof(DataBlockHeader);
    //check that the trkDataBlock is not empty

    assert(trackerData.first.PacketCount == 0 || trackerData.first.PacketCount == 2); // Tracker DataBlocks have 0 or 2 DataPackets
    if (trackerData.first.PacketCount > 0) {
      sz += sizeof(TrackerDataPacket);
    }


    if (sz + dataStream.back().second >= sizeof(mu2e_databuff_t)) {
      closeDataBuffer(dataStream);
    }

    if (_diagLevel > 3) {
      std::cout << "[ArtBinaryPacketsFromDigis::fillTrackerDataStream] 1) sizeof(dataStream.back().second) = "<<sizeof(dataStream.back().second)<<std::endl;
      std::cout << "[ArtBinaryPacketsFromDigis::fillTrackerDataStream] 1a) dataStream.back().second = "<< dataStream.back().second << std::endl;
    }

    auto pos = dataStream.back().second;
    memcpy(&dataStream.back().first[pos], &trackerData.first, sizeof(DataBlockHeader));
    pos += sizeof(DataBlockHeader);

    if (trackerData.first.PacketCount > 0) {
      assert(sizeof(TrackerDataPacket) % 16 == 0); // Make sure that TrackerDataPacket is an even number of DataPackets!
      memcpy(&dataStream.back().first[pos], &trackerData.second, sizeof(TrackerDataPacket));
      pos += sizeof(TrackerDataPacket);
    }

    dataStream.back().second = pos;


  }

  void  ArtBinaryPacketsFromDigis::fillTrackerDMABlocks(raw_data_list_t& dataStream, tracker_data_block_list_t const& trkData) {

    auto curDTCID = trkData.front().first.DTCID;
    bool first = true;
    if (_diagLevel > 1) {
      std::cout << "[ArtBinaryPacketsFromDigis::fillTrackerDMABlocks] trkData.size() = " << trkData.size() << std::endl;
    }
    for (auto& dataBlock : trkData) {

      fillTrackerDataStream(dataStream, dataBlock);

      if (_diagLevel > 1) {
	if (dataBlock.first.DTCID != curDTCID || first) {
	  std::cout << "================================================" << std::endl;
	  //std::cout << "\t\tTimestamp: " << ts << std::endl;
	  std::cout << "\t\tDTCID: " << (int)dataBlock.first.DTCID << std::endl;
	  std::cout << "\t\tSYSID: " << (int)dataBlock.first.SubsystemID << std::endl;
	  curDTCID = dataBlock.first.DTCID;
	  first = false;
	}
	if (dataBlock.first.PacketCount > 0) {
	  printHeader(dataBlock.first);
	  printTrackerData(dataBlock.second);
	}
      }


    } // End loop over DataBlocks
  }

  //--------------------------------------------------------------------------------
  //
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillEmptyTrackerDataPacket(TrackerDataPacket& TrkData) {
    bzero(&TrkData, sizeof(TrkData));
  }



  //--------------------------------------------------------------------------------  
  // 
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillEmptyHeaderDataPacket(DataBlockHeader& headerData, uint64_t& EventNum,
							      uint8_t& ROCId, uint8_t& DTCId, uint8_t Subsys) {

    // Fill in the byte count field of the header packet
    // Word 0
    headerData.ByteCount = sizeof(DataBlockHeader);
    // Word 1
    headerData.Hopcount = 0;//ask Eric!!!//FIX ME!
    headerData.PacketType = 5;//PacketType::Dataheader; 

    headerData.ROCID = ROCId;
    headerData.unused1 = 0;//ask Eric!

    headerData.SubsystemID = Subsys;//DTCLib::DTC_Subsystem_Tracker; //: 3;

    headerData.Valid = 1;
    // Word 2
    headerData.PacketCount = 0;
    headerData.unused2 = 0;// : 5;
    // Word 3
    uint64_t timestamp = EventNum;
    headerData.TimestampLow = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    headerData.TimestampMed = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    headerData.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    headerData.Status = 0; // 0 corresponds to "TimeStamp had valid data"
    headerData.FormatVersion = format_version;
    // Word 7
    headerData.DTCID = DTCId;
    uint8_t  evbMode = 0;//maybe off-spill vs on-spill?
    headerData.EVBMode = evbMode;
  }

  //--------------------------------------------------------------------------------
  // create the header for the StrawPacket
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillTrackerHeaderDataPacket(const StrawDigi& SD, DataBlockHeader& HeaderData,
								uint64_t& EventNum) {
    // Word 0
    adc_t   nBytes = sizeof(DataBlockHeader) + sizeof(TrackerDataPacket);//ask Eric! //FIX ME!
    HeaderData.ByteCount = nBytes;
    // Word 1
    HeaderData.Hopcount = 0;//currently unused
    HeaderData.PacketType = 5;//PacketType::Dataheader;

    // 96 straws per ROC/panel
    // 6 panels / plane
    //      // 36 planes => 216 panels/ROCs
    // There are actually 40 planes in the input file => 240 panels/ROCs    
    int panel = SD.strawId().getPanel();
    int plane = SD.strawId().getPlane();

    // ROC ID, counting from 0 across all DTCs (for the tracker)
    //    uint8_t localROCID = panel;
    uint8_t globalROCID = (plane * 6) + panel;

    // 240 ROCs total
    if (globalROCID >= number_of_rocs) {
      throw cet::exception("DATA") << " Global ROC ID " << globalROCID
				   << " exceeds limit of " << number_of_rocs;
    }
    HeaderData.ROCID = panel;
    HeaderData.unused1 = 0;
    HeaderData.SubsystemID = DTCLib::DTC_Subsystem_Tracker;
    HeaderData.Valid = 1;
    // Word 2
    HeaderData.PacketCount = 2;
    HeaderData.unused2 = 0;
    // Word 3
    uint64_t timestamp = EventNum;
    HeaderData.TimestampLow = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    HeaderData.TimestampMed = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    HeaderData.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    HeaderData.Status = 0; // 0 corresponds to "TimeStamp had valid data"
    HeaderData.FormatVersion = format_version;
    // Word 7
    HeaderData.DTCID = static_cast<uint8_t>(globalROCID / number_of_rocs_per_dtc);
    uint8_t  evbMode = 0;//ask Eric
    HeaderData.EVBMode = evbMode;

  }

  void   ArtBinaryPacketsFromDigis::fillTrackerDataPacket(const StrawDigi& SD, TrackerDataPacket& TrkData) {

    TrkData.StrawIndex = SD.strawId().asUint16();
    TrkData.TDC0 = SD.TDC(StrawEnd::cal);
    TrkData.TDC1 = SD.TDC(StrawEnd::hv);
    TrkData.TOT0 = SD.TOT(StrawEnd::cal);
    TrkData.TOT1 = SD.TOT(StrawEnd::hv);

    TrkTypes::ADCWaveform const& theWaveform = SD.adcWaveform();

    TrkData.SetWaveform(0,  theWaveform[0]);
    TrkData.SetWaveform(1,  theWaveform[1]);
    TrkData.SetWaveform(2,  theWaveform[2]);
    TrkData.SetWaveform(3,  theWaveform[3]);
    TrkData.SetWaveform(4,  theWaveform[4]);
    TrkData.SetWaveform(5,  theWaveform[5]);
    TrkData.SetWaveform(6,  theWaveform[6]);
    TrkData.SetWaveform(7,  theWaveform[7]);  
    TrkData.SetWaveform(8,  theWaveform[8]);  
    TrkData.SetWaveform(9,  theWaveform[9]);  
    TrkData.SetWaveform(10, theWaveform[10]);  
    TrkData.SetWaveform(11, theWaveform[11]); 
    TrkData.SetWaveform(12, theWaveform[12]); 
    TrkData.SetWaveform(13, theWaveform[13]); 
    TrkData.SetWaveform(14, theWaveform[14]); 

    TrkData.PreprocessingFlags = 0;

    //    TrkData.unused1 = ;
    if (_enableFPGAEmulation) {
      // Note: Eventually, there will be a calibration database to provide these sorts of values
      // For now, these are just placeholders
      calib_constant_type clockstart = 320; //10 / .03125 (tdclsb)
      calib_constant_type maxenergyls8 = 583, minenergyls8 = 0;
      calib_constant_type gainrs15 = 1389;
      calib_constant_type inverseionizationenergyls26 = 633;
      calib_constant_type panelTDCoffset = 0, hvoffset = 0, caloffset = 0;


      flag_mask_type f = filter(TrkData,
				clockstart, panelTDCoffset, hvoffset, caloffset,
				maxenergyls8, minenergyls8, gainrs15, inverseionizationenergyls26);

      TrkData.PreprocessingFlags = f & 0xFF;
    }
  }


  ArtBinaryPacketsFromDigis::ArtBinaryPacketsFromDigis(fhicl::ParameterSet const& pset) :
    art::EDProducer{ pset },
		       _generateTimestampTable(pset.get<size_t>("generateTimestampTable", 0)),
		       _tableFile(pset.get<std::string>("tableFile", "tsTable.bin")),
		       _timestampOffset(pset.get<size_t>("timestampOffset", 0)),
		       _numWordsWritten(0),
		       _numEventsProcessed(0),
		       _includeTracker(pset.get<int>("includeTracker", 1)),
		       _includeCalorimeter(pset.get<int>("includeCalorimeter", 1)),
		       _includeCosmicRayVeto(pset.get<int>("includeCosmicRayVeto", 0)),
		       _includeDMAHeaders(pset.get<int>("includeDMAHeaders", 1)),
		       _generateBinaryFile(pset.get<int>("generateBinaryFile", 1)),
		       _outputFile(pset.get<string>("outputFile", "DTC_packets.bin")),
		       _generateTextFile(pset.get<int>("generateTextFile", 0)),
		       _diagLevel(pset.get<int>("diagLevel", 0)),
		       _enableFPGAEmulation(pset.get<int>("enableFPGAEmulation", 0)),
		       _maxFullPrint(pset.get<int>("maxFullPrint", 5)),
		       _sdtoken{ consumes<StrawDigiCollection>(pset.get<std::string>("StrawDigiCollection","makeSD")) },
		       _cdtoken{ consumes<CaloDigiCollection>(pset.get<std::string>("CaloDigiCollection","CaloDigiFromShower")) }{

			 produces<timestamp>();
			 //		produces< raw_data_list_t >();

			 if (_generateBinaryFile == 1) {
			   outputStream.open(_outputFile, std::ios::out | std::ios::binary);
			 }
		       }

  void ArtBinaryPacketsFromDigis::beginJob() {

    if (_diagLevel > 0) {
      cout << "ArtBinaryPacketsFromDigis Diaglevel: "
	   << _diagLevel << " "
	   << _maxFullPrint
	   << endl;

      if (_enableFPGAEmulation) {
	cout << "FPGA Emulation Enabled" << endl;
      }

    }

  }

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

	if (_diagLevel > 3) {
	  std::cout << "TIMESTAMP_MAPPING: timestamp: "
		    << tsTable[idx].first
		    << " uniqueid: "
		    << tsTable[idx].second
		    << std::endl;
	}

      }
      tsTableStream << std::flush;
      tsTableStream.close();
    }

    if (_diagLevel > 0) {
      std::cout << "BinaryPacketsFromDataBlocks: "
		<< "Finished writing "
		<< _numWordsWritten
		<< " words from "
		<< _numEventsProcessed
		<< " events to "
		<< _outputFile
		<< std::endl;
    }

  }


  void ArtBinaryPacketsFromDigis::flushBuffer(raw_data_list_t& buffers) {
    closeDataBuffer(buffers, false);
    for (size_t idx = 0; idx < buffers.size(); idx++) {
      outputStream.write(reinterpret_cast<const char*>(&(buffers[idx].first)), buffers[idx].second);
      _numWordsWritten += buffers[idx].second;
    }
    outputStream << std::flush;
  }

  void ArtBinaryPacketsFromDigis::produce(art::Event& evt) {

    // unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

    uint64_t eventNum = evt.id().event();//is not unique! internal counter???//FIXME!
    uint64_t ts = _numEventsProcessed + _timestampOffset;

    if (_diagLevel > 2) {
      cout << "ArtBinaryPacketsFromDigis: eventNum: " << eventNum << endl;
    }

    tracker_data_block_list_t trackerData;
    calo_data_block_list_t caloData;

    if (_includeTracker > 0) {
      processTrackerData(evt, ts,
			 trackerData);
    }

    if (_includeCalorimeter > 0) {
      processCalorimeterData(evt, ts,
			     caloData);
    }

    // Break the DataBlocks into DMABlocks and add DMABlock headers
    raw_data_list_t dataStream;
    dataStream.push_back(std::pair<mu2e_databuff_t,size_t>());
    dataStream.back().second = 16;

    if (_includeTracker > 0) {

      fillTrackerDMABlocks(dataStream, trackerData);
    }

    if (_includeCalorimeter > 0) {

      fillCalorimeterDMABlocks(dataStream, caloData);
    }

    // Write all values, including superblock header and DMA header values, to output buffer
    if (_generateBinaryFile == 1) {
      flushBuffer(dataStream);
    }

    _numEventsProcessed += 1;

    // Store the timestamp and DataBlockCollection in the event
    evt.put(std::unique_ptr<timestamp>(new timestamp(ts)));
    //		evt.put(std::make_unique< raw_data_list_t >(dataStream));

  } // end of ::produce


  // method....
  void ArtBinaryPacketsFromDigis::processCalorimeterData(art::Event& evt, uint64_t& eventNum, calo_data_block_list_t& caloDataBlocks) {
    auto  const& cdH = evt.getValidHandle(_cdtoken);
    const CaloDigiCollection& hits_CD(*cdH);

    calo_data_block_list_t tmpCaloDataBlockList;

    for (size_t i = 0; i < hits_CD.size(); ++i) {
      CaloDigi const& CD = hits_CD.at(i);

      // Fill struct with info for current hit
      DataBlockHeader   headerData;
      fillCalorimeterHeaderDataPacket(CD, headerData, eventNum);
      CaloDataPacket    caloData;
      fillCalorimeterDataPacket(CD, caloData);

      tmpCaloDataBlockList.push_back(std::pair<DataBlockHeader, CaloDataPacket>(headerData, caloData));
    }

    if (_diagLevel > 1) {
      std::cout << "[ArtBinaryPacketsFromDigis::processCalorimeterData ] Total number of calorimeter non-empty DataBlocks = " <<
	tmpCaloDataBlockList.size() << std::endl;
    }

    uint8_t max_dtc_id = number_of_calo_rocs / number_of_calo_rocs_per_dtc - 1;
    if (number_of_calo_rocs % number_of_calo_rocs_per_dtc > 0) {
      max_dtc_id += 1;
    }

    // Loop over the DTC/ROC pairs and generate datablocks for each ROC
    for (uint8_t dtcID = 0; dtcID < max_dtc_id; dtcID++) {

      for (uint8_t rocID = 0; rocID < number_of_calo_rocs_per_dtc; ++rocID) {
	// Find all hits for this event coming from the specified DTC/ROC combination
	bool   is_first(true);
	for (size_t curHitIdx = 0; curHitIdx < tmpCaloDataBlockList.size(); curHitIdx++) {
	  if (tmpCaloDataBlockList[curHitIdx].first.DTCID == dtcID &&
	      tmpCaloDataBlockList[curHitIdx].first.ROCID == rocID) {
	    if (is_first) {
	      is_first = false;
	      caloDataBlocks.push_back(tmpCaloDataBlockList[curHitIdx]);
	    }
	    else {
	      addCaloHitToCaloPacket(caloDataBlocks.back(), tmpCaloDataBlockList[curHitIdx].second);
	    }
	  }
	}

	if (is_first) {
	  // No hits, so just fill a header packet and no data packets
	  DataBlockHeader   headerData;
	  CaloDataPacket    caloData;

	  fillEmptyHeaderDataPacket(headerData, eventNum, rocID, dtcID, DTCLib::DTC_Subsystem_Calorimeter);

	  caloDataBlocks.push_back(std::pair<DataBlockHeader, CaloDataPacket>(headerData, caloData));
	}
	else {
	  fillHeaderByteAndPacketCounts(caloDataBlocks.back());
	}


      } //Done looping over the ROCs in a given DTC
    }

  }


  //--------------------------------------------------------------------------------
  // Fix header ByteCount and PacketCount fields
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillHeaderByteAndPacketCounts(calo_data_block_t& caloData)
  {
    caloData.first.ByteCount = 16 /*header packet*/ + sizeof(uint16_t) /* num hits */ + sizeof(CalorimeterBoardID) + (sizeof(uint16_t) + sizeof(CalorimeterHitReadoutPacket)) * caloData.second.hitPacketVec.size();

    auto idxPos = sizeof(uint16_t) + sizeof(CalorimeterBoardID) + sizeof(uint16_t) * caloData.second.hitPacketVec.size();
    for (auto& vec : caloData.second.waveformVec) {
      caloData.first.ByteCount += sizeof(adc_t) * vec.size();
      caloData.second.hitIndex.push_back(idxPos);

      idxPos += sizeof(CalorimeterHitReadoutPacket) + sizeof(adc_t) * vec.size();
    }

    while (caloData.first.ByteCount % 16 != 0) caloData.first.ByteCount++;

    caloData.first.PacketCount = (caloData.first.ByteCount - 16) / 16;
  }

  //--------------------------------------------------------------------------------
  // crate a caloPacket from the digi
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillCalorimeterDataPacket(const CaloDigi& CD,
							      CaloDataPacket& CaloData) {
    CaloData.dataPacket.NumberOfHits = 1;

    CalorimeterBoardID      ccBoardID;
    // ROC ID, counting from 0, across all (for the calorimeter)
    size_t crystalId = _calorimeter->caloInfo().crystalByRO(CD.roId());
    size_t globalROCID = crystalId / number_of_crystals_per_roc;

    ccBoardID.BoardID = globalROCID % number_of_calo_rocs_per_dtc;
    ccBoardID.ChannelStatusFlagsA = 0;
    ccBoardID.ChannelStatusFlagsB = 0;
    ccBoardID.unused = 0;

    CaloData.boardID = ccBoardID;

    CalorimeterHitReadoutPacket   hitPacket;
    hitPacket.ChannelNumber = CD.roId();
    hitPacket.DIRACA = 0;
    hitPacket.DIRACB = (((CD.roId() % 2) << 12) | (crystalId));
    hitPacket.ErrorFlags = 0;
    hitPacket.Time = CD.t0();
    std::vector<adc_t>      theWaveform;
    for (size_t i = 0; i < CD.waveform().size(); ++i) { theWaveform.push_back((adc_t)CD.waveform().at(i)); }
    hitPacket.NumberOfSamples = theWaveform.size();
    hitPacket.IndexOfMaxDigitizerSample = waveformMaximumIndex(theWaveform);
    CaloData.hitPacketVec.push_back(hitPacket);

    CaloData.waveformVec.push_back(theWaveform);

  }


  //--------------------------------------------------------------------------------
  // add a caloHit to a caloPacketVector
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::addCaloHitToCaloPacket(calo_data_block_t& caloDataBlock,
							   CaloDataPacket& caloHit) {
    caloDataBlock.second.dataPacket.NumberOfHits += 1;

    caloDataBlock.second.hitPacketVec.push_back(caloHit.hitPacketVec[0]);//hitPacket);
    caloDataBlock.second.waveformVec.push_back(caloHit.waveformVec[0]);

    //increase the size of the block in the header
    caloDataBlock.first.ByteCount += sizeof(uint16_t) * (caloHit.hitPacketVec[0].NumberOfSamples + 1) + sizeof(CalorimeterHitReadoutPacket);
    caloDataBlock.first.PacketCount = std::ceil((caloDataBlock.first.ByteCount - 16) / 16);
  }

  //--------------------------------------------------------------------------------
  //
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillCalorimeterDMABlocks(raw_data_list_t& dataStream, calo_data_block_list_t& caloData) {


    bool first = true;
    auto curDTCID = caloData.front().first.DTCID;
    for (size_t dataBlockIdx = 0; dataBlockIdx < caloData.size(); dataBlockIdx++) {

      // Add the current DataBlock to the current SuperBlock
      //curDataBlock.setTimestamp(ts); // Overwrite the timestamp

      fillCalorimeterDataStream(dataStream, caloData[dataBlockIdx]);

      if (_diagLevel > 1) {
	if (first || curDTCID != caloData[dataBlockIdx].first.DTCID) {
	  std::cout << "================================================" << std::endl;
	  //std::cout << "\t\tTimestamp: " << ts << std::endl;
	  std::cout << "\t\tDTCID: " << (int)caloData[dataBlockIdx].first.DTCID << std::endl;
	  std::cout << "\t\tSYSID: " << (int)caloData[dataBlockIdx].first.SubsystemID << std::endl;
	  first = false;
	  curDTCID = caloData[dataBlockIdx].first.DTCID;
	}
	if (caloData[dataBlockIdx].first.PacketCount > 0){
	  printHeader(caloData[dataBlockIdx].first);
	  printCalorimeterData(caloData[dataBlockIdx].second);
	}
      }

    } // End loop over DataBlocks
  }

  //--------------------------------------------------------------------------------
  //  method to fill the datastream with the calorimeter packets
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillCalorimeterDataStream(raw_data_list_t& dataStream, calo_data_block_t& caloData) {

    size_t  sz = sizeof(DataBlockHeader);
    //check that the trkDataBlock is not empty
    if (caloData.second.hitPacketVec.size() != 0) {
      sz += sizeof(CalorimeterDataPacket) + caloData.second.hitPacketVec.size() * sizeof(uint16_t) + sizeof(CalorimeterBoardID);
      for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {
	auto   nSamples = caloData.second.hitPacketVec[i].NumberOfSamples;
	sz += sizeof(uint16_t) * nSamples + sizeof(CalorimeterHitReadoutPacket);
      }
    }
    while (sz % 16 != 0) sz++;

    assert(sz < sizeof(mu2e_databuff_t));
    if (dataStream.back().second + sz >= sizeof(mu2e_databuff_t)) {
      closeDataBuffer(dataStream);
    }

    assert(sz == caloData.first.ByteCount);

    auto pos = dataStream.back().second;
    memcpy(&dataStream.back().first[pos], &caloData.first, sizeof(DataBlockHeader));
    pos += sizeof(DataBlockHeader);

    if (caloData.second.hitPacketVec.size() != 0) {
      //memcpy(&dataStream.back().first[pos], &(caloData.second.dataPacket), sizeof(CalorimeterDataPacket));//HERE!
      //      pos += sizeof(CalorimeterDataPacket);

      uint16_t hitCount = caloData.second.hitPacketVec.size();
      memcpy(&dataStream.back().first[pos], &hitCount, sizeof(uint16_t));
      pos += sizeof(uint16_t);

      memcpy(&dataStream.back().first[pos], &caloData.second.hitIndex[0], sizeof(uint16_t) * caloData.second.hitIndex.size());
      pos += sizeof(uint16_t) * caloData.second.hitIndex.size();

      memcpy(&dataStream.back().first[pos], &(caloData.second.boardID), sizeof(CalorimeterBoardID));
      pos += sizeof(CalorimeterBoardID);

      for (size_t i = 0; i < caloData.second.hitPacketVec.size(); ++i) {

	memcpy(&dataStream.back().first[pos], &(caloData.second.hitPacketVec[i]), sizeof(CalorimeterHitReadoutPacket));
	pos += sizeof(CalorimeterHitReadoutPacket);

	auto waveform_size = sizeof(uint16_t) * (caloData.second.waveformVec[i].size());
	memcpy(&dataStream.back().first[pos], &(caloData.second.waveformVec[i][0]), waveform_size);
	pos += waveform_size;
      }//end loop over the calorimeterHitReadoutPacketVector
    }

    while (pos % 16 != 0) pos++;
    
    dataStream.back().second = pos;
  }


  //--------------------------------------------------------------------------------
  // create the header for the caloPacket
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillCalorimeterHeaderDataPacket(const CaloDigi& CD,
								    DataBlockHeader& HeaderData,
								    uint64_t& EventNum) {
    // Word 0
    adc_t   nBytes = sizeof(DataBlockHeader) + sizeof(CalorimeterDataPacket) + sizeof(CalorimeterBoardID);//this needs to be increased every time a new hit is addeded!
    HeaderData.ByteCount = nBytes;
    // Word 1
    HeaderData.Hopcount = 0;//currently unused
    HeaderData.PacketType = 5;//PacketType::Dataheader;

    // ROC ID, counting from 0, across all (for the calorimeter)
    size_t crystalId = _calorimeter->caloInfo().crystalByRO(CD.roId());
    size_t globalROCID = crystalId / number_of_crystals_per_roc;

    HeaderData.ROCID = globalROCID % number_of_rocs_per_dtc;//currently unknown. FIXME!
    HeaderData.unused1 = 0;
    HeaderData.SubsystemID = DTCLib::DTC_Subsystem_Calorimeter;
    HeaderData.Valid = 1;
    // Word 2
    HeaderData.PacketCount = 1;//NEEDS TO BE INCREASED EVERY TIME A NEW HIT IS ADDED!
    HeaderData.unused2 = 0;
    // Word 3
    uint64_t timestamp = EventNum;
    HeaderData.TimestampLow = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    HeaderData.TimestampMed = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    HeaderData.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    HeaderData.Status = 0; // 0 corresponds to "TimeStamp had valid data"
    HeaderData.FormatVersion = format_version;
    // Word 7
    HeaderData.DTCID = static_cast<uint8_t>(globalROCID / number_of_rocs_per_dtc);
    uint8_t  evbMode = 0;//ask Eric
    HeaderData.EVBMode = evbMode;

  }

  //--------------------------------------------------------------------------------
  //  method that process the tracker data 
  //--------------------------------------------------------------------------------
  void ArtBinaryPacketsFromDigis::processTrackerData(art::Event& evt, uint64_t& eventNum,
						     tracker_data_block_list_t &trackerData) {
    auto  const& sdH = evt.getValidHandle(_sdtoken);
    const StrawDigiCollection& hits_SD(*sdH);

    tracker_data_block_list_t tmpTrackerData;

    for (size_t i = 0; i < hits_SD.size(); ++i) {
      StrawDigi const& SD = hits_SD.at(i);

      // Fill struct with info for current hit
      TrackerDataPacket trkData;
      fillTrackerDataPacket(SD, trkData);
      DataBlockHeader   headerData;
      fillTrackerHeaderDataPacket(SD, headerData, eventNum);

      tmpTrackerData.push_back(std::pair<DataBlockHeader, TrackerDataPacket>(headerData, trkData));
    }

    if (_diagLevel > 1) {
      std::cout << "[ArtBinaryPacketsFromDigis::processTrackerData ] Total number of tracker non-empty DataBlocks = " <<
	tmpTrackerData.size() << std::endl;
    }

    uint8_t max_dtc_id = number_of_rocs / number_of_rocs_per_dtc - 1;
    if (number_of_rocs % number_of_rocs_per_dtc > 0) {
      max_dtc_id += 1;
    }

    // Loop over the DTC/ROC pairs and generate datablocks for each ROC
    for (uint8_t dtcID = 0; dtcID < max_dtc_id; dtcID++) {

      for (uint8_t rocID = 0; rocID < number_of_rocs_per_dtc; ++rocID) {
	bool data_found = false;
	// Find all hits for this event coming from the specified DTC/ROC combination
	// if (_diagLevel > 1) {
	//   std::cout << "[ArtBinaryPacketsFromDigis::processTrackerData ] dtcId ="<< (int)dtcID << " rcoID = "<< (int)rocID << std::endl;
	//   std::cout << " DTC    ROC" <<std::endl;
	// }

	for (size_t curHitIdx = 0; curHitIdx < tmpTrackerData.size(); curHitIdx++) {
	  // if (_diagLevel > 1) {
	  //   std::cout << (int)tmpTrackerData[curHitIdx].first.DTCID << "    "<< (int)tmpTrackerData[curHitIdx].first.ROCID <<std::endl;
	  // }

	  if (tmpTrackerData[curHitIdx].first.DTCID == dtcID &&
	      tmpTrackerData[curHitIdx].first.ROCID == rocID) {
	    // if (_diagLevel > 1) {
	    //   std::cout <<"\t data found!" << std::endl;
	    // }
	    trackerData.push_back(tmpTrackerData[curHitIdx]);
	    data_found = true;
	  }
	}

	if (!data_found) {
	  // No hits, so just fill a header packet and no data packets
	  DataBlockHeader   headerData;
	  TrackerDataPacket trkData;

	  fillEmptyHeaderDataPacket(headerData, eventNum, rocID, dtcID, DTCLib::DTC_Subsystem_Tracker);
	  fillEmptyTrackerDataPacket(trkData);

	  trackerData.push_back(std::pair<DataBlockHeader, TrackerDataPacket>(headerData, trkData));
	}
      } //Done looping over the ROCs in a given DTC
    }
  }



  flag_mask_type ArtBinaryPacketsFromDigis::filter(TrackerDataPacket& TrkData,
						   calib_constant_type clockstart,
						   calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
						   calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
						   calib_constant_type gain_RSHIFT15,
						   calib_constant_type inverse_ionization_energy_LSHIFT26) {

    int failed_time = 0;
    int failed_energy = 0;
    int thistdc = TrkData.TDC0 < TrkData.TDC1 ? TrkData.TDC0 + hvoffset : TrkData.TDC1 + caloffset;
    thistdc = thistdc + clockstart + panelTDCoffset;
    if (thistdc < LOWER_TDC || thistdc > UPPER_TDC) {
      failed_time = 1;
    }
    std::array<adc_t, 15> const waveform = TrkData.Waveform();
    //sum up presamples
    int pedsum = 0;
    for (int i = 0; i < NUM_PRESAMPLES; i++) {
      pedsum += waveform[i];
    }
    int pedestal = pedsum / NUM_PRESAMPLES;

    int peak = 0;
    for (int i = START_SAMPLES; i < NUM_SAMPLES; i++) {
      if (TrkData.Waveform().at(i) > peak) {
	peak = waveform[i];
      }
      else {
	break;
      }
    }

    int energy = peak - pedestal;

    int energy_max_adjusted = ((((energy_max_LSHIFT8 * gain_RSHIFT15) >> 9)* inverse_ionization_energy_LSHIFT26) >> 10);
    int energy_min_adjusted = ((((energy_min_LSHIFT8 * gain_RSHIFT15) >> 9)* inverse_ionization_energy_LSHIFT26) >> 10);
    if (energy > energy_max_adjusted || energy < energy_min_adjusted) {
      failed_energy = 1;
    }

    return (failed_energy << 1) | failed_time;

  }

}


using mu2e::ArtBinaryPacketsFromDigis;
DEFINE_ART_MODULE(ArtBinaryPacketsFromDigis);
