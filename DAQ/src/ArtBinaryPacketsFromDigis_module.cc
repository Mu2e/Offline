//
// EDProducer module for converting tracker/calo/crv digis 
// into DTC formatted packets 
//
//
// Original author G. Pezzullo, E. Flumerfelt, and R. Ehrlich
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
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

// mu2e-artdaq-core includes
#include "mu2e-artdaq-core/Overlays/ArtFragment.hh"
#include "mu2e-artdaq-core/Overlays/ArtFragmentReader.hh"

//pci_linux_kernel_module includes
#include "dtcInterfaceLib/DTC.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
//#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>


// Typedefs needed for compatibility with HLS codeblock
using tdc_type = uint16_t;
using tot_type = uint8_t;
using adc_type = uint16_t;
using calib_constant_type = uint16_t;
using flag_mask_type = uint8_t;
using timestamp = uint64_t;

using namespace std;

using  DataBlockHeader = mu2e::ArtFragmentReader::DataBlockHeader;
using  TrackerDataPacket = mu2e::ArtFragmentReader::TrackerDataPacket;
using  adc_t = mu2e::ArtFragment::adc_t;
using  CalorimeterDataPacket = mu2e::ArtFragmentReader::CalorimeterDataPacket;
using  CalorimeterBoardID = mu2e::ArtFragmentReader::CalorimeterBoardID;
using  CalorimeterHitReadoutPacket = mu2e::ArtFragmentReader::CalorimeterHitReadoutPacket;
using  CRVROCStatusPacket = mu2e::ArtFragmentReader::CRVROCStatusPacket;
using  CRVHitReadoutPacket = mu2e::ArtFragmentReader::CRVHitReadoutPacket;

//data struct for the calorimeter
struct CaloDataPacket {
  CalorimeterDataPacket                    dataPacket;
  CalorimeterBoardID                       boardID;
  std::vector<CalorimeterHitReadoutPacket> hitPacketVec;
  std::vector< std::vector<adc_t> >        waveformVec;
  std::vector<uint16_t> hitIndex;
};

using raw_data_list_t      = std::deque<std::pair<mu2e_databuff_t, size_t>>; 
using tracker_data_block_t = std::pair<DataBlockHeader, TrackerDataPacket>;  
using calo_data_block_t    = std::pair<DataBlockHeader, CaloDataPacket>;     

using tracker_data_block_list_t = std::deque<tracker_data_block_t> ;
using calo_data_block_list_t    = std::deque<calo_data_block_t>;

//data struct for the crv
struct CrvDataPacket
{
  DataBlockHeader                  header;
  CRVROCStatusPacket               rocStatus;
  std::vector<CRVHitReadoutPacket> hits;
};

using crv_data_block_list_t = std::map<int,CrvDataPacket>;  //the map key is the CRV ROC ID

namespace mu2e {

  constexpr int format_version = 6;
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

    struct Config 
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<size_t>         generateTimestampTable{ Name("generateTimestampTable"), Comment("generate Timestamp table")};
      fhicl::Atom<std::string>    tableFile             { Name("tableFile"),              Comment("table file name")};
      fhicl::Atom<size_t>         timestampOffset       { Name("timestampOffset"),        Comment("timestamp Offset")};
      fhicl::Atom<int>            includeTracker        { Name("includeTracker"),         Comment("include Tracker digis")};
      fhicl::Atom<int>            includeCalorimeter    { Name("includeCalorimeter"),     Comment("include Calorimeter digis")};
      fhicl::Atom<int>            includeCrv            { Name("includeCrv"),             Comment("include Crv digis")};
      fhicl::Atom<int>            includeDMAHeaders     { Name("includeDMAHeaders"),      Comment("include DMA Headers")};
      fhicl::Atom<int>            generateBinaryFile    { Name("generateBinaryFile"),     Comment("generate BinaryFile")};
      fhicl::Atom<std::string>    outputFile            { Name("outputFile"),             Comment("output File name")};
      fhicl::Atom<int>            generateTextFile      { Name("generateTextFile"),       Comment("generate Text File")};
      fhicl::Atom<int>            diagLevel             { Name("diagLevel"),              Comment("diagnostic Level")};
      fhicl::Atom<int>            maxFullPrint          { Name("maxFullPrint"),           Comment("maxFullPrint")};
      fhicl::Atom<art::InputTag>  sdtoken               { Name("strawDigiCollection"),    Comment("Straw digi collection name") };
      fhicl::Atom<art::InputTag>  cdtoken               { Name("caloDigiCollection"),     Comment("Calo digi collection name") };
      fhicl::Atom<art::InputTag>  crvtoken              { Name("crvDigiCollection"),      Comment("Crv digi collection name") };
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
    std::vector< std::pair<timestamp, timestamp> > tsTable;

    size_t  _timestampOffset;
    int     _includeTracker;
    int     _includeCalorimeter;
    int     _includeCrv;

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
    art::ProductToken<StrawDigiCollection> const  _sdtoken;
    art::ProductToken<CaloDigiCollection>  const  _cdtoken;
    art::ProductToken<CrvDigiCollection>   const  _crvtoken;

    size_t  _numWordsWritten;
    size_t  _numEventsProcessed;


    const Calorimeter* _calorimeter; // cached pointer to the calorimeter geometry
    const CosmicRayShield* _crv;     // cached pointer to the crv geometry

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

    const size_t waveformMaximumIndex(std::vector<adc_t>const& waveform);

    //--------------------------------------------------------------------------------
    //  methods used to handle the crv data
    //-------------------------------------------------------------------------------- 

    void processCrvData(art::Event& evt, uint64_t& eventNum, crv_data_block_list_t& crvDataBlocks);
    uint8_t compressCrvDigi(int adc);
    void fillCrvDataPacket(const CrvDigi& digi, CRVHitReadoutPacket& hit, int& globalRocID); 
    void fillCrvHeaderPacket(CrvDataPacket& crvData, uint8_t globalRocID, uint64_t eventNum);
    void fillCrvDMABlocks(raw_data_list_t& dataStream, const crv_data_block_list_t& crvData);
    void fillCrvDataStream(raw_data_list_t& dataStream, const CrvDataPacket& crvData);
    void printCrvData(const CrvDataPacket &curDataBlock);

    //-------------------------------------------------------------------------------- 

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
  const size_t ArtBinaryPacketsFromDigis::waveformMaximumIndex(std::vector<adc_t>const & waveform) {
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

  void   ArtBinaryPacketsFromDigis::printCrvData(CrvDataPacket const& crvData) 
  {
    size_t  nHits = crvData.hits.size();
    printf("[ArtBinaryPacketsFromDigis::printCrvData] START crv-data print \n");
    printf("[ArtBinaryPacketsFromDigis::printCrvData] ROC controller ID   : %i \n", (int)crvData.rocStatus.ControllerID);
    printf("[ArtBinaryPacketsFromDigis::printCrvData] Errors              : %i \n", (int)crvData.rocStatus.Errors);
    printf("[ArtBinaryPacketsFromDigis::printCrvData] NHits               : %i \n", (int)nHits);

    for(size_t i = 0; i < nHits; ++i) 
    {
      printf("[ArtBinaryPacketsFromDigis::printCrvData] hit : %i \n", (int)i);
      printf("[ArtBinaryPacketsFromDigis::printCrvData] Channel       : %i \n", (int)(crvData.hits[i].SiPMID&0x7F));
      printf("[ArtBinaryPacketsFromDigis::printCrvData] FEB           : %i \n", (int)(crvData.hits[i].SiPMID>>7));
      printf("[ArtBinaryPacketsFromDigis::printCrvData] Time          : %i \n", (int)crvData.hits[i].HitTime);
      printf("[ArtBinaryPacketsFromDigis::printCrvData] NumOfSamples  : %i \n", (int)crvData.hits[i].NumSamples);
    }
  }


  void   ArtBinaryPacketsFromDigis::fillTrackerDataStream(raw_data_list_t& dataStream,
							  tracker_data_block_t const& trackerData) {

    auto sz = sizeof(DataBlockHeader);
    //check that the trkDataBlock is not empty

    if(trackerData.first.PacketCount != 0 && trackerData.first.PacketCount != 2){ // Tracker DataBlocks have 0 or 2 DataPackets
      throw cet::exception("Online-RECO")<<"ArtBinaryPacketsFromDigis::fillTrackerDataStream : trackerData.first.PacketCount == 0 || trackerData.first.PacketCount == 2)" << std::endl;
    }
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
      if(sizeof(TrackerDataPacket) % 16 != 0){ // Make sure that TrackerDataPacket is an even number of DataPackets!
	throw cet::exception("Online-RECO")<<"ArtBinaryPacketsFromDigis::fillTrackerDataStream : sizeof(TrackerDataPacket) % 16 == 0" << std::endl;
      }
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

    for (size_t i=0; i<TrkTypes::NADC; ++i){
      TrkData.SetWaveform(i,  theWaveform[i]);
    }

    TrkData.PreprocessingFlags = 0;

  }


  ArtBinaryPacketsFromDigis::ArtBinaryPacketsFromDigis(const art::EDProducer::Table<Config>& config):
    art::EDProducer{ config },
    _generateTimestampTable(config().generateTimestampTable()),
    _tableFile             (config().tableFile()),
    _timestampOffset       (config().timestampOffset()),
    _includeTracker        (config().includeTracker()),
    _includeCalorimeter    (config().includeCalorimeter()),
    _includeCrv            (config().includeCrv()),
    _includeDMAHeaders     (config().includeDMAHeaders()),
    _generateBinaryFile    (config().generateBinaryFile()),
    _outputFile            (config().outputFile()),
    _generateTextFile      (config().generateTextFile()),
    _diagLevel             (config().diagLevel()),
    _maxFullPrint          (config().maxFullPrint()),
    _sdtoken               { consumes<mu2e::StrawDigiCollection>(config().sdtoken())},
    _cdtoken               { consumes<mu2e::CaloDigiCollection> (config().cdtoken())},
    _crvtoken              { consumes<mu2e::CrvDigiCollection>  (config().crvtoken())},
    _numWordsWritten(0),
    _numEventsProcessed(0){

      produces<timestamp>();

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
    crv_data_block_list_t crvData;

    if (_includeTracker > 0) {
      processTrackerData(evt, ts,
			 trackerData);
    }

    if (_includeCalorimeter > 0) {
      processCalorimeterData(evt, ts,
			     caloData);
    }

    if (_includeCrv > 0) {
      processCrvData(evt, ts, crvData);
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

    if (_includeCrv > 0) {
      fillCrvDMABlocks(dataStream, crvData);
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

    if(sz >= sizeof(mu2e_databuff_t)){
      throw cet::exception("Online-RECO")<<"ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz < sizeof(mu2e_databuff_t)" << std::endl;
    }
    if (dataStream.back().second + sz >= sizeof(mu2e_databuff_t)) {
      closeDataBuffer(dataStream);
    }

    if(sz != caloData.first.ByteCount){
      throw cet::exception("Online-RECO")<<"ArtBinaryPacketsFromDigis::fillCalorimeterDataStream : sz == caloData.first.ByteCount" << std::endl;
    }

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

  //------------------------------------
  // Crv Methods
  //------------------------------------
  void ArtBinaryPacketsFromDigis::processCrvData(art::Event& evt, uint64_t& eventNum, crv_data_block_list_t& crvDataBlocks) 
  {
    auto  const& crvdH = evt.getValidHandle(_crvtoken);
    const CrvDigiCollection& digis(*crvdH);

    for(size_t i = 0; i < digis.size(); ++i) 
    {
      CrvDigi const& digi = digis.at(i);

      // Fill struct with info for current hit
      CRVHitReadoutPacket hit;
      int                 globalRocID;
      fillCrvDataPacket(digi, hit, globalRocID);
      crvDataBlocks[globalRocID].hits.push_back(hit);
    }

    if(_diagLevel > 1) 
    {
      std::cout << "[ArtBinaryPacketsFromDigis::processCrvData ] Total number of CRV digis = " << digis.size() << std::endl;
    }

    // Loop over all ROCs, fill headers for each ROC - even for ROCs without hits
    for(uint8_t globalRocID=0; globalRocID<number_of_crv_rocs; globalRocID++) 
    {
      fillCrvHeaderPacket(crvDataBlocks[globalRocID],globalRocID,eventNum); //this will create a new entry for ROCs without hits
    }
  }


  //--------------------------------------------------------------------------------
  // crate a crvPacket from the digi
  //--------------------------------------------------------------------------------
  uint8_t ArtBinaryPacketsFromDigis::compressCrvDigi(int adc)
  {
    //TODO: Temporary implementation until we have the real compression used at the FEBs
    adc-=95;
    if(adc<0) adc=0;
    uint8_t toReturn=adc;
    if(adc>50 && adc<=100) toReturn=50+(adc-50)/2;
    if(adc>100 && adc<=200) toReturn=75+(adc-100)/4;
    if(adc>200 && adc<=400) toReturn=100+(adc-200)/8;
    if(adc>400 && adc<=2480) toReturn=125+(adc-400)/16;
    if(adc>2480) toReturn=255;
    return toReturn;
  }

  void ArtBinaryPacketsFromDigis::fillCrvDataPacket(const CrvDigi& digi, CRVHitReadoutPacket& hit, int& globalRocID) 
  {
    //TODO: This is a temporary implementation.
    //There will be a major change on the barIndex+SiPMNumber system,
    //which will be replaced by a channel ID system
    int crvSiPMNumber = digi.GetSiPMNumber();
    mu2e::CRSScintillatorBarIndex crvBarIndex = digi.GetScintillatorBarIndex();
    //Only a toy model is used here. The real implementation will follow.
    int channel = (crvBarIndex.asUint()*4 + crvSiPMNumber) % 64;  //channel within an FEB
    int FEB     = (crvBarIndex.asUint()*4 + crvSiPMNumber) / 64;  //globale FEBId
    uint16_t SiPMID = (FEB<<7) | channel;
    globalRocID = FEB / 24; //global ROCId

    hit.SiPMID = SiPMID;
    hit.HitTime = digi.GetStartTDC();
    hit.NumSamples = 8;
    hit.WaveformSample0 = compressCrvDigi(digi.GetADCs().at(0));  //TODO: There should be a better way of filling the waveform
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
  void ArtBinaryPacketsFromDigis::fillCrvHeaderPacket(CrvDataPacket& crvData, uint8_t globalRocID, uint64_t eventNum) 
  {
    size_t nHits = crvData.hits.size();

    //--------------
    //DataBlocHeader
    //--------------
    // Word 0
    adc_t nBytes = sizeof(DataBlockHeader) + sizeof(CRVROCStatusPacket) + sizeof(CRVHitReadoutPacket)*nHits;
    while(nBytes % 16 != 0) nBytes++;
    crvData.header.ByteCount = nBytes;
    // Word 1
    crvData.header.Hopcount = 0;//currently unused
    crvData.header.PacketType = 5;//PacketType::Dataheader;

    crvData.header.ROCID = globalRocID % number_of_crv_rocs_per_dtc;    //TODO: Is this correct?
    crvData.header.unused1 = 0;
    crvData.header.SubsystemID = DTCLib::DTC_Subsystem_CRV;
    crvData.header.Valid = 1;
    // Word 2
    crvData.header.PacketCount = (crvData.header.ByteCount - 16) / 16;  //TODO: That's how pcie_linux_kernel_module/dtcInterfaceLib/DTC.cpp
                                                                        //interpretes it, but it seems redundant
    crvData.header.unused2 = 0;
    // Word 3
    uint64_t timestamp = eventNum;  //TODO: Is this correct?
    crvData.header.TimestampLow = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    crvData.header.TimestampMed = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    crvData.header.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    crvData.header.Status = 0; // 0 corresponds to "TimeStamp had valid data"
    crvData.header.FormatVersion = format_version;
    // Word 7
    crvData.header.DTCID = globalRocID / number_of_crv_rocs_per_dtc;
    uint8_t  evbMode = 0;//ask Eric
    crvData.header.EVBMode = evbMode;

    //------------------
    //CRVROCStatusPacket
    //------------------
    // Word 0
    crvData.rocStatus.unused1 = 0;
    crvData.rocStatus.PacketType = 0x06;
    crvData.rocStatus.ControllerID = globalRocID % number_of_crv_rocs_per_dtc;    //TODO: Is this correct?
    // Word 1
    crvData.rocStatus.ControllerEventWordCount = sizeof(CRVROCStatusPacket) + sizeof(CRVHitReadoutPacket)*nHits; //TODO: ArtFragmentReader::GetCRVHitCount() seems to interpret this as byte counter and not as word count
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
    crvData.rocStatus.TriggerCount = nHits;  //TODO: Is this is what is meant by TriggerCount? Why isn't this number used in ArtFragmentReader::GetCRVHitCount()?
    // Word 5
    crvData.rocStatus.unused5 = 0;
    crvData.rocStatus.unused6 = 0;
    // Word 6
    crvData.rocStatus.Errors = 0x0;
    crvData.rocStatus.EventType = 0;  //TODO: How is this defined?
  }

  //--------------------------------------------------------------------------------
  //
  //--------------------------------------------------------------------------------
  void ArtBinaryPacketsFromDigis::fillCrvDMABlocks(raw_data_list_t& dataStream, const crv_data_block_list_t& crvDataBlocks)
  {
    // Loop over all ROCs
    uint8_t currentDTCID=0;
    for(uint8_t globalRocID=0; globalRocID<number_of_crv_rocs; globalRocID++) 
    {
      //Add the current DataBlock to the current SuperBlock
      //curDataBlock.setTimestamp(ts); // Overwrite the timestamp
      const CrvDataPacket &crvData = crvDataBlocks.at(globalRocID);
      fillCrvDataStream(dataStream, crvData);

      if(_diagLevel > 1) 
      {
 	if(globalRocID==0 || currentDTCID != crvData.header.DTCID) 
        {
	  std::cout << "================================================" << std::endl;
	  //std::cout << "\t\tTimestamp: " << ts << std::endl;
	  std::cout << "\t\tDTCID: " << (int)crvData.header.DTCID << std::endl;
	  std::cout << "\t\tSYSID: " << (int)crvData.header.SubsystemID << std::endl;
	  currentDTCID = crvData.header.DTCID;
	}
	if(crvData.header.PacketCount > 0)
        {
	  printHeader(crvData.header);
	  printCrvData(crvData);
	}
      }

    } // End loop over DataBlocks
  }

  //--------------------------------------------------------------------------------
  //  method to fill the datastream with the crv packets
  //--------------------------------------------------------------------------------
  void ArtBinaryPacketsFromDigis::fillCrvDataStream(raw_data_list_t& dataStream, const CrvDataPacket& crvData)
  {
    size_t sz = crvData.header.ByteCount;  //byte count was increased to get full chunks of 16 bytes

    assert(sz < sizeof(mu2e_databuff_t));
    if(dataStream.back().second + sz >= sizeof(mu2e_databuff_t)) 
    {
      closeDataBuffer(dataStream);
    }

    auto pos = dataStream.back().second;
    memcpy(&dataStream.back().first[pos], &crvData.header, sizeof(DataBlockHeader));
    pos += sizeof(DataBlockHeader);
    memcpy(&dataStream.back().first[pos], &crvData.rocStatus, sizeof(CRVROCStatusPacket));
    pos += sizeof(CRVROCStatusPacket);

    uint16_t hitCount = crvData.hits.size();

    for(size_t i = 0; i < hitCount; i++) 
    {
      memcpy(&dataStream.back().first[pos], &(crvData.hits[i]), sizeof(CRVHitReadoutPacket));
      pos += sizeof(CRVHitReadoutPacket);
    }

    while (pos % 16 != 0) pos++;  //move forward to get full chunks of 16 bytes
    
    dataStream.back().second = pos;
  }

}


DEFINE_ART_MODULE(mu2e::ArtBinaryPacketsFromDigis);
