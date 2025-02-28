// CrvGlobalRunDataFromFragments:  Print out CRV global run data

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <string>
#include <memory>

namespace mu2e
{

class CrvGlobalRunDataFromFragments : public art::EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<int>           diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"), fhicl::Comment("crv Fragments Tag")};
    fhicl::Atom<std::string>   csvFileName{fhicl::Name("csvFileName"), fhicl::Comment("csv file name")};
    fhicl::Atom<bool>          writeCsvFile{fhicl::Name("writeCsvFile"), fhicl::Comment("write csv file")};
  };

  explicit CrvGlobalRunDataFromFragments(const art::EDProducer::Table<Config>& config);
  ~CrvGlobalRunDataFromFragments() override;
  void produce(art::Event&) override;

  private:
  int                                      _diagLevel;
  art::InputTag                            _CRVDataDecodersTag;
  std::string                              _csvFileName;
  bool                                     _writeCsvFile;
  std::ofstream                            _csvFile;

};

CrvGlobalRunDataFromFragments::CrvGlobalRunDataFromFragments(const art::EDProducer::Table<Config>& config) : 
	                                                     art::EDProducer{config}, 
							      _diagLevel(config().diagLevel()),
							     _CRVDataDecodersTag(config().CRVDataDecodersTag()),
							     _csvFileName(config().csvFileName()),
							     _writeCsvFile(config().writeCsvFile())
{
  produces<mu2e::CRVDataDecoder::CRVGlobalRunDataCollection>();
  produces<mu2e::CrvDAQerrorCollection>();

  if(_writeCsvFile)
  {
    _csvFile.open(_csvFileName);

    _csvFile << "event#,subEvent#,dataBlock#,";
    _csvFile << "EWT(subEventHeader),packetCount,byteCount(subEventHeader),";
    _csvFile << "ROCID,wordCount(ROCstatus),triggerCount,EWT(ROCstatus),";
    _csvFile << "#EWTs,#markers,lastEWTs,CRC,PLL,lock,injectionTime,injectionWindow";
    _csvFile << std::endl;
  }
}

CrvGlobalRunDataFromFragments::~CrvGlobalRunDataFromFragments()
{
  if(_writeCsvFile)
  {
    _csvFile.close();
  }
}

void CrvGlobalRunDataFromFragments::produce(art::Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  auto CRVDataDecoders = event.getValidHandle<std::vector<mu2e::CRVDataDecoder> >(_CRVDataDecodersTag);
  size_t nSubEvents = CRVDataDecoders->size();

  if(_diagLevel>0) std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;

  size_t totalSize = 0;
  for(const auto &frag : *CRVDataDecoders)
  {
    for(size_t i = 0; i < frag.block_count(); ++i)
    {
      auto size = frag.blockSizeBytes(i);
      totalSize += size;
    }
  }

  if(_diagLevel>1)
  {
    std::cout << "#SubEvents: " << nSubEvents << std::endl;
    std::cout << "Total Size: " << totalSize << " bytes." << std::endl;
  }

  std::unique_ptr<mu2e::CRVDataDecoder::CRVGlobalRunDataCollection> crv_globalRun(new mu2e::CRVDataDecoder::CRVGlobalRunDataCollection);
  std::unique_ptr<mu2e::CrvDAQerrorCollection> crv_daq_errors(new mu2e::CrvDAQerrorCollection);

  for(size_t iSubEvent = 0; iSubEvent < nSubEvents; ++iSubEvent)
  {
    const mu2e::CRVDataDecoder& CRVDataDecoder((*CRVDataDecoders)[iSubEvent]);
    CRVDataDecoder.setup_event();

    for(size_t iDataBlock = 0; iDataBlock < CRVDataDecoder.block_count(); ++iDataBlock)
    {
      auto block = CRVDataDecoder.dataAtBlockIndex(iDataBlock);
      if(block == nullptr)
      {
        std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cerr << "Unable to retrieve data block." << std::endl;
        std::cerr << std::endl;
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::unableToGetDataBlock,iSubEvent,iDataBlock,0);
        continue;
      }
      auto header = block->GetHeader();

      if(!header->isValid())
      {
        std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cerr << "CRV packet is not valid." << std::endl;
        std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        std::cerr << std::endl;
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::invalidPacket,iSubEvent,iDataBlock,header->GetPacketCount());
        continue;
      }

      if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV)
      {
        if(_diagLevel>2)
        {
          std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
          std::cout << "CRV packet does not have subsystem ID 2." << std::endl;
          std::cout << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
          std::cout << std::endl;
        }
        crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::wrongSubsystemID,iSubEvent,iDataBlock,header->GetPacketCount());
        continue;
      }

      if(_diagLevel>0)
      {
        std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cout << "packet count: " << header->GetPacketCount() << std::endl;
        std::cout << std::endl;
      }

      if(header->GetPacketCount() > 0)
      {
	if(_diagLevel>0)
        {
          std::cout << "********** DTC Data Header **************" << std::endl;
          std::cout << "ValidFlag: " << (uint16_t)header->isValid() << std::endl;
          std::cout << "Status: " << (uint16_t)header->GetStatus() << std::endl;
          std::cout << "EventWindowTag (TDC Header): " << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
          std::cout << "SubsystemID: " << (uint16_t)header->GetSubsystemID() << std::endl;
          std::cout << "Subsystem: " << (uint16_t)header->GetSubsystem() << std::endl;
          std::cout << "DTCID: " << (uint16_t)header->GetID() << std::endl;
          std::cout << "ROCID: " << (uint16_t)header->GetLinkID() << std::endl;
          std::cout << "HopCount: " << (uint16_t)header->GetHopCount() << std::endl;
          std::cout << "PacketType: " << (uint16_t)header->GetPacketType() << std::endl;
          std::cout << "PacketCount: " << header->GetPacketCount() << std::endl;
          std::cout << "ByteCount: " << header->GetByteCount() << std::endl;
          std::cout << "EVB mode: " << (uint16_t)header->GetEVBMode() << std::endl;
          std::cout << "Version: " << (uint16_t)header->GetVersion() << std::endl;
          std::cout << std::endl;
        }

        std::unique_ptr<mu2e::CRVDataDecoder::CRVROCStatusPacket> crvRocHeader = CRVDataDecoder.GetCRVROCStatusPacket(iDataBlock);
        if(crvRocHeader==nullptr)
        {
          std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
          std::cerr << "Error retrieving CRV ROC Status Packet" << std::endl;
          crv_daq_errors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingStatusPacket,iSubEvent,iDataBlock,header->GetPacketCount());
          continue;
        }

        if(crvRocHeader!=nullptr && _diagLevel>0)
        {
          std::cout << "**** PACKET 0 **** ROC Status Header ****" << std::endl;
          std::cout << "ROCID (ROC Status): "<< (uint16_t)crvRocHeader->ControllerID << std::endl;
          std::cout << "WordCount: "<< crvRocHeader->ControllerEventWordCount << std::endl;
          std::cout << "ActiveFEBFlags: " << crvRocHeader->GetActiveFEBFlags() << std::endl;
          std::cout << "TriggerCount: " << crvRocHeader->TriggerCount << std::endl;
          std::cout << "EventWindowTag (ROC Status): " << crvRocHeader->GetEventWindowTagGlobalRun() << std::endl;
          std::cout<<std::endl;
        }

        mu2e::CRVDataDecoder::CRVGlobalRunInfo globalRunInfo;
        if(CRVDataDecoder.GetCRVGlobalRunInfo(iDataBlock, globalRunInfo) && _diagLevel>0)
        {
          std::cout << "**** PACKET 1 **** Global Run Info ******" << std::endl;
          std::cout << "Word 0: "<< std::hex << globalRunInfo.word0 << std::dec << std::endl;
          std::cout << "#EWT: "<< globalRunInfo.EWTCount << std::endl;
          std::cout << "#Markers: "<< globalRunInfo.markerCount << std::endl;
          std::cout << "LastEWTs: "<< globalRunInfo.lastEWT << std::endl;
          std::cout << "CRC: " << globalRunInfo.CRC << std::endl;
          std::cout << "PLL: " << globalRunInfo.PLL << std::endl;
          std::cout << "Lock: " << globalRunInfo.lock << std::endl;
          std::cout << "InjectionTime: "<< globalRunInfo.injectionTime << std::endl;
          std::cout << "InjectionWindow: "<< globalRunInfo.injectionWindow << std::endl;
          std::cout << "Word 7: "<< std::hex << globalRunInfo.word7 << std::dec << std::endl;
          std::cout<<std::endl;
        }

        mu2e::CRVDataDecoder::CRVGlobalRunPayload globalRunPayload;
        if(CRVDataDecoder.GetCRVGlobalRunPayload(iDataBlock, globalRunPayload) && _diagLevel>0)
        {
          for(size_t iWord=0; iWord<globalRunPayload.size(); ++iWord)
          {
            if(iWord%8==0) std::cout << "**** PACKET " << iWord/8 + 2 <<" **** Payload **************" << std::endl;
            std::cout << std::hex << globalRunPayload.at(iWord) << std::dec << " ";
            if(iWord%8==7) {std::cout<<std::endl; std::cout<<std::endl;}
          }
        }

        crv_globalRun->emplace_back(*crvRocHeader,globalRunInfo,globalRunPayload);

        if(_writeCsvFile)
	{
          _csvFile << event.event() << "," << iSubEvent << "," << iDataBlock << ","; //event number, sub event, data block

          //from subEvent header
          _csvFile << header->GetEventWindowTag().GetEventWindowTag(true) << ",";    //EWT
          _csvFile << header->GetPacketCount() << ",";                               //packet count
          _csvFile << header->GetByteCount() << ",";                                 //byte count

          //from ROC status header
          _csvFile << (uint16_t)crvRocHeader->ControllerID << ",";                   //ROC ID
          _csvFile << crvRocHeader->ControllerEventWordCount << ",";                 //word count
          _csvFile << crvRocHeader->TriggerCount << ",";                             //trigger count
          _csvFile << crvRocHeader->GetEventWindowTagGlobalRun() << ",";             //EWT

          //from Global Run Info packet
          _csvFile << globalRunInfo.EWTCount << ",";                                 //#EWTs
          _csvFile << globalRunInfo.markerCount << ",";                              //#markers
          _csvFile << globalRunInfo.lastEWT << ",";                                  //last EWTs
          _csvFile << globalRunInfo.CRC << ",";                                      //CRC
          _csvFile << globalRunInfo.PLL << ",";                                      //PLL
          _csvFile << globalRunInfo.lock << ",";                                     //lock
          _csvFile << globalRunInfo.injectionTime << ",";                            //injection time
          _csvFile << globalRunInfo.injectionWindow;                                 //injection window

          _csvFile << std::endl;
	}   // write CSV file

      }     // end parsing CRV DataBlocks
    }       // loop over DataBlocks within CRVDataDecoders
  }         // Close loop over fragments

  event.put(std::move(crv_globalRun));
  event.put(std::move(crv_daq_errors));

  if(_diagLevel>0) std::cout << "=========================================" << std::endl;

} //produce

} //namespace mu2e

using mu2e::CrvGlobalRunDataFromFragments;
DEFINE_ART_MODULE(CrvGlobalRunDataFromFragments)
