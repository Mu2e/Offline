// ======================================================================
//
// Extract CRV global run data from CRVDataDecoders
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <iostream>
#include <string>
#include <memory>

namespace mu2e
{

class CrvGRdataFromArtdaqFragments : public art::EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<int>           diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<std::string>   csvFileName{fhicl::Name("csvFileName"), fhicl::Comment("csv file name")};
    fhicl::Atom<bool>          writeCsvFile{fhicl::Name("writeCsvFile"), fhicl::Comment("write csv file")};
    //for currently wrongly encoded subevent headers
    fhicl::Atom<bool>          useSubsystem0{fhicl::Name("useSubsystem0"), fhicl::Comment("consider subevents encoded with subsystem 0")};
  };

  explicit CrvGRdataFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  ~CrvGRdataFromArtdaqFragments() override;
  void produce(art::Event&) override;

  private:
  int                                      _diagLevel;
  std::string                              _csvFileName;
  bool                                     _writeCsvFile;
  std::ofstream                            _csvFile;
  bool                                     _useSubsystem0;
};

CrvGRdataFromArtdaqFragments::CrvGRdataFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
                                                             art::EDProducer{config},
                                                             _diagLevel(config().diagLevel()),
                                                             _csvFileName(config().csvFileName()),
                                                             _writeCsvFile(config().writeCsvFile()),
                                                             _useSubsystem0(config().useSubsystem0())
{
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
CrvGRdataFromArtdaqFragments::~CrvGRdataFromArtdaqFragments()
{
  if(_writeCsvFile)
  {
    _csvFile.close();
  }
}

void CrvGRdataFromArtdaqFragments::produce(art::Event& event)
{
  art::EventNumber_t eventNumber = event.event();
  if(_diagLevel>1)
  {
    std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
  }

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

  size_t totalSize = 0;
  size_t iFragment = 0;
  for(const auto& fragmentHandle : fragmentHandles)
  {
    if(!fragmentHandle.isValid()) continue;

    for(auto artFragment : *fragmentHandle)
    {
      if(artFragment.type()!=mu2e::FragmentType::DTCEVT) continue;

      mu2e::DTCEventFragment dtcEventFragment(artFragment);
      ++iFragment;

      //collect sub events
      auto dtcSubEvents = dtcEventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
      if(_useSubsystem0)
      {
        auto dtcSubEventsTmp = dtcEventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);  //currently wrongly encoded in the DTC Subevent header
        dtcSubEvents.insert(dtcSubEvents.end(),dtcSubEventsTmp.begin(),dtcSubEventsTmp.end()); //temporarily add fragments with wrongly encoded headers
      }

      //check for errors
      const DTCLib::DTC_Event &dtcEvent =  dtcEventFragment.getData();
      size_t expectedSize = dtcEvent.GetEventByteCount();
      size_t actualSize = sizeof(DTCLib::DTC_EventHeader);
      for(size_t iSubEvent=0; iSubEvent<dtcSubEvents.size(); ++iSubEvent) actualSize+=dtcSubEvents.at(iSubEvent).GetSubEventByteCount();
      if(_diagLevel>1)
      {
        std::cout << "Fragment index: " << iFragment << "      expected event size: " << expectedSize << ", actual event size: " << actualSize << std::endl;
      }
      if(expectedSize!=actualSize)
      {
        std::cerr << "mismatch between expected event size and actual event size!" << std::endl;
      }

      //loop through subevents
      for(size_t iSubEvent=0; auto& dtcSubEvent : dtcSubEvents)
      {
        mu2e::CRVDataDecoder decoder(dtcSubEvent);
        for(size_t iDataBlock = 0; iDataBlock < decoder.block_count(); ++iDataBlock)
        {
          totalSize+=decoder.blockSizeBytes(iDataBlock);

          auto block = decoder.dataAtBlockIndex(iDataBlock);
          if(block == nullptr)
          {
            std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cerr << "Unable to retrieve data block." << std::endl;
            continue;
          }

          auto header = block->GetHeader();
          if(!header->isValid())
          {
            std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cerr << "CRV packet is not valid." << std::endl;
            std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
            continue;
          }

          if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV)
          {
            if(_diagLevel>0)
            {
              std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cout << "CRV packet does not have subsystem ID 2." << std::endl;
              std::cout << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
            }
            continue;
          }

          if(_diagLevel>1)
          {
            std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
            std::cout << "packet count: " << header->GetPacketCount() << std::endl;
          }
          if(header->GetPacketCount()>0)
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

            std::unique_ptr<mu2e::CRVDataDecoder::CRVROCStatusPacket> crvRocHeader = decoder.GetCRVROCStatusPacket(iDataBlock);
            if(crvRocHeader==nullptr)
            {
              std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cerr << "Error retrieving CRV ROC Status Packet" << std::endl;
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
            if(decoder.GetCRVGlobalRunInfo(iDataBlock, globalRunInfo) && _diagLevel>0)
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
            if(decoder.GetCRVGlobalRunPayload(iDataBlock, globalRunPayload) && _diagLevel>0)
            {
              for(size_t iWord=0; iWord<globalRunPayload.size(); ++iWord)
              {
                if(iWord%8==0) std::cout << "**** PACKET " << iWord/8 + 2 <<" **** Payload **************" << std::endl;
                std::cout << std::hex << globalRunPayload.at(iWord) << std::dec << " ";
                if(iWord%8==7) {std::cout<<std::endl; std::cout<<std::endl;}
              }
            }

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
          }     // packet count > 0
        }       // loop over data blocks
      }         // loop over subEvents(=decoders)
    }           // loop over fragments
  }             // loop over fragment handles

  if(_diagLevel>1) std::cout << "Total Size: " << totalSize << std::endl;

  if(_diagLevel>0) std::cout << "=========================================" << std::endl;
}

} //namespace mu2e

using mu2e::CrvGRdataFromArtdaqFragments;
DEFINE_ART_MODULE(CrvGRdataFromArtdaqFragments)
