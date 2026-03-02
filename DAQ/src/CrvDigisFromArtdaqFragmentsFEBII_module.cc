// ======================================================================
//
// Make CRVDigis from CRVDataDecoders
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <iostream>
#include <string>
#include <memory>

namespace mu2e
{

  class CrvDigisFromArtdaqFragmentsFEBII : public art::EDProducer
  {
    public:
      struct Config
      {
        fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
        //for currently wrongly encoded subevent headers
        fhicl::Atom<bool> produceZS{fhicl::Name("produceZS"), fhicl::Comment("produce NZS digi collection"), true};
        fhicl::Atom<bool> produceNZS{fhicl::Name("produceNZS"), fhicl::Comment("produce NZS digi collection"), true};
      };

      explicit CrvDigisFromArtdaqFragmentsFEBII(const art::EDProducer::Table<Config>& config);
      ~CrvDigisFromArtdaqFragmentsFEBII() override {}
      void produce(art::Event&) override;

    private:
      int                                      _diagLevel;
      bool                                     _produceZS;
      bool                                     _produceNZS;
      mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
  };

  CrvDigisFromArtdaqFragmentsFEBII::CrvDigisFromArtdaqFragmentsFEBII(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel(config().diagLevel()),
    _produceZS(config().produceZS()),
    _produceNZS(config().produceNZS())
    {
      if(_produceZS) produces<mu2e::CrvDigiCollection>();
      if(_produceNZS) produces<mu2e::CrvDigiCollection>("NZS");
      produces<mu2e::CrvDAQerrorCollection>();
      produces<mu2e::CrvStatusCollection>();
    }

  void CrvDigisFromArtdaqFragmentsFEBII::produce(art::Event& event)
  {
    // Collection of CrvDigis for the event
    // Digis belonging to the same channel are grouped together and ordered by timestamp
    // This is needed by the reconstruction so that subsequent digis can be merged
    std::unique_ptr<mu2e::CrvDigiCollection> crvDigis(new mu2e::CrvDigiCollection);
    std::unique_ptr<mu2e::CrvDigiCollection> crvDigisNZS(new mu2e::CrvDigiCollection);
    std::unique_ptr<mu2e::CrvDAQerrorCollection> crvDaqErrors(new mu2e::CrvDAQerrorCollection);
    std::unique_ptr<mu2e::CrvStatusCollection> crvStatus(new mu2e::CrvStatusCollection);

    // Temporary collections for unordered digis
    std::map<int, std::vector<mu2e::CrvDigi>> crvDigisTmp;
    std::map<int, std::vector<mu2e::CrvDigi>> crvDigisNZSTmp;

    auto const& channelMap = _channelMap_h.get(event.id());

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

        if(_diagLevel>3)
        {
          std::cout << std::dec << "Fragment index: " << iFragment << std::endl;
        }

        //collect sub events
        auto dtcSubEvents = dtcEventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
        if(dtcSubEvents.empty())
        {
          if(_diagLevel>3)
          {
            std::cout << "This fragment has no CRV subEvents. Size: " << dtcEventFragment.getData().GetEventByteCount() << std::endl;
            size_t nSubEvents = dtcEventFragment.getData().GetSubEvents().size();
            if(nSubEvents>0)
            {
              std::cout << "It has " << nSubEvents << " subEvents of subsystem ID " << (int)dtcEventFragment.getData().GetSubEvents().at(0).GetSubsystem() << std::endl;
            }
          }
          continue;
        }

        //check for errors
        const DTCLib::DTC_Event &dtcEvent =  dtcEventFragment.getData();
        size_t expectedSize = dtcEvent.GetEventByteCount();
        size_t actualSize = sizeof(DTCLib::DTC_EventHeader);
        for(size_t iSubEvent=0; iSubEvent<dtcSubEvents.size(); ++iSubEvent) actualSize+=dtcSubEvents.at(iSubEvent).GetSubEventByteCount();
        if(_diagLevel>3)
        {
          std::cout << "Number of subEvents: " << dtcSubEvents.size() << " Size: " << expectedSize << std::endl;
        }
        if(expectedSize!=actualSize)
        {
          std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
          std::cout << "Fragment index: " << iFragment << "      expected event size: " << expectedSize << ", actual event size: " << actualSize << "      ";
          std::cerr << "mismatch between expected event size and actual event size!" << std::endl;
          crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::byteCountMismatch,iFragment,0,0,0);
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
              std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
              std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cerr << "Unable to retrieve data block." << std::endl;
              crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::unableToGetDataBlock,iFragment,iSubEvent,iDataBlock,0);
              continue;
            }

            auto header = block->GetHeader();
            if(!header->isValid())
            {
              std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
              std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cerr << "CRV packet is not valid." << std::endl;
              std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
              crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::invalidPacket,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
              continue;
            }

            if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV)
            {
              if(_diagLevel>2)
              {
                std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
                std::cout << "CRV packet does not have subsystem ID 2." << std::endl;
                std::cout << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
              }
              crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::wrongSubsystemID,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
              continue;
            }
            auto subeventHeader = dtcSubEvent.GetHeader();
            crvStatus->emplace_back(*header, *subeventHeader);

            if(_diagLevel>1)
            {
              std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
              std::cout << "packet count (from data header): " << header->GetPacketCount() << std::endl;
            }
            if(header->GetPacketCount()>0)
            {
              auto crvRocHeader = decoder.GetCRVROCStatusPacketFEBII(iDataBlock);
              if(crvRocHeader == nullptr)
              {
                std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
                std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
                std::cerr << "Error retrieving CRV ROC Status Packet" << std::endl;
                crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingStatusPacket,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
                continue;
              }
              crvStatus->back().GetROCHeader().push_back(*crvRocHeader);

              auto crvHits = decoder.GetCRVHitRangeFEBII(iDataBlock);
              if(crvHits.error())
              {
                std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
                std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
                std::cerr << "Error unpacking of CRV Hits" << std::endl;
                decoder.PrintBlockFEBII(iDataBlock);
                crvDaqErrors->emplace_back(mu2e::CrvDAQerrorCode::errorUnpackingCrvHits,iFragment,iSubEvent,iDataBlock,header->GetPacketCount());
                break;
              }
              for(auto const& crvHit : crvHits)
              {
                uint16_t dtcID = header->GetID();
                uint16_t linkID = header->GetLinkID();
                uint16_t rocID = dtcID*CRVId::nROCPerDTC + linkID + 1; //ROC IDs are between 1 and 18
                uint16_t rocPort = crvHit.getPortNumber(); //Port numbers beween 1 and 24
                uint16_t febChannel = (crvHit.getFpgaNumber()<<4) + (crvHit.getFpgaChannel() & 0xF);  //use only 4 lowest bits of the fpgaChannel
                //the 5th bit indicates special situations
                //e.g. fake pulses
                if((crvHit.getFpgaChannel() & 0x10) != 0) continue;  //special situation, if the 5th bit of the fpgaChannel is non-zero
                //don't decode them, since there is no match to any offline channel.
                if(rocPort==0) //corrupted data
                {
                  std::cout << std::dec << "Run/Subrun/Event: " << event.run() << "/" << event.subRun() << "/" << eventNumber << std::endl;
                  std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
                  std::cerr << "ROC-port-0 error!" << std::endl;
                  decoder.PrintBlockFEBII(iDataBlock);
                  //TODO: Add to crvDaqErrors
                  continue;
                }

                mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

                uint16_t offlineChannel = channelMap.offline(onlineChannel);
                int crvBarIndex = offlineChannel / CRVId::nChanPerBar;
                int SiPMNumber = offlineChannel % CRVId::nChanPerBar;

                //time stamps coming from FEB-II are in units of 6.25ns, but only the even time stamps are used.
                //convert them into time stamps in units of 12.5ns - same period as the ADC samples.
                bool oddTimestamp=(crvHit.getHitTime()%2==1);

                //extract waveform only when needed (lazy evaluation)
                //auto waveform = crvHit.getWaveform();
                if(_produceZS)
                {
                  crvDigisTmp[offlineChannel].emplace_back(crvHit.getWaveform(), crvHit.getHitTime()/2, false, oddTimestamp,
                      mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber,
                      rocID, rocPort, febChannel);
                }

                //NZS digis - only if requested
                //this is a temporary implementation (simply using ZS data as NZS data) until we get real NZS data
                if(_produceNZS)
                {
                  crvDigisNZSTmp[offlineChannel].emplace_back(crvHit.getWaveform(), crvHit.getHitTime()/2, true, oddTimestamp,
                      mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber,
                      rocID, rocPort, febChannel);
                }
              } // loop over all crvHits

              if(_diagLevel>2)
              {
                std::cout << "EventWindowTag (data header): " << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
                std::cout << "EventWindowTag (ROC header): " << crvRocHeader->GetEventWindowTag() << std::endl;
                std::cout << "SubsystemID (data header): " << (uint16_t)header->GetSubsystemID() << std::endl;
                std::cout << "DTCID (data header): " << (uint16_t)header->GetID() << std::endl;
                std::cout << "LinkID (data header): " << (uint16_t)header->GetLinkID() << std::endl;
                std::cout << "EVB mode (data header): " << (uint16_t)header->GetEVBMode() << std::endl;
                std::cout << "TriggerCount (ROC header): " << crvRocHeader->TriggerCount << std::endl;
                std::cout << "ActiveFEBFlags (ROC header): " << crvRocHeader->GetActiveFEBFlags() << std::endl;
                std::cout << "MicroBunchStatus (ROC header): 0x" << std::hex << (uint32_t)crvRocHeader->GetMicroBunchStatus() << std::dec << std::endl;
              }

              if(_diagLevel>3)
              {
                for(auto const& crvHit : crvHits)
                {
                  // Access hit info directly from raw memory (no copy)
                  uint16_t rocPort = crvHit.getPortNumber(); //Port numbers beween 1 and 24
                  uint16_t febChannel = (crvHit.getFpgaNumber()<<4) + (crvHit.getFpgaChannel() & 0xF);  //use only 4 lowest bits of the fpgaChannel
                  //the 5th bit indicates special situations
                  //e.g. fake pulses
                  if(rocPort==0) continue; //corrupted data

                  uint16_t dtcID = header->GetID();
                  uint16_t linkID = header->GetLinkID();
                  uint16_t rocID = dtcID*CRVId::nROCPerDTC + linkID + 1; //ROC IDs are between 1 and 18

                  if((crvHit.getFpgaChannel() & 0x10) == 0)  //special situation, if the 5th bit of the fpgaChannel is non-zero (see below)
                  {
                    mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

                    uint16_t offlineChannel = channelMap.offline(onlineChannel);
                    int crvBarIndex = offlineChannel / CRVId::nChanPerBar;
                    int SiPMNumber = offlineChannel % CRVId::nChanPerBar;

                    std::cout << "ROCID (increased by 1 to match the Online/Offline-Channel Map) " << rocID
                      << "   rocPort " << rocPort << "   febChannel " << febChannel
                      << "   crvBarIndex " << crvBarIndex << "   SiPMNumber " << SiPMNumber
                      << std::endl;
                  }
                  else
                  {
                    std::cout << "Special hits (fake pulses, etc.) "
                      << "   ROCID (increased by 1 to match the Online/Offline-Channel Map) " << rocID
                      << "   rocPort " << rocPort << "   fpgaNumber " << crvHit.getFpgaNumber() << "   fpgaChannel " << crvHit.getFpgaChannel()
                      << std::endl;
                  }
                  std::cout << "TDC (units of 6.25ns): " << crvHit.getHitTime() << std::endl;

                  // Extract waveform only when needed for diagnostics
                  auto waveform = crvHit.getWaveform();
                  std::cout << "nSamples " << waveform.size() << "  ";
                  std::cout << "Waveform: {";
                  for(size_t i = 0; i < waveform.size(); ++i) std::cout << "  " << waveform.at(i);
                  std::cout << "}" << std::endl;
                  std::cout << std::endl;
                } // loop over hits
              }   // debug output
            }     // packet count > 0
          }       // loop over data blocks
        }         // loop over subEvents(=decoders)
      }           // loop over fragments
    }             // loop over fragment handles

    if(_diagLevel>1) std::cout << "Total Size: " << totalSize << std::endl;

    // Order and store ZS digis
    if(_produceZS)
    {
      for(auto digis=crvDigisTmp.begin(); digis!=crvDigisTmp.end(); ++digis)
      {
        std::sort(digis->second.begin(),digis->second.end(), [](const mu2e::CrvDigi &d1, const mu2e::CrvDigi &d2) {return d1.GetStartTDC()<d2.GetStartTDC();});  //sort by TDC
        crvDigis->insert(crvDigis->end(),digis->second.begin(),digis->second.end());
      }
    }

    // Order and store NZS digis
    if(_produceNZS)
    {
      for(auto digis=crvDigisNZSTmp.begin(); digis!=crvDigisNZSTmp.end(); ++digis)
      {
        std::sort(digis->second.begin(),digis->second.end(), [](const mu2e::CrvDigi &d1, const mu2e::CrvDigi &d2) {return d1.GetStartTDC()<d2.GetStartTDC();});  //sort by TDC
        crvDigisNZS->insert(crvDigisNZS->end(),digis->second.begin(),digis->second.end());
      }
    }

    // Store the crv digis in the event
    if(_produceZS) event.put(std::move(crvDigis));
    if(_produceNZS) event.put(std::move(crvDigisNZS),"NZS");
    event.put(std::move(crvDaqErrors));
    event.put(std::move(crvStatus));
  }

} //namespace mu2e

using mu2e::CrvDigisFromArtdaqFragmentsFEBII;
DEFINE_ART_MODULE(CrvDigisFromArtdaqFragmentsFEBII)
