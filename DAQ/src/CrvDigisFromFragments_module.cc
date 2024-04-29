// ======================================================================
//
// CrvDigisFromFragments_plugin:  Add CRV data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class CrvDigisFromFragments;
}

using art::CrvDigisFromFragments;

// ======================================================================

class art::CrvDigisFromFragments : public EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"),
                                               fhicl::Comment("crv Fragments Tag")};
  };

  // --- C'tor/d'tor:
  explicit CrvDigisFromFragments(const art::EDProducer::Table<Config>& config);
  ~CrvDigisFromFragments() override {}

  // --- Production:
  void produce(Event&) override;

  private:
  int                                      _diagLevel;
  art::InputTag                            _CRVDataDecodersTag;
  mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;

}; // CrvDigisFromFragments

// ======================================================================

CrvDigisFromFragments::CrvDigisFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _diagLevel(config().diagLevel()), _CRVDataDecodersTag(config().CRVDataDecodersTag())
{
  produces<mu2e::CrvDigiCollection>();
}

// ----------------------------------------------------------------------

void CrvDigisFromFragments::produce(Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  auto CRVDataDecoders = event.getValidHandle<std::vector<mu2e::CRVDataDecoder> >(_CRVDataDecodersTag);
  size_t nSubEvents = CRVDataDecoders->size();

  if(_diagLevel>1)
  {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << nSubEvents << " CRV SubEvents." << std::endl;

    size_t totalSize = 0;
    for(const auto &frag : *CRVDataDecoders)
    {
      for(size_t i = 0; i < frag.block_count(); ++i)
      {
        auto size = frag.blockSizeBytes(i);
        totalSize += size;
      }
    }

    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;
  }

  // Collection of CrvDigis for the event
  std::unique_ptr<mu2e::CrvDigiCollection> crv_digis(new mu2e::CrvDigiCollection);
  auto const& channelMap = _channelMap_h.get(event.id());

  // Loop over the CRV fragments
  for(size_t iSubEvent = 0; iSubEvent < nSubEvents; ++iSubEvent)
  {
    const mu2e::CRVDataDecoder& CRVDataDecoder((*CRVDataDecoders)[iSubEvent]);
    CRVDataDecoder.setup_event();

    for(size_t iDataBlock = 0; iDataBlock < CRVDataDecoder.block_count(); ++iDataBlock)
    {
      if(_diagLevel>0) std::cout << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;

      auto block = CRVDataDecoder.dataAtBlockIndex(iDataBlock);
      if(block == nullptr)
      {
        std::cerr << "Unable to retrieve block in " << std::endl;
        continue;
      }
      auto header = block->GetHeader();
/*
FIXME: This function will be available in a new release of artdaq_core_mu2e
      if(!header->isValid())
      {
        std::cerr << "CRV packet is not valid." << std::endl;
        std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        continue;
      }
*/
      if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV)
      {
        std::cerr << "CRV packet does not have system ID 2." << std::endl;
        std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        continue;
      }

      if(_diagLevel>0) std::cout << "packet count: " << header->GetPacketCount() << std::endl;
      if(header->GetPacketCount() > 0)
      {
        auto crvRocHeader = CRVDataDecoder.GetCRVROCStatusPacket(iDataBlock);
        if(crvRocHeader == nullptr)
        {
          std::cerr << "Error retrieving CRV ROC Status Packet from DataBlock in " << iDataBlock << std::endl;
          continue;
        }

        auto crvHits = CRVDataDecoder.GetCRVHits(iDataBlock);
        for(auto const& crvHit : crvHits)
        {
          const auto& crvHitInfo = crvHit.first;
          const auto& waveform = crvHit.second;

          uint16_t rocID = crvHitInfo.controllerNumber + 1; // FIXME ROC IDs between 1 and 17   //also: header->GetLinkID()+1;
          uint16_t rocPort = crvHitInfo.portNumber;
          uint16_t febChannel = crvHitInfo.febChannel;
          mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

          uint16_t offlineChannel = channelMap.offline(onlineChannel);
          int crvBarIndex = offlineChannel / 4;
          int SiPMNumber = offlineChannel % 4;

          for(size_t i = 0; i < crvHitInfo.NumSamples; i += mu2e::CrvDigi::NSamples)
          {
            std::array<int16_t, mu2e::CrvDigi::NSamples> adc = {0};
            for(size_t j = 0; j < mu2e::CrvDigi::NSamples && i+j < crvHitInfo.NumSamples; ++j)
              adc[j] = waveform.at(i+j).ADC;

            // CrvDigis use a constant array size of 8 samples
            // waveforms with more than 8 samples need to be written to multiple CrvDigis
            // the TDC increases by 8 for every subsequent CrvDigi
            crv_digis->emplace_back(adc, crvHitInfo.HitTime + i, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);
          }
        } // loop over all crvHits

        if(_diagLevel>1)
        {
          std::cout << "EventWindowTag (TDC header): "
                    << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
          std::cout << "SubsystemID: " << (uint16_t)header->GetSubsystemID() << std::endl;
          std::cout << "DTCID: " << (uint16_t)header->GetID() << std::endl;
          std::cout << "ROCID: " << (uint16_t)header->GetLinkID() << std::endl;
          std::cout << "packetCount: " << header->GetPacketCount() << std::endl;
          std::cout << "EVB mode: " << (uint16_t)header->GetEVBMode() << std::endl;
          std::cout << "TriggerCount: " << crvRocHeader->TriggerCount << std::endl;
          std::cout << "ActiveFEBFlags: " << crvRocHeader->GetActiveFEBFlags() << std::endl;
          std::cout << "ROCID (ROC header): " << (uint16_t)crvRocHeader->ControllerID  << std::endl;
          std::cout << "EventWindowTag (ROC header): " << crvRocHeader->GetEventWindowTag() << std::endl;
        }

        if(_diagLevel>0)
        {
          for(auto const& crvHit : crvHits)
          {
            const auto& crvHitInfo = crvHit.first;
            const auto& waveform = crvHit.second;

            uint16_t rocID = crvHitInfo.controllerNumber + 1; // FIXME  //ROC IDs are between 1 and 17
            uint16_t rocPort = crvHitInfo.portNumber;
            uint16_t febChannel = crvHitInfo.febChannel;
            mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

            uint16_t offlineChannel = channelMap.offline(onlineChannel);
            int crvBarIndex = offlineChannel / 4;
            int SiPMNumber = offlineChannel % 4;

            std::cout << "ROCID (increased by 1 to match the Online/Offline-Channel Map) " << rocID
                      << "   rocPort " << rocPort << "   febChannel " << febChannel
                      << "   crvBarIndex " << crvBarIndex << "   SiPMNumber " << SiPMNumber
                      << std::endl;
            std::cout << "TDC: " << crvHitInfo.HitTime << std::endl;
            std::cout << "nSamples " << crvHitInfo.NumSamples << "  ";
            std::cout << "Waveform: {";
            for(size_t i = 0; i < crvHitInfo.NumSamples; i++)
              std::cout << "  " << waveform.at(i).ADC;
            std::cout << "}" << std::endl;
            std::cout << std::endl;
          } // loop over hits
        }   // debug output
      }     // end parsing CRV DataBlocks
    }       // loop over DataBlocks within CRVDataDecoders
  }         // Close loop over fragments

  // Store the straw digis and calo digis in the event
  event.put(std::move(crv_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(CrvDigisFromFragments)
