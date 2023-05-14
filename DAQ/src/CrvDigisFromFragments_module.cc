// ======================================================================
//
// CrvDigisFromFragments_plugin:  Add CRV data products to the event
//
// ======================================================================

#include "Offline/DAQ/inc/FragmentType.hh"
#include "Offline/DAQ/inc/DTCEventFragment.hh"
#include "Offline/DAQ/inc/DTC_Packets.h"
#include "Offline/DAQ/inc/CRVFragment.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class CrvDigisFromFragments;
}

using art::CrvDigisFromFragments;

class art::CrvDigisFromFragments : public EDProducer
{

  public:
  struct Config
  {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> fragmentsTag{fhicl::Name("fragmentsTag"),
                                            fhicl::Comment("Fragments Tag")};
  };

  explicit CrvDigisFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CrvDigisFromFragments() {}

  virtual void produce(Event&);

  private:
//  int decompressCrvDigi(uint8_t adc);
  int16_t decompressCrvDigi(int16_t adc);

  int           _diagLevel;
  art::InputTag _fragmentsTag;

  mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;

}; // CrvDigisFromFragments

CrvDigisFromFragments::CrvDigisFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel(config().diagLevel()),
    _fragmentsTag(config().fragmentsTag())
{
  produces<mu2e::CrvDigiCollection>();
}

int16_t CrvDigisFromFragments::decompressCrvDigi(int16_t adc)
{
  //TODO: This is a temporary implementation.
  return adc;
}

void CrvDigisFromFragments::produce(Event& event)
{
  std::unique_ptr<mu2e::CrvDigiCollection> crvDigis(new mu2e::CrvDigiCollection);

  auto const& channelMap = _channelMap_h.get(event.id());

  auto fragments = event.getValidHandle<artdaq::Fragments>(_fragmentsTag);
  for(auto fragment=fragments->begin(); fragment!=fragments->end(); ++fragment)
  {
    if(fragment->type()==mu2e::FragmentType::DTCEVT)
    {
      mu2e::DTCEventFragment dtcEventFragment(*fragment);
      auto dtcEvent = dtcEventFragment.getData();
      dtcEvent.SetupEvent();

      auto crvBlocks = dtcEvent.GetSubsystemDataBlocks(DTCLib::DTC_Subsystem_CRV);  //FIXME: removed check for DTC_Subsystem_CRV in DTC_Packets.h
      for(size_t iCrvBlock=0; iCrvBlock<crvBlocks.size(); ++iCrvBlock)
      {
        mu2e::CRVFragment crvFragment(crvBlocks.at(iCrvBlock).blockPointer, crvBlocks.at(iCrvBlock).byteSize);

        for(size_t iDataBlock = 0; iDataBlock < crvFragment.block_count(); ++iDataBlock)
        {
          auto block = crvFragment.dataAtBlockIndex(iDataBlock);
          if(block == nullptr)
          {
            std::cerr << "Unable to retrieve block " << iDataBlock << "!" << std::endl;
            continue;
          }
          auto header = block->GetHeader();
          if(header->GetSubsystemID() != 2)
          {
            throw cet::exception("DATA") << " CRV packet does not have system ID 2";
          }

          if(header->GetPacketCount() > 0)
          {
            auto crvRocHeader = crvFragment.GetCRVROCStatusPacket(iDataBlock);
            if (crvRocHeader == nullptr)
            {
              std::cerr << "Error retrieving CRV ROC Status Packet from DataBlock " << iDataBlock << std::endl;
              continue;
            }

            auto crvHits = crvFragment.GetCRVHits(iDataBlock);
            for(auto const& crvHit : crvHits)
            {
              const auto &crvHitInfo = crvHit.first;
              const auto &waveform   = crvHit.second;

              uint16_t rocID      = header->GetLinkID()+1;  //FIXME
              uint16_t rocPort    = crvHitInfo.portNumber;
              uint16_t febChannel = crvHitInfo.febChannel;
              mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

              uint16_t offlineChannel = channelMap.offline(onlineChannel);
              int crvBarIndex = offlineChannel / 4;
              int SiPMNumber  = offlineChannel % 4;

              for(int i=0; i<crvHitInfo.NumSamples; i+=8)
              {
                std::array<int16_t, 8> adc={0};
                for(int j=i; j<i+8 && j<crvHitInfo.NumSamples; ++j) adc[j] = decompressCrvDigi(waveform.at(j).ADC);

                //CrvDigis use a constant array size of 8 samples
                //waveforms with more than 8 samples need to be written to multiple CrvDigis
                //the TDC increases by 8 for every subsequent CrvDigi
                crvDigis->emplace_back(adc, crvHitInfo.HitTime+i, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);
              }
            } //loop over all crvHits

            if(_diagLevel > 0)
            {
              for(auto const& crvHit : crvHits)
              {
                std::cout << "iCrvBlock/iDataBlock: " << iCrvBlock<<"/"<<iDataBlock << std::endl;
                if(_diagLevel > 1)
                {
                  std::cout << "EventWindowTag (TDC header): " << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
                  std::cout << "SubsystemID: " << (uint16_t)header->GetSubsystemID() << std::endl;
                  std::cout << "DTCID: " << (uint16_t)header->GetID() << std::endl;
                  std::cout << "ROCID: " << (uint16_t)header->GetLinkID() << std::endl;
                  std::cout << "packetCount: " << header->GetPacketCount() << std::endl;
                  std::cout << "EVB mode: " << header->GetEVBMode() << std::endl;
                  std::cout << "TriggerCount: " << crvRocHeader->TriggerCount << std::endl;
                  std::cout << "ActiveFEBFlags: " << crvRocHeader->GetActiveFEBFlags() << std::endl;
                  std::cout << "EventWindowTag (ROC header): " << crvRocHeader->GetEventWindowTag() << std::endl;
                }

                const auto &crvHitInfo = crvHit.first;
                const auto &waveform   = crvHit.second;

                uint16_t rocID      = header->GetLinkID()+1; //FIXME  //ROC IDs are between 1 and 17
                uint16_t rocPort    = crvHitInfo.portNumber;
                uint16_t febChannel = crvHitInfo.febChannel;
                mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

                uint16_t offlineChannel = channelMap.offline(onlineChannel);
                int crvBarIndex = offlineChannel / 4;
                int SiPMNumber  = offlineChannel % 4;

                std::cout << "ROCID "<<rocID
                          << "   rocPort "<<rocPort
                          << "   febChannel "<< febChannel
                          << "   crvBarIndex "<< crvBarIndex
                          << "   SiPMNumber "<< SiPMNumber << std::endl;
                std::cout << "TDC: " << crvHitInfo.HitTime << std::endl;
                std::cout << "nSamples "<<crvHitInfo.NumSamples<<"  ";
                std::cout << "Waveform: {";
                for (size_t i = 0; i < crvHitInfo.NumSamples; i++) std::cout << "  " << waveform.at(i).ADC;
                std::cout << "}" << std::endl;
                std::cout << "Waveform decompressed: {";
                for (size_t i = 0; i < crvHitInfo.NumSamples; i++) std::cout << "  " << decompressCrvDigi(waveform.at(i).ADC);
                std::cout << "}" << std::endl;
                std::cout << std::endl;
              } //loop over hits
            } //debug output
          } //end parsing CRV DataBlocks
        } //loop over DataBlocks within CRVFragments
      } //loop over SubsystemDataBlocks for the CRV
    } //fragment type: DTCEVT
  } //loop over fragments

  event.put(std::move(crvDigis));

} // produce()

DEFINE_ART_MODULE(CrvDigisFromFragments)
