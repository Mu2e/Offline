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

  int _diagLevel;

  art::InputTag _fragmentsTag;

}; // CrvDigisFromFragments

CrvDigisFromFragments::CrvDigisFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel(config().diagLevel()),
    _fragmentsTag(config().fragmentsTag())
{
  produces<EventNumber_t>();
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

            auto crvHits = crvFragment.GetCRVHitReadoutPackets(iDataBlock);
            for(auto const& crvHit : crvHits)
            {
              //TODO: This is a temporary implementation.
              int channel = crvHit.SiPMID & 0x7F; // right 7 bits
              int FEB = crvHit.SiPMID >> 7;
              int crvBarIndex = (FEB * 64 + channel) / 4;
              int SiPMNumber = (FEB * 64 + channel) % 4;
              if(crvHit.NumSamples!=8)
              {
                std::cerr<<"Number of samples is not 8!"<<std::endl;
                continue;
              }

              //TODO: This is a temporary implementation.
              std::array<int16_t, 8> adc;
              for(int j = 0; j < 8; j++) adc[j] = decompressCrvDigi(crvHit.WaveformSamples[j].ADC);
              crvDigis->emplace_back(adc, crvHit.HitTime, mu2e::CRSScintillatorBarIndex(crvBarIndex), SiPMNumber);
            } //loop over all crvHits

            if(_diagLevel > 0)
            {
              for(auto const& crvHit : crvHits)
              {
                std::cout << "iCrvBlock/iDataBlock: " << iCrvBlock<<"/"<<iDataBlock << std::endl;
                if(_diagLevel > 1)
                {
                  std::cout << "timestamp: " << header->GetEventWindowTag().GetEventWindowTag(true) << std::endl;
                  std::cout << "SubsystemID: " << (uint16_t)header->GetSubsystemID() << std::endl;
                  std::cout << "DTCID: " << (uint16_t)header->GetID() << std::endl;
                  std::cout << "ROCID: " << (uint16_t)header->GetLinkID() << std::endl;
                  std::cout << "packetCount: " << header->GetPacketCount() << std::endl;
                  std::cout << "EVB mode: " << header->GetEVBMode() << std::endl;
                }

                // TODO: This is a temporary implementation.
                int channel = crvHit.SiPMID & 0x7F; // right 7 bits
                int FEB = crvHit.SiPMID >> 7;
                int crvBarIndex = (FEB * 64 + channel) / 4;
                int SiPMNumber = (FEB * 64 + channel) % 4;
                if(crvHit.NumSamples!=8)
                {
                  std::cerr<<"Number of samples is not 8!"<<std::endl;
                  continue;
                }
                std::cout << "SiPMID "<<crvHit.SiPMID
                          << "   channel "<< channel
                          << "   FEB "<< FEB
                          << "   crvBarIndex "<< crvBarIndex
                          << "   SiPMNumber "<< SiPMNumber << std::endl;
                std::cout << "TDC: " << crvHit.HitTime << std::endl;
                std::cout << "nSamples "<<crvHit.NumSamples<<"  ";
                std::cout << "Waveform: {";
                for (size_t j = 0; j < 8; j++) std::cout << "  " << crvHit.WaveformSamples[j].ADC;
                std::cout << "}" << std::endl;
                std::cout << "Waveform decompressed: {";
                for (size_t j = 0; j < 8; j++) std::cout << "  " << decompressCrvDigi(crvHit.WaveformSamples[j].ADC);
                std::cout << "}" << std::endl;
                std::cout << std::endl;
              } //loop over hits
            } //debug output
          } //end parsing CRV DataBlocks
        } //loop over DataBlocks within CRVFragments
      } //loop over SubsystemDataBlocks for the CRV
    } //fragment type: DTCEVT
  } //loop over fragments

  event.put(std::unique_ptr<EventNumber_t>(new EventNumber_t(event.event())));

  event.put(std::move(crvDigis));

} // produce()

DEFINE_ART_MODULE(CrvDigisFromFragments)
