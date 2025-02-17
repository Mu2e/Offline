// ======================================================================
//
// CrvGlobalRunDataFromFragments_plugin:  Print out CRV global run data
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <ctime>
#include <filesystem>

namespace mu2e
{

class CrvGlobalRunDataFromFragments : public art::EDAnalyzer
{
  public:
  struct Config
  {
    fhicl::Atom<int>           diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> CRVDataDecodersTag{fhicl::Name("crvTag"), fhicl::Comment("crv Fragments Tag")};
    fhicl::Atom<std::string>   cvsFileName{fhicl::Name("cvsFileName"), fhicl::Comment("cvs file name")};
  };

  explicit CrvGlobalRunDataFromFragments(const art::EDAnalyzer::Table<Config>& config);
  ~CrvGlobalRunDataFromFragments() override;
  void analyze(const art::Event&) override;

  private:
  int                                      _diagLevel;
  art::InputTag                            _CRVDataDecodersTag;
  std::string                              _cvsFileName;
  std::ofstream                            _cvsFile;

};

CrvGlobalRunDataFromFragments::CrvGlobalRunDataFromFragments(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _diagLevel(config().diagLevel()),
    _CRVDataDecodersTag(config().CRVDataDecodersTag()),
    _cvsFileName(config().cvsFileName())
{
  time_t rawtime;
  time(&rawtime);

  struct tm *timeinfo = localtime(&rawtime);

  char buffer[80];
  strftime(buffer,80,"%Y%m%d%H%M%S",timeinfo);

  std::filesystem::path path(_cvsFileName);
  _cvsFile.open(path.stem().string()+"_"+std::string(buffer)+path.extension().string());

  _cvsFile << "event#,subEvent#,dataBlock#,";
  _cvsFile << "EWT(subEventHeader),packetCount,byteCount,";
  _cvsFile << "ROCID,wordCount,triggerCount,EWT(ROCstatus),";
  _cvsFile << "#EWTs,#markers,lastEWTs,CRC,PLL,lock,injectionTime,injectionWindow";
  _cvsFile << std::endl;
}

CrvGlobalRunDataFromFragments::~CrvGlobalRunDataFromFragments()
{
  _cvsFile.close();
}

void CrvGlobalRunDataFromFragments::analyze(const art::Event& event)
{
  auto CRVDataDecoders = event.getValidHandle<std::vector<mu2e::CRVDataDecoder> >(_CRVDataDecodersTag);

  for(size_t iSubEvent = 0; iSubEvent < CRVDataDecoders->size(); ++iSubEvent)
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
        continue;
      }
      auto header = block->GetHeader();

      if(!header->isValid())
      {
        std::cerr << "iSubEvent/iDataBlock: " << iSubEvent << "/" << iDataBlock << std::endl;
        std::cerr << "CRV packet is not valid." << std::endl;
        std::cerr << "sub system ID: "<<(uint16_t)header->GetSubsystemID()<<" packet count: "<<header->GetPacketCount() << std::endl;
        std::cerr << std::endl;
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
        continue;
      }

      if(header->GetPacketCount() > 1)
      {
        _cvsFile << event.event() << "," << iSubEvent << "," << iDataBlock << ",";   //event number, sub event, data block

        //from subEvent header
        _cvsFile << header->GetEventWindowTag().GetEventWindowTag(true) << ",";      //EWT
        _cvsFile << header->GetPacketCount() << ",";                                 //packet count
        _cvsFile << header->GetByteCount() << ",";                                   //byte count

        //from ROC status header
        const uint16_t *dataPtr=reinterpret_cast<const uint16_t*>(block->GetData());
        _cvsFile << ((*(dataPtr+0)) & 0xF) << ",";                                   //ROC ID
        _cvsFile << (*(dataPtr+1)) << ",";                                           //word count
        //_cvsFile << std::bitset<24>( (*(dataPtr+2)) + ((*(dataPtr+3))<<16) ) << ","; //active FEB flags
        //_cvsFile << std::bitset<24>( (*(dataPtr+2)<<16) + (*(dataPtr+3)) ) << ",";   //active FEB flags
        _cvsFile << (*(dataPtr+4)) << ",";                                           //trigger count
        _cvsFile << ((uint64_t)(*(dataPtr+5))) + ((uint64_t)(*(dataPtr+6))<<16) + ((uint64_t)(*(dataPtr+7))<<32) << ",";  //EWT

        //from Global Run Info packet
        _cvsFile << (*(dataPtr+8+1)) << ",";               //#EWTs
        _cvsFile << (*(dataPtr+8+2)) << ",";               //#markers
        _cvsFile << (*(dataPtr+8+3)) << ",";               //last EWTs
        _cvsFile << (((*(dataPtr+8+4))>>8)&0xFF) << ",";   //CRC
        _cvsFile << (((*(dataPtr+8+4))>>4)&0xF) << ",";    //PLL
        _cvsFile << ((*(dataPtr+8+4))&0x1) << ",";         //lock
        _cvsFile << (*(dataPtr+8+5)) << ",";               //injection time
        _cvsFile << (*(dataPtr+8+6));                      //injection window

        _cvsFile << std::endl;
      }
    }
  }

}

} //namespace mu2e
using mu2e::CrvGlobalRunDataFromFragments;
DEFINE_ART_MODULE(CrvGlobalRunDataFromFragments)
