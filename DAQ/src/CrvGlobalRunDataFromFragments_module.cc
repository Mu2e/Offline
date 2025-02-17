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
    fhicl::Atom<std::string>   csvFileName{fhicl::Name("csvFileName"), fhicl::Comment("csv file name")};
  };

  explicit CrvGlobalRunDataFromFragments(const art::EDAnalyzer::Table<Config>& config);
  ~CrvGlobalRunDataFromFragments() override;
  void analyze(const art::Event&) override;

  private:
  int                                      _diagLevel;
  art::InputTag                            _CRVDataDecodersTag;
  std::string                              _csvFileName;
  std::ofstream                            _csvFile;

};

CrvGlobalRunDataFromFragments::CrvGlobalRunDataFromFragments(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _diagLevel(config().diagLevel()),
    _CRVDataDecodersTag(config().CRVDataDecodersTag()),
    _csvFileName(config().csvFileName())
{
  time_t rawtime;
  time(&rawtime);

  struct tm *timeinfo = localtime(&rawtime);

  char buffer[80];
  strftime(buffer,80,"%Y%m%d%H%M%S",timeinfo);

  std::filesystem::path path(_csvFileName);
  _csvFile.open(path.stem().string()+"_"+std::string(buffer)+path.extension().string());

  _csvFile << "event#,subEvent#,dataBlock#,";
  _csvFile << "EWT(subEventHeader),packetCount,byteCount,";
  _csvFile << "ROCID,wordCount,triggerCount,EWT(ROCstatus),";
  _csvFile << "#EWTs,#markers,lastEWTs,CRC,PLL,lock,injectionTime,injectionWindow";
  _csvFile << std::endl;
}

CrvGlobalRunDataFromFragments::~CrvGlobalRunDataFromFragments()
{
  _csvFile.close();
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
        _csvFile << event.event() << "," << iSubEvent << "," << iDataBlock << ",";   //event number, sub event, data block

        //from subEvent header
        _csvFile << header->GetEventWindowTag().GetEventWindowTag(true) << ",";      //EWT
        _csvFile << header->GetPacketCount() << ",";                                 //packet count
        _csvFile << header->GetByteCount() << ",";                                   //byte count

        //from ROC status header
        const uint16_t *dataPtr=reinterpret_cast<const uint16_t*>(block->GetData());
        _csvFile << ((*(dataPtr+0)) & 0xF) << ",";                                   //ROC ID
        _csvFile << (*(dataPtr+1)) << ",";                                           //word count
        //_csvFile << std::bitset<24>( (*(dataPtr+2)) + ((*(dataPtr+3))<<16) ) << ","; //active FEB flags
        //_csvFile << std::bitset<24>( (*(dataPtr+2)<<16) + (*(dataPtr+3)) ) << ",";   //active FEB flags
        _csvFile << (*(dataPtr+4)) << ",";                                           //trigger count
        _csvFile << ((uint64_t)(*(dataPtr+5))) + ((uint64_t)(*(dataPtr+6))<<16) + ((uint64_t)(*(dataPtr+7))<<32) << ",";  //EWT

        //from Global Run Info packet
        _csvFile << (*(dataPtr+8+1)) << ",";               //#EWTs
        _csvFile << (*(dataPtr+8+2)) << ",";               //#markers
        _csvFile << (*(dataPtr+8+3)) << ",";               //last EWTs
        _csvFile << (((*(dataPtr+8+4))>>8)&0xFF) << ",";   //CRC
        _csvFile << (((*(dataPtr+8+4))>>4)&0xF) << ",";    //PLL
        _csvFile << ((*(dataPtr+8+4))&0x1) << ",";         //lock
        _csvFile << (*(dataPtr+8+5)) << ",";               //injection time
        _csvFile << (*(dataPtr+8+6));                      //injection window

        _csvFile << std::endl;
      }
    }
  }

}

} //namespace mu2e
using mu2e::CrvGlobalRunDataFromFragments;
DEFINE_ART_MODULE(CrvGlobalRunDataFromFragments)
