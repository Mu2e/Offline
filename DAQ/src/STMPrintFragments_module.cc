// ======================================================================
//
// STMPrintFragments_plugin: Print STM art file payload and data
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>
#include "artdaq-core/Data/ContainerFragment.hh"

#include <iostream>
#include <string>
#include <memory>

namespace art
{
  class STMPrintFragments;
}

using art::STMPrintFragments;

// ======================================================================

class art::STMPrintFragments : public EDAnalyzer
{
  public:
  struct Config
   {
     fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"),
       fhicl::Comment("stmTag for new file")};
     fhicl::Atom<bool> printEverything {fhicl::Name("printEverything"),
       fhicl::Comment("Switch to print everything in STM File ")};
     fhicl::Atom<bool> printPayloads {fhicl::Name("printPayloads"),
       fhicl::Comment("Switch to print only fragment payloads ")};
     fhicl::Atom<bool> printInnerFrags {fhicl::Name("printInnerFrags"),
       fhicl::Comment("Switch to print inner frags begining from header to end of payload ")};
  };

  // --- C'tor/d'tor:
  explicit STMPrintFragments(const art::EDAnalyzer::Table<Config>& config);

  // --- Production:
  virtual void analyze(const Event&);

  private:

  art::InputTag _stmFragmentsTag;
  bool _printEverything{false};
  bool _printPayloads{false};
  bool _printInnerFrags{false};
}; // STMPrintFragments

// ======================================================================

STMPrintFragments::STMPrintFragments(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}
    ,_stmFragmentsTag(config().stmTag())
    ,_printEverything(config().printEverything())
    ,_printPayloads(config().printPayloads())
    ,_printInnerFrags(config().printInnerFrags())
{}

// ----------------------------------------------------------------------

void STMPrintFragments::analyze(const Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  auto STMContainerFragments = event.getValidHandle<artdaq::Fragments>(_stmFragmentsTag); //Changed from STMFragments -> STMContainerFragments

  std::cout << std::dec << "Analyzer: Run " << event.run() << ", subrun " << event.subRun()
            << ", event " << eventNumber << " has " << std::endl;
  std::cout << STMContainerFragments->size() << " STM fragments." << std::endl;

  int frag_counter = 0;
  for (const auto& frag : *STMContainerFragments) {
    ++frag_counter;
    artdaq::ContainerFragment cont_frag(frag);
    std::cout << "N blocks in the container = " << cont_frag.block_count() << std::endl;

    if (_printEverything){//Print Eveyrthing Block-----------------------------------------
      mu2e::STMFragment container_stm_frag(frag);

      const auto dataBegin = frag.dataBegin();
      const auto dataEnd = frag.dataEnd();
      const auto containerDataBegin = reinterpret_cast<int16_t const*>(dataBegin);
      const auto containerDataEnd = reinterpret_cast<int16_t const*>(dataEnd);
      auto frag_id = frag.fragmentID();

      std::cout << "\n==== FULL CONTAINER DUMP ====\n"
                << "Container Frag_ID      : " << frag_id << "\n"
                << "Container Frag #       : " << frag_counter << "\n"
                << std::boolalpha
                << "HPGe Container         : " << container_stm_frag.isHPGeContainer() << "\n"
                << "LaBr Container         : " << container_stm_frag.isLaBrContainer() << "\n"
                << "Container block_count  : " << cont_frag.block_count() << std::endl;

      for (auto ii = containerDataBegin ; ii !=containerDataEnd; ++ii){
        std::cout << "Container frag # " << frag_counter << ": *(containerDataBegin+" << ii - containerDataBegin << ") = " << *ii << std::endl;
      }


    } //Print Everything Option -----------------------------------------

    if (_printInnerFrags){ //Print Inner Fragment Block----------------------------------------
      for (size_t i = 0; i < cont_frag.block_count(); ++i){
        auto inner_frag = cont_frag.at(i);
        mu2e::STMFragment stm_frag(*inner_frag);
        const auto frag_id = inner_frag->fragmentID();
        const auto stmDataBegin = stm_frag.dataBegin();
        const auto stmDataWords = stm_frag.dataWords();
        const auto stmDataEnd = stmDataBegin + stmDataWords;

        const auto stmPayloadBegin = stm_frag.payloadBegin();
        const auto stmPayloadWords = stm_frag.payloadWords();
        const auto physicalWords = inner_frag->dataSizeBytes() / sizeof(int16_t);

        //Always print
        std::cout << "Frag_ID = " << frag_id << std::endl;
        std::cout << "Container block count = " << cont_frag.block_count() << std::endl;
        std::cout << "stmdataWords = " << stmDataWords
                  << ", physicalWords = " << physicalWords
                  << ", payloadWords = " << stmPayloadWords << std::endl;
        if ( stm_frag.isRaw() ){
          std::cout << "Following information comes from Raw Header:\n" << std::boolalpha
                    << ", hasValidAnchors = " << stm_frag.hasValidAnchors()
                    << ", rawLength  = " << stm_frag.rawLength()
                    << ", badData = " << stm_frag.badData()
                    << ", missingData = " << stm_frag.missing()
                    << ", zsPrescale = " << stm_frag.zsPrescaled()
                    << ", rawPrescale = " << stm_frag.rawPrescaled()
                    << ", zsPrescaleValue = " << unsigned(stm_frag.zsPrescaleValue())
                    << ", rawPrescaleValue = " << unsigned(stm_frag.rawPrescaleValue())
                    << "\n";
        }

        for (auto ii = stmDataBegin; ii != stmDataEnd; ++ii){
          std::cout << "Frag #" << frag_counter
                    << ", inner #" << i
                    << ": *(stmDataBegin+" << ii - stmDataBegin << ") = "
                    << *ii << std::endl;
        }//should print the whole data set

        size_t physicalPayloadWords{0};
        if ( stm_frag.isRaw() ){
          physicalPayloadWords = physicalWords > stm::RawHeader::WORDS ? physicalWords - stm::RawHeader::WORDS : 0;
          size_t nPrintLimit = std::min<size_t>(physicalPayloadWords,20);
          std::cout << "\n";

        //Print preview of the payload to use as cross reference
        //with dataset printed out for raw frags
          for (auto ii = stmPayloadBegin; ii != stmPayloadBegin + nPrintLimit; ++ii){
            std::cout << "Frag #" << frag_counter
                      << ", inner #" << i
                      << ": *(stmPayloadBegin+" << ii - stmPayloadBegin << ")="
                      << *ii << std::endl;
          }
        }

      }//inner frag loop
    }// DataSets Block Option-------------------------------------------------------
    
    if (_printPayloads){
      for (size_t i = 0; i < cont_frag.block_count(); ++i){
	auto inner_frag = cont_frag.at(i);
	mu2e::STMFragment stm_frag(*inner_frag);
	size_t physicalPayloadWords{0};

	if (stm_frag.isRaw()) {
	  const auto stmDataWords = stm_frag.dataWords();

	  const auto stmPayloadBegin = stm_frag.payloadBegin();
	  const auto physicalWords = inner_frag->dataSizeBytes() / sizeof(int16_t);
	  physicalPayloadWords = physicalWords > stm::RawHeader::WORDS ? physicalWords - stm::RawHeader::WORDS : 0;
	  
	  std::cout << "\n=== Raw Payload Dump ===\n"
		    << "Frag # : " << frag_counter << "\n"
		    << ", inner # : " << i << "\n"
		    << ", PhysicalWords : " << physicalWords << "\n"
		    << ", stmDataWords : " << stmDataWords << "\n"
		    << ", rawLength from Header : " << stm_frag.rawLength() << std::endl;
	  for (auto ii = stmPayloadBegin ; ii != stmPayloadBegin + physicalPayloadWords ; ++ii){
	    std::cout << "Raw payload[" << ii - stmPayloadBegin << "] = " << *ii << "\n"; 
	  }
	  
	}
	else if (stm_frag.isZS()){
          const auto stmDataWords = stm_frag.dataWords();
	  const auto physicalWords = inner_frag->dataSizeBytes() / sizeof(int16_t);
	  
	  auto dataPtr = stm_frag.dataBegin();
          auto dataWords = stm_frag.dataWords();
          auto dataEnd = dataPtr + dataWords;

	  std::cout << "\n=== ZS Payload Dump ===\n"
		    << "Frag # : " << frag_counter << "\n"
                    << ", inner # : " << i << "\n"
                    << ", PhysicalWords : " << physicalWords << "\n"
                    << ", stmDataWords : " << stmDataWords << "\n";
	  size_t seg = 0;
	  size_t totalLen= 0;
	  uint16_t lastZSindex = 0;
	  uint16_t lastLen = 0;
	  
	  while (dataPtr+2 <= dataEnd){
	    uint16_t current_zs_location = static_cast<uint16_t>(dataPtr[0]);
            uint16_t current_zs_size = static_cast<uint16_t>(dataPtr[1]);
	    std::cout << "Currnent ZS Index : " << current_zs_location
		      << ", current ZS Length : " << current_zs_size << "\n" ;
            auto adc = dataPtr + 2;
	    if (adc + current_zs_size > dataEnd){
	      std::cout << "Malformed or truncated ZS Region: Declared size extends beyond fragment data " << "\n";
	      break;
	    }
	    for (size_t j = 0; j < current_zs_size; ++j) {
	      std::cout << "ZS Payload[" << j << "] = " << adc[j] << ", region = " << seg <<  "\n";
	    }
	    
	    std::cout << "Region = " << seg
		      << " , zs_index = " << current_zs_location
		      << " , zs_size = " << current_zs_size << "\n" ;
	    
	    // Update Variables
            ++seg;
            dataPtr = adc + current_zs_size;
	  } // end of while loop
	}
	  
	else if (stm_frag.isPH()){
	  size_t digiWords = stm_frag.payloadWords();
	  auto const* digiPtr = stm_frag.payloadBegin();  
	  std::cout << "\n=== PH Payload Dump ===\n" << "\n";
	  if (digiWords %2 !=0) {
	    std::cout << "Found a malformed PH frag " << "\n";
	  }
	  
	  for (size_t i_PH = 0; i_PH + 1 < digiWords ; i_PH +=2){
	    uint32_t time = static_cast<uint16_t>(digiPtr[i_PH]);
            int16_t PH = digiPtr[i_PH + 1];
	    
	    std::cout << "PH PayloadTime[" << i_PH << "] = " << time << "\n"
		      << "PH PayloadHit[" << i_PH+1 << "] = " << PH << "\n";
	  }
	} 
	else {
	  std::cout << "Unread Inner fragment" << "\n"
		    << "Inner frag i : " << i << "\n"
		    << "Frag ID      : " << inner_frag->fragmentID() << "\n";
	}
	  
      }// End of container frag loop

    }// Payload Block Option------------------------------------------------------

  } //container frag loop
}//analyze
// ======================================================================
DEFINE_ART_MODULE(STMPrintFragments)
