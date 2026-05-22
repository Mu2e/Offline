// ======================================================================
//
// STMPrintFragments_plugin:  Add CRV data products to the event
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
     fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("stmTag for new file")};
  };

  // --- C'tor/d'tor:
  explicit STMPrintFragments(const art::EDAnalyzer::Table<Config>& config);

  // --- Production:
  virtual void analyze(const Event&);

  private:

  art::InputTag _stmFragmentsTag;
}; // STMPrintFragments

// ======================================================================

STMPrintFragments::STMPrintFragments(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}
    ,_stmFragmentsTag(config().stmTag())
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
  for (auto& frag : *STMContainerFragments) {
    ++frag_counter;
    
    //New lines
    artdaq::ContainerFragment contf(frag); // interpret the fragment as a ContainerFragemnt (Will look inside here)
    std::cout<<"N Blocks in the container = " << contf.block_count() << std::endl; //Should be 3 for the 3 STM Fragments
  
    for (size_t ii = 0; ii< contf.block_count(); ++ii){
      auto inner = contf.at(ii);
      const auto dataBegin = inner->dataBegin();
      const auto dataEnd = inner->dataEnd();
      auto frag_id = inner->fragmentID();
      const auto stmDataBegin = reinterpret_cast<int16_t const*>(dataBegin);
      const auto stmDataEnd = reinterpret_cast<int16_t const*>(dataEnd);
      std::cout << "Frag_ID = " << frag_id << std::endl;
      std::cout << "Container block_count = "<<contf.block_count()<<std::endl;
      
      for (auto i = stmDataBegin; i != stmDataEnd; ++i) {
	std::cout << "Frag #" << frag_counter << ", inner #" << ii  << ": *(stmDataBegin+" << i - stmDataBegin << ") = " << *i << std::endl;
      }

    }
  }
}// produce()

// ======================================================================

DEFINE_ART_MODULE(STMPrintFragments)
