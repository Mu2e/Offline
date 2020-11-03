//
//  Utility class to print StrawDigiADCWaveform
// 
#ifndef Print_inc_StrawDigiADCWaveformPrinter_hh
#define Print_inc_StrawDigiADCWaveformPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class StrawDigiADCWaveformPrinter : public ProductPrinter {
  public:

    StrawDigiADCWaveformPrinter() { }
    StrawDigiADCWaveformPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<StrawDigiADCWaveformCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<StrawDigiADCWaveformCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const StrawDigiADCWaveformCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<StrawDigiADCWaveform>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::StrawDigiADCWaveform& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
