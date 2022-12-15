//
//  Utility class to print STMWaveformDigi
//
#ifndef Print_inc_STMWaveformDigiPrinter_hh
#define Print_inc_STMWaveformDigiPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class STMWaveformDigiPrinter : public ProductPrinter {
 public:
  STMWaveformDigiPrinter() {}
  STMWaveformDigiPrinter(const Config& conf) : ProductPrinter(conf) {}

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<STMWaveformDigiCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<STMWaveformDigiCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const STMWaveformDigiCollection& coll,
             std::ostream& os = std::cout);
  void Print(const art::Ptr<STMWaveformDigi>& ptr, int ind = -1,
             std::ostream& os = std::cout);
  void Print(const mu2e::STMWaveformDigi& obj, int ind = -1,
             std::ostream& os = std::cout);

  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintListHeader(std::ostream& os = std::cout);
};

}  // namespace mu2e
#endif
