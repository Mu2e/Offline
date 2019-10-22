//
//  Utility class to print ComboHit
// 
#ifndef Print_inc_ComboHitPrinter_hh
#define Print_inc_ComboHitPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class ComboHitPrinter : public ProductPrinter {
  public:

    ComboHitPrinter() { }
    ComboHitPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<ComboHitCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<ComboHitCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const ComboHitCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<ComboHit>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::ComboHit& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
