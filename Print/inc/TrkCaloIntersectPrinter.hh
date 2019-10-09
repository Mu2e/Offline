//
//  Utility class to print TrkCaloIntersect
// 
#ifndef Print_inc_TrkCaloIntersectPrinter_hh
#define Print_inc_TrkCaloIntersectPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class TrkCaloIntersectPrinter : public ProductPrinter {
  public:

    TrkCaloIntersectPrinter() { }
    TrkCaloIntersectPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<TrkCaloIntersectCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<TrkCaloIntersectCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const TrkCaloIntersectCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<TrkCaloIntersect>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::TrkCaloIntersect& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
