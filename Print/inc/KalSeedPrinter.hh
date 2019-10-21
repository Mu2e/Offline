//
//  Utility class to print KalSeed
// 
#ifndef Print_inc_KalSeedPrinter_hh
#define Print_inc_KalSeedPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class KalSeedPrinter : public ProductPrinter {
  public:

    KalSeedPrinter() { }
    KalSeedPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<KalSeedCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<KalSeedCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const KalSeedCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<KalSeed>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::KalSeed& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
