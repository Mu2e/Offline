//
//  Utility class to print SimParticleTimeMap
// 
#ifndef Print_inc_SimParticleTimeMapPrinter_hh
#define Print_inc_SimParticleTimeMapPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class SimParticleTimeMapPrinter : public ProductPrinter {
  public:

    SimParticleTimeMapPrinter() { }
    SimParticleTimeMapPrinter(const ProductPrinter::Config& conf):
      ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<SimParticleTimeMap>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<SimParticleTimeMap>& handle, 
	       std::ostream& os = std::cout);
    void Print(const SimParticleTimeMap& coll, 
	       std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
