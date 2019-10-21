//
//  Utility class to print SimParticlePtr
// 
#ifndef Print_inc_SimParticlePtrPrinter_hh
#define Print_inc_SimParticlePtrPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class SimParticlePtrPrinter : public ProductPrinter {
  public:


    SimParticlePtrPrinter() { }
    SimParticlePtrPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<SimParticlePtrCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<SimParticlePtrCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const SimParticlePtrCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<SimParticle>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
