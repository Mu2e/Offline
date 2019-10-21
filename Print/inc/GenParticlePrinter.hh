//
//  Utility class to print GenParticle
// 
#ifndef Print_inc_GenParticlePrinter_hh
#define Print_inc_GenParticlePrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class GenParticlePrinter : public ProductPrinter {
  public:

    GenParticlePrinter() {}
    GenParticlePrinter(const Config& conf):ProductPrinter(conf) {}

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<GenParticleCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<GenParticleCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const GenParticleCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<GenParticle>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::GenParticle& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
