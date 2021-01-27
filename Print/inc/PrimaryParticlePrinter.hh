
//
//  Utility class to print PrimaryParticle
// 
#ifndef Print_inc_PrimaryParticlePrinter_hh
#define Print_inc_PrimaryParticlePrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "Print/inc/GenParticlePrinter.hh"
#include "Print/inc/SimParticlePrinter.hh"
#include "art/Framework/Principal/Handle.h"

namespace mu2e {

  class PrimaryParticlePrinter : public ProductPrinter {
  public:

    PrimaryParticlePrinter() {}
    PrimaryParticlePrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<PrimaryParticle>& handle, 
               std::ostream& os = std::cout);
    void Print(const art::ValidHandle<PrimaryParticle>& handle, 
               std::ostream& os = std::cout);
    void Print(const mu2e::PrimaryParticle& obj, 
	       int ind = -1, std::ostream& os = std::cout);
    void PrintHeader(const std::string& tag, 
                     std::ostream& os = std::cout);
  private:
    GenParticlePrinter gprint_;
    SimParticlePrinter sprint_;
  };

}
#endif
