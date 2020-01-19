//
//  Utility class to print StepPointMC
// 
#ifndef Print_inc_StepPointMCPrinter_hh
#define Print_inc_StepPointMCPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class StepPointMCPrinter : public ProductPrinter {
  public:

    struct Config : public ProductPrinter::Config {
      fhicl::Atom<float> pCut{ fhicl::Name("pCut"), 
	  fhicl::Comment("minimum momentum"),-1};
    };

    StepPointMCPrinter():_pCut(-1) { }
    StepPointMCPrinter(const Config& conf):ProductPrinter(conf) { 
      _pCut = conf.pCut();
    }

    // methods to setup parameters
    // do not print if p is below this cut
    void setPCut(double p) { _pCut = p; }
    double pCut() const { return _pCut; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<StepPointMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<StepPointMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const StepPointMCCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<StepPointMC>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::StepPointMC& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:
    double _pCut;

  };

}
#endif
