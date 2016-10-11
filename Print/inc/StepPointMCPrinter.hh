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

    typedef std::vector<std::string> vecstr;

    StepPointMCPrinter() { set( fhicl::ParameterSet() ); }
    StepPointMCPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }
    // usually customized in the subclass for items relevant to that product

    // methods to setup parameters
    // do not print if p is below this cut
    void setPCut(double p) { _pCut = p; }

    // pset should contain a table called StepPointMCPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

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
    vecstr _tags;

  };

}
#endif
