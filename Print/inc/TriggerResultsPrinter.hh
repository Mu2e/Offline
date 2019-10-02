//
//  Utility class to print TriggerResults
// 
#ifndef Print_inc_TriggerResultsPrinter_hh
#define Print_inc_TriggerResultsPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class TriggerResultsPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    TriggerResultsPrinter() { set( fhicl::ParameterSet() ); }
    TriggerResultsPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called TriggerResultsPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<art::TriggerResults>& handle, 
               std::ostream& os = std::cout);
    void Print(const art::ValidHandle<art::TriggerResults>& handle, 
               std::ostream& os = std::cout);
    void Print(const art::TriggerResults& obj, 
	       int ind = -1, std::ostream& os = std::cout);
    void PrintHeader(const std::string& tag, 
                     std::ostream& os = std::cout);

  private:

    vecstr _tags;

  };

}
#endif
