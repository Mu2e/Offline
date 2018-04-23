//
//  Utility class to print CrvDigiMC
// 
#ifndef Print_inc_CrvDigiMCPrinter_hh
#define Print_inc_CrvDigiMCPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CrvDigiMCPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    CrvDigiMCPrinter() { set( fhicl::ParameterSet() ); }
    CrvDigiMCPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called CrvDigiMCPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CrvDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CrvDigiMCCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CrvDigiMCCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CrvDigiMC>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CrvDigiMC& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    vecstr _tags;

  };

}
#endif
