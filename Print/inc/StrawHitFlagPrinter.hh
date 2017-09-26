//
//  Utility class to print StrawHitFlag
// 
#ifndef Print_inc_StrawHitFlagPrinter_hh
#define Print_inc_StrawHitFlagPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class StrawHitFlagPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    StrawHitFlagPrinter() { set( fhicl::ParameterSet() ); }
    StrawHitFlagPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called StrawHitFlagPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<StrawHitFlagCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<StrawHitFlagCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const StrawHitFlagCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const mu2e::StrawHitFlag& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    vecstr _tags;

  };

}
#endif
