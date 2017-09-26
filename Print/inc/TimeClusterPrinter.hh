//
//  Utility class to print TimeCluster
// 
#ifndef Print_inc_TimeClusterPrinter_hh
#define Print_inc_TimeClusterPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class TimeClusterPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    TimeClusterPrinter() { set( fhicl::ParameterSet() ); }
    TimeClusterPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called TimeClusterPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<TimeClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<TimeClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const TimeClusterCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<TimeCluster>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::TimeCluster& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    vecstr _tags;

  };

}
#endif
