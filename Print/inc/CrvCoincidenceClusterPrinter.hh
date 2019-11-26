//
//  Utility class to print CrvCoincidenceCluster
// 
#ifndef Print_inc_CrvCoincidenceClusterPrinter_hh
#define Print_inc_CrvCoincidenceClusterPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CrvCoincidenceClusterPrinter : public ProductPrinter {
  public:

    CrvCoincidenceClusterPrinter() { }
    CrvCoincidenceClusterPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CrvCoincidenceClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CrvCoincidenceClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CrvCoincidenceClusterCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CrvCoincidenceCluster>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CrvCoincidenceCluster& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
