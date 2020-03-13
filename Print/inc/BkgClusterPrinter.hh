//
//  Utility class to print BkgCluster
// 
#ifndef Print_inc_BkgClusterPrinter_hh
#define Print_inc_BkgClusterPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgClusterHit.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class BkgClusterPrinter : public ProductPrinter {
  public:

    BkgClusterPrinter() { }
    BkgClusterPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<BkgClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<BkgClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const BkgClusterCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<BkgCluster>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::BkgCluster& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
