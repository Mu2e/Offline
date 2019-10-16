//
//  Utility class to print KalRep
// 
#ifndef Print_inc_KalRepPrinter_hh
#define Print_inc_KalRepPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "Print/inc/TrackSummaryPrinter.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class KalRepPrinter : public ProductPrinter {
  public:

    KalRepPrinter() { }
    KalRepPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<KalRepPtrCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<KalRepPtrCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const KalRepPtrCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<KalRep>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const KalRep& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    TrackSummaryPrinter _tsprinter;

  };

}
#endif
