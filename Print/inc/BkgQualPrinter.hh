//
//  Utility class to print BkgQual
// 
#ifndef Print_inc_BkgQualPrinter_hh
#define Print_inc_BkgQualPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/BkgQual.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class BkgQualPrinter : public ProductPrinter {
  public:

    BkgQualPrinter() { }
    BkgQualPrinter(const ConfigE& conf):ProductPrinter(conf) { 
      _eCut = conf.eCut();
    }

    // do not print if e is below this cut
    void setECut(double e) { _eCut = e; }
    double eCut() const { return _eCut; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<BkgQualCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<BkgQualCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const BkgQualCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<BkgQual>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::BkgQual& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    double _eCut;

  };

}
#endif
