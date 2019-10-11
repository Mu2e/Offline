//
//  Utility class to print CaloShowerStep
// 
#ifndef Print_inc_CaloShowerStepPrinter_hh
#define Print_inc_CaloShowerStepPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloShowerStepPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    CaloShowerStepPrinter():_eCut(-1) { }
    CaloShowerStepPrinter(const ConfigE& conf):ProductPrinter(conf) { 
      _eCut = conf.eCut();
    }

    // do not print if e is below this cut
    void setECut(double e) { _eCut = e; }
    double eCut() const { return _eCut; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CaloShowerStepCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CaloShowerStepCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CaloShowerStepCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CaloShowerStep>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CaloShowerStep& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    double _eCut;

  };

}
#endif
