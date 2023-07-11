//
//  Utility class to print CaloHitMC
//
#ifndef Print_inc_CaloHitMCPrinter_hh
#define Print_inc_CaloHitMCPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

class CaloHitMCPrinter : public ProductPrinter {
 public:
  CaloHitMCPrinter() {}
  CaloHitMCPrinter(const ConfigE& conf) : ProductPrinter(conf) {
    _eCut = conf.eCut();
  }

  // do not print if p is below this cut
  void setECut(double e) { _eCut = e; }
  double eCut() const { return _eCut; }

  // all the ways to request a printout
  void Print(art::Event const& event, std::ostream& os = std::cout) override;
  void Print(const art::Handle<CaloHitMCCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const art::ValidHandle<CaloHitMCCollection>& handle,
             std::ostream& os = std::cout);
  void Print(const CaloHitMCCollection& coll, std::ostream& os = std::cout);
  void Print(const art::Ptr<CaloHitMC>& ptr, int ind = -1,
             std::ostream& os = std::cout);
  void Print(const mu2e::CaloHitMC& obj, int ind = -1,
             std::ostream& os = std::cout);

  void PrintHeader(const std::string& tag, std::ostream& os = std::cout);
  void PrintListHeader(std::ostream& os = std::cout);

 private:
  double _eCut;
};

}  // namespace mu2e
#endif
