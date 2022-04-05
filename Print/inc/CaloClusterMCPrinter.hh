//
//  Utility class to print CaloClusterMC
//
#ifndef Print_inc_CaloClusterMCPrinter_hh
#define Print_inc_CaloClusterMCPrinter_hh

#include <cstring>
#include <iostream>

#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloClusterMCPrinter : public ProductPrinter {
  public:

    CaloClusterMCPrinter() { }
    CaloClusterMCPrinter(const ConfigE& conf):ProductPrinter(conf) {
      _eCut = conf.eCut();
    }

    // do not print if e is below this cut
    void setECut(double e) { _eCut = e; }
    double eCut() const { return _eCut; }

    // all the ways to request a printout
    void Print(art::Event const& event,
               std::ostream& os = std::cout) override;
    void Print(const art::Handle<CaloClusterMCCollection>& handle,
               std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CaloClusterMCCollection>& handle,
               std::ostream& os = std::cout);
    void Print(const CaloClusterMCCollection& coll,
               std::ostream& os = std::cout);
    void Print(const art::Ptr<CaloClusterMC>& ptr,
               int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CaloClusterMC& obj,
               int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag,
                     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    double _eCut;

  };

}
#endif
