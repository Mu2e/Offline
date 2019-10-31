//
//  Utility class to print MCTrajectory
// 
#ifndef Print_inc_MCTrajectoryPrinter_hh
#define Print_inc_MCTrajectoryPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class MCTrajectoryPrinter : public ProductPrinter {
  public:


    MCTrajectoryPrinter() { }
    MCTrajectoryPrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<MCTrajectoryCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<MCTrajectoryCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const MCTrajectoryCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<MCTrajectory>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::MCTrajectory& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
