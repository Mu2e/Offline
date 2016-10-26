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

    typedef std::vector<std::string> vecstr;

    MCTrajectoryPrinter() { set( fhicl::ParameterSet() ); }
    MCTrajectoryPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }
    // usually customized in the subclass for items relevant to that product

    // pset should contain a table called MCTrajectoryPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

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

  private:
    vecstr _tags;

  };

}
#endif
