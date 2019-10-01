//
//  Utility class to print TrackClusterMatch
// 
#ifndef Print_inc_TrackClusterMatchPrinter_hh
#define Print_inc_TrackClusterMatchPrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class TrackClusterMatchPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    TrackClusterMatchPrinter() { set( fhicl::ParameterSet() ); }
    TrackClusterMatchPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called TrackClusterMatchPrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<TrackClusterMatchCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<TrackClusterMatchCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const TrackClusterMatchCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<TrackClusterMatch>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::TrackClusterMatch& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    vecstr _tags;

  };

}
#endif
