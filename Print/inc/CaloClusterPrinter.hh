//
//  Utility class to print CaloCluster
// 
#ifndef Print_inc_CaloClusterPrinter_hh
#define Print_inc_CaloClusterPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CaloClusterPrinter : public ProductPrinter {
  public:

    typedef std::vector<std::string> vecstr;

    CaloClusterPrinter() { set( fhicl::ParameterSet() ); }
    CaloClusterPrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }

    // pset should contain a table called CaloClusterPrinter
    void set(const fhicl::ParameterSet& pset);

    // do not print if p is below this cut
    void setECut(double e) { _eCut = e; }

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CaloClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CaloClusterCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CaloClusterCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<CaloCluster>& ptr, 
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::CaloCluster& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:

    double _eCut;
    vecstr _tags;

  };

}
#endif
