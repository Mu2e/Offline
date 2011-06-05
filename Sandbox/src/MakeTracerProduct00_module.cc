//
// Test of producing another type of transient data product.
//
// $Id: MakeTracerProduct00_module.cc,v 1.1 2011/06/05 16:41:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/05 16:41:14 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"

#include <memory>

using namespace std;

namespace mu2e {

  class MakeTracerProduct00 : public art::EDProducer {
  public:

    explicit MakeTracerProduct00(fhicl::ParameterSet const& pset);
    virtual ~MakeTracerProduct00() { }

    void produce( art::Event& e);

  private:

  };

  MakeTracerProduct00::MakeTracerProduct00(fhicl::ParameterSet const& pset){
    produces<TracerProduct>();
  }

  void MakeTracerProduct00::produce(art::Event& event) {


    auto_ptr<TracerProduct> prod(new TracerProduct(100+event.id().event()));
    mf::LogVerbatim("Tracing") << "Before put: " << endl;
    event.put(prod);
    mf::LogVerbatim("Tracing") << "After put: " << endl;

    /*
    TracerProduct q;
    TracerProduct r(123);
    TracerProduct t(r);
    q=r;
    boost::scoped_ptr<TracerProduct> pp(new TracerProduct(789));

    auto_ptr<TracerProduct> qq(new TracerProduct(-987));
    */

    // Apparently I cannot make an auto_ptr to a scoped_ptr.???
    //auto_ptr<boost::scoped_ptr<TracerProduct>  > rr( new boost::scoped_ptr<TracerProduct>(TracerProduct(999)));

  } // end MakeTracerProduct00::analyze

}  // end namespace mu2e

using mu2e::MakeTracerProduct00;
DEFINE_ART_MODULE(MakeTracerProduct00);
