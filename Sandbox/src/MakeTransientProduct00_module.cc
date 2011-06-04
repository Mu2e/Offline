//
// Test of producing one type of transient data product.
//
// $Id: MakeTransientProduct00_module.cc,v 1.2 2011/06/04 20:37:27 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/04 20:37:27 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "Sandbox/inc/TransientProduct00Collection.hh"
#include "Sandbox/inc/TracerProduct.hh"

//#include <boost/smart_ptr/scoped_ptr.hpp>

#include <memory>

using namespace std;

namespace mu2e {

  class MakeTransientProduct00 : public art::EDProducer {
  public:

    explicit MakeTransientProduct00(fhicl::ParameterSet const& pset);
    virtual ~MakeTransientProduct00() { }

    void produce( art::Event& e);

  private:

  };

  MakeTransientProduct00::MakeTransientProduct00(fhicl::ParameterSet const& pset){
    produces<TransientProduct00Collection>();
  }

  void MakeTransientProduct00::produce(art::Event& event) {

    art::Handle<StrawHitCollection> hitsHandle;
    event.getByLabel("makeSH",hitsHandle);
    StrawHitCollection const& hits = *hitsHandle;

    auto_ptr<TransientProduct00Collection> prod(new TransientProduct00Collection);
    TransientProduct00Collection& p = *prod;

    for ( size_t i=0; i<hits.size(); ++i ){
      p.push_back( hits.at(i) );
    }

    event.put(prod);

    TracerProduct q;
    TracerProduct r(123);
    TracerProduct t(r);
    q=r;
    boost::scoped_ptr<TracerProduct> pp(new TracerProduct(789));

    auto_ptr<TracerProduct> qq(new TracerProduct(-987));

    // Apparently I cannot make an auto_ptr to a scoped_ptr.???
    //auto_ptr<boost::scoped_ptr<TracerProduct>  > rr( new boost::scoped_ptr<TracerProduct>(TracerProduct(999)));

  } // end MakeTransientProduct00::analyze

}  // end namespace mu2e

using mu2e::MakeTransientProduct00;
DEFINE_ART_MODULE(MakeTransientProduct00);
