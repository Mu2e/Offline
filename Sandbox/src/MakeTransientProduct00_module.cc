//
// Test of producing one type of transient data product.
//
// $Id: MakeTransientProduct00_module.cc,v 1.7 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "Sandbox/inc/TransientProduct00Collection.hh"

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

  MakeTransientProduct00::MakeTransientProduct00(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset}
  {
    produces<TransientProduct00Collection>();
  }

  void MakeTransientProduct00::produce(art::Event& event) {

    art::Handle<StrawHitCollection> hitsHandle;
    event.getByLabel("makeSH",hitsHandle);
    StrawHitCollection const& hits = *hitsHandle;

    unique_ptr<TransientProduct00Collection> prod(new TransientProduct00Collection);
    TransientProduct00Collection& p = *prod;

    // The transient product holds pointers to all of the StrawHits
    for ( size_t i=0; i<hits.size(); ++i ){
      p.push_back( hits.at(i) );
    }

    event.put(std::move(prod));

  } // end MakeTransientProduct00::analyze

}  // end namespace mu2e

using mu2e::MakeTransientProduct00;
DEFINE_ART_MODULE(MakeTransientProduct00)
