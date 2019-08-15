//
// Add a TracerProduct to the event and watch what happens.
//
// $Id: MakeTracerProduct00_module.cc,v 1.6 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// Other infrastructure includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"

// C++ includes
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

  MakeTracerProduct00::MakeTracerProduct00(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset}
  {
    produces<TracerProduct>();
  }

  void MakeTracerProduct00::produce(art::Event& event) {


    unique_ptr<TracerProduct> prod(new TracerProduct(100+event.id().event()));
    mf::LogVerbatim("Tracing") << "Before put: " << endl;
    event.put(std::move(prod));
    mf::LogVerbatim("Tracing") << "After put: " << endl;

  } // end MakeTracerProduct00::analyze

}  // end namespace mu2e

using mu2e::MakeTracerProduct00;
DEFINE_ART_MODULE(MakeTracerProduct00)
