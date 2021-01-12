//
// Add a std::vector<TracerProduct> to the event and watch what happens.
//
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
#include "Sandbox/inc/TracerProductCollection.hh"

// C++ includes
#include <memory>

using namespace std;

namespace mu2e {

  class MakeVectorTracerProduct00 : public art::EDProducer {
  public:

    explicit MakeVectorTracerProduct00(fhicl::ParameterSet const& pset);
    virtual ~MakeVectorTracerProduct00() { }

    void produce( art::Event& e);

  private:

  };

  MakeVectorTracerProduct00::MakeVectorTracerProduct00(fhicl::ParameterSet const& pset): 
    art::EDProducer{pset}
  {
    produces<TracerProductCollection>();
  }

  void MakeVectorTracerProduct00::produce(art::Event& event) {

    mf::LogVerbatim("Tracing") << "Start produce: ";
    unique_ptr<TracerProductCollection > prod(new TracerProductCollection );
    prod->push_back( new TracerProduct( 100*event.id().event() + 1) );
    prod->push_back( new TracerProduct( 100*event.id().event() + 2) );

    mf::LogVerbatim("Tracing") << "Before put: " << endl;
    event.put(std::move(prod));
    mf::LogVerbatim("Tracing") << "After put: " << endl;

  } // end MakeVectorTracerProduct00::analyze

}  // end namespace mu2e

using mu2e::MakeVectorTracerProduct00;
DEFINE_ART_MODULE(MakeVectorTracerProduct00)
