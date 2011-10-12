//
// Extract all StatusG4 objects from the event, compute the overall status and
// put it into the event.
//
// $Id: SummarizeStatusG4_module.cc,v 1.1 2011/10/12 20:09:27 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/10/12 20:09:27 $
//
// Contact person Rob Kutschke
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/MixingSummary.hh"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/Selector.h"
#include "art/Persistency/Common/Handle.h"

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <string>
#include <memory>

using namespace std;

namespace mu2e {

  class SummarizeStatusG4 : public art::EDProducer {
  public:

    explicit SummarizeStatusG4(fhicl::ParameterSet const& pset);
    virtual ~SummarizeStatusG4() { }

    virtual void produce(art::Event& e);

  private:

    // Start: run time parameters

    // Diagnostics printout level
    int diagLevel_;

    // Module label of the g4 module that made the hits.
    std::string g4ModuleLabel_;

    // End: run time parameters

  };

  SummarizeStatusG4::SummarizeStatusG4(fhicl::ParameterSet const& pset) :

    // Run time parameters
    diagLevel_(pset.get<int>("diagLevel",0)),
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel","")){

    produces<StatusG4>();
  }

  void SummarizeStatusG4::produce( art::Event& event) {

    static bool firstEvent(true);

    // Default construct the output object.
    auto_ptr<StatusG4> summaryStatus( new StatusG4() );

    auto_ptr<mf::LogInfo> log(0);
    if (  firstEvent || diagLevel_ > 0 ) {
      log = auto_ptr<mf::LogInfo>(new mf::LogInfo("MIXING"));
      (*log) << "Creating Summary of StatusG4 from the following products: \n";
    }

    // Incorporate at the completion status of G4, if present.
    if ( g4ModuleLabel_.size() > 0 ){
      art::Handle<StatusG4> g4StatusHandle;
      event.getByLabel( g4ModuleLabel_, g4StatusHandle);

      summaryStatus->add(*g4StatusHandle);
      if ( firstEvent || diagLevel_ > 0 ){
        if ( log.get() != 0 ) (*log)<< "   " << g4StatusHandle.provenance()->branchName() << "\n";
      }
    }

    // Get all of the tracker MixingSummary collections from the event:
    art::ProductInstanceNameSelector selector("");

    typedef std::vector< art::Handle<MixingSummary> > HandleVector;
    HandleVector statHandles;

    event.getMany( selector, statHandles);
    for ( HandleVector::const_iterator i=statHandles.begin(), e=statHandles.end();
          i != e; ++i ){

      if ( log.get() != 0 ) (*log)<< "   " << i->provenance()->branchName() << "\n";
      MixingSummary const& sum(**i);
      summaryStatus->add(sum.status());
    }
    event.put(summaryStatus);

    firstEvent = false;

  } // end produce

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SummarizeStatusG4);
