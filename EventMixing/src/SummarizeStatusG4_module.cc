//
// Extract all StatusG4 objects from the event, compute the overall status and
// put it into the event.  Status G4 objects can be found in two places:
//  - as a StatusG4 object that is a top level data product.
//  - within a MixingSummary object.
//
// $Id: SummarizeStatusG4_module.cc,v 1.6 2013/03/15 18:20:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 18:20:22 $
//
// Contact person Rob Kutschke
//
// Notes:
// 1) To do:
//    Can we do better to avoid cross-stitching of the summary.  What happens
//    if there is more than one bare StatusG4 object?  Can we pick the wrong one?
//    Is there a loophole that would cause us to double count anything?
//    Does is make sense to generate a mixing summary object for a simple input file?
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/MixingSummary.hh"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Handle.h"

#include "cetlib_except/exception.h"
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
    art::EDProducer{pset},
    // Run time parameters
    diagLevel_(pset.get<int>("diagLevel",0)),
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel","")){

    produces<StatusG4>();
  }

  void SummarizeStatusG4::produce( art::Event& event) {

    static bool firstEvent(true);

    // Default construct the output object.
    unique_ptr<StatusG4> summaryStatus( new StatusG4() );

    // Printout on some events.
    unique_ptr<mf::LogInfo> log(nullptr);
    if (  firstEvent || diagLevel_ > 0 ) {
      log = unique_ptr<mf::LogInfo>(new mf::LogInfo("MIXING"));
      (*log) << "Creating Summary of StatusG4 from the following products: \n";
    }

    // Incorporate at most one top level StatusG4 object.
    // If the name is present, it is an error for the StatusG4 object to be absent
    // and the handle will throw on dereference.
    if ( g4ModuleLabel_.size() > 0 ){
      art::Handle<StatusG4> g4StatusHandle;
      event.getByLabel( g4ModuleLabel_, g4StatusHandle);

      summaryStatus->add(*g4StatusHandle);
      if ( log.get() != 0 ) (*log)<< "   " << g4StatusHandle.provenance()->branchName() << "\n";
    }

    // Get all of the tracker MixingSummary collections from the event:
    art::ProductInstanceNameSelector selector("");

    typedef std::vector< art::Handle<MixingSummary> > HandleVector;
    HandleVector statHandles;
    event.getMany( selector, statHandles);

    // Incorporate each of the summary StatusG4 objects from each MixingSummary.
    for ( HandleVector::const_iterator i=statHandles.begin(), e=statHandles.end();
          i != e; ++i ){

      if ( log.get() != 0 ) (*log)<< "   " << i->provenance()->branchName() << "\n";
      summaryStatus->add((**i).status());
    }

    event.put(std::move(summaryStatus));
    firstEvent = false;

  } // end produce

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SummarizeStatusG4);
