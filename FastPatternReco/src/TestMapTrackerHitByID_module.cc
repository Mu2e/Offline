//
// example of a module to read Data of the Electrons tracks that came from the targets
//
// $Id: TestMapTrackerHitByID_module.cc,v 1.2 2013/10/21 21:01:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:22 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <map>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
//#include "TrackerGeom/inc/Tracker.hh"
//#include "TrackerGeom/inc/Straw.hh"
//#include "ITrackerGeom/inc/Cell.hh"
//#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
//#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"

using namespace std;

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  class TestMapTrackerHitByID : public art::EDAnalyzer {
  public:

    explicit TestMapTrackerHitByID(fhicl::ParameterSet const& pset);
    virtual ~TestMapTrackerHitByID() {}

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event const& e);

  private:

    // Start: run time parameters

//    // The module label of this module.
//    std::string _moduleLabel;
//
//    // Label of the G4 module
//    std::string _g4ModuleLabel;
//
//    // Name of the tracker StepPoint collection
//    std::string _trackerStepPoints;
//
//    // Label of the module that made the hits.
//    std::string _makerModuleLabel;
//
//    // Label of the generator.
//    std::string _generatorModuleLabel;

    // Label of the module that made the hits.
    std::string _mapTrackerHitByID;

  };

  TestMapTrackerHitByID::TestMapTrackerHitByID(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _mapTrackerHitByID(pset.get<string>("mapTrackerHitByID"))
{
}

  void TestMapTrackerHitByID::beginJob(){

	  cout<<"Starting Test jos!"<<endl;
  }

  void TestMapTrackerHitByID::analyze(art::Event const& event ) {


//    const Tracker& tracker = getTrackerOrThrow();
//    const TTracker &ttr = static_cast<const TTracker&>( tracker );
//    const std::vector<Device> ttrdev = ttr.getDevices();

//    art::Handle<StrawHitCollection> pdataHandle;
//    event.getByLabel(_makerModuleLabel,pdataHandle);
//    StrawHitCollection const* hits = pdataHandle.product();


    art::Handle<TrackerHitByID> hitByIDHandle;
    event.getByLabel(_mapTrackerHitByID,hitByIDHandle);
    TrackerHitByID const* hitByID = hitByIDHandle.product();

    TrackerHitByID::const_iterator hitByID_it;

    cout<<"In Run "<<event.run()<<" event "<<event.event()<<endl;

    for ( hitByID_it = hitByID->begin(); hitByID_it!= hitByID->end(); ++hitByID_it ){
            cout<<"Cell/Straw id "<<hitByID_it->first<<" hit Data :"<<endl;
            hitByID_it->second->print();
    }


  } // end analyze

  void TestMapTrackerHitByID::endJob(){
  }


}  // end namespace mu2e

using mu2e::TestMapTrackerHitByID;
DEFINE_ART_MODULE(TestMapTrackerHitByID);
