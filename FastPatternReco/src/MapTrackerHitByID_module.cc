//
// map Tracker Hits to be accessed by Cell/Straw ID
//
// $Id: MapTrackerHitByID_module.cc,v 1.2 2013/03/14 19:47:45 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/14 19:47:45 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <utility>
#include <limits>

#include <boost/shared_ptr.hpp>


// Framework includes.
//#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
//#include "GeometryService/inc/GeometryService.hh"
//#include "GeometryService/inc/GeomHandle.hh"
//#include "GeometryService/inc/getTrackerOrThrow.hh"
//#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
//#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"

using namespace std;

namespace mu2e {

class MapTrackerHitByID : public art::EDProducer {
public:

        explicit MapTrackerHitByID(fhicl::ParameterSet const& pset);
        virtual ~MapTrackerHitByID() {
                //if (_fakeCanvas)        delete _fakeCanvas;
        }

        virtual void beginJob();
        void endJob();

        // This is called for each event.
        //void analyze(art::Event const& e);
        void produce(art::Event & e);

private:

        // Start: run time parameters

        // Label of the module that made the hits.
        std::string _makerModuleLabel;

};

MapTrackerHitByID::MapTrackerHitByID(fhicl::ParameterSet const& pset) :

      // Run time parameters
      _makerModuleLabel(pset.get<string>("makerModuleLabel"))
{
        // Tell the framework what we make.
        produces<TrackerHitByID>();

}

void MapTrackerHitByID::beginJob(){

        cout<<"Starting mapping Tracker hit to by accessed by Cell/Straw ID!"<<endl;

}

//  void MapTrackerHitByID::analyze(art::Event const& event ) {
void MapTrackerHitByID::produce(art::Event & event ) {

        auto_ptr<TrackerHitByID> trkHits(new TrackerHitByID);

//        const Tracker& tracker = getTrackerOrThrow();
//        const TTracker &ttr = static_cast<const TTracker&>( tracker );
//        const std::vector<Device> ttrdev = ttr.getDevices();

        art::Handle<StrawHitCollection> pdataHandle;
        event.getByLabel(_makerModuleLabel,pdataHandle);
        StrawHitCollection const* hits = pdataHandle.product();

        size_t nStrawPerEvent = hits->size();
        unsigned long int tmpId;

        for (size_t i=0; i<nStrawPerEvent; ++i) {

                // Access data
                StrawHit        const&      hit(hits->at(i));

                //Get hit straw
                //StrawIndex si = hit.strawIndex();
                tmpId = hit.strawIndex().asUint();
                trkHits->insert( TrackerHitByID::value_type(tmpId,art::Ptr<mu2e::StrawHit>( pdataHandle, i) ) );

        }
        event.put(std::move(trkHits));


  } // end produce

  void MapTrackerHitByID::endJob(){
  }


}  // end namespace mu2e

using mu2e::MapTrackerHitByID;
DEFINE_ART_MODULE(MapTrackerHitByID);
