//
// $Id: MakeStrawCluster_module.cc,v 1.7 2011/06/01 21:11:24 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/01 21:11:24 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <cmath>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Provenance/Provenance.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawCluster.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "DataProducts/inc/DPIndexVector.hh"
#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class MakeStrawCluster : public art::EDProducer {
  public:
    explicit MakeStrawCluster(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _messageCategory("StrawClusterMaker"){

      // Tell the framework what we make.
      produces<StrawClusterCollection>();
    }
    virtual ~MakeStrawCluster() { }

    virtual void beginJob();
    void produce( art::Event& e);
    //  void analyze( art::Event const& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void MakeStrawCluster::beginJob(){

    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
         << endl;

    art::ServiceHandle<art::TFileService> tfs;
  }
   void MakeStrawCluster::produce(art::Event& evt)
   {
     if ( _diagLevel > 0 ) cout << "MakeStrawCluster: produce() begin" << endl;
     static int ncalls(0);
     ++ncalls;
     // A container to hold the output hits.
     auto_ptr<StrawClusterCollection> listofClusterspointer(new StrawClusterCollection);

     StrawCluster Cluster;
     StrawCluster tmpCluster;
     StrawClusterCollection&  listofClusters = *listofClusterspointer;

     //     DPIndexVector ptrtoHits;

     std::vector<DPIndex> ptrtoHits;
     //     StrawCluster::const_iterator ostrawIter;
     //    StrawClusterCollection::const_iterator oClusterIter;
     //StrawCluster::const_iterator istrawIter;
     //StrawClusterCollection::const_iterator iClusterIter;

     // Geometry info for the TTracker.
     // Get a reference to one of the L or T trackers.
     // Throw exception if not successful.
     const Tracker& tracker = getTrackerOrThrow();
     // source of the info:
     art::Handle<StrawHitCollection> pdataHandle;
     evt.getByLabel(_makerModuleLabel,pdataHandle);
     bool haveStrawHits = pdataHandle.isValid();
     if ( _diagLevel > -1 && !haveStrawHits) cout << __func__ << ": No StrawHits" << endl;
    
     if( !haveStrawHits) return;
     art::ProductID const& id(pdataHandle.id());
     StrawHitCollection const* hits = pdataHandle.product();

     cout << "Number of StrawHits in Event:  "<< hits->size()<<endl;
     for ( size_t i=0; i<hits->size(); ++i ) {
       // Access data
       StrawHit        const&      hit(hits->at(i));
       //StrawIndex si = hit.strawIndex();
       //Straw str = tracker.getStraw(si);
       //StrawId sid = str.id();    
       ptrtoHits.push_back(DPIndex(id,i));

     }// loop over all straws that fired.
     cout << "size of index ptr vector:  " << ptrtoHits.size()<<endl;
     Cluster = StrawCluster(ptrtoHits);
     listofClusterspointer->push_back(Cluster);
    // Add the output hit collection to the event
    evt.put(listofClusterspointer);

    if (_diagLevel>2){
      cout << " Nr of Hits:  "<< hits->size()<<endl;
      //cout << " nr of clusters:  " <<listofClusters.size()<<endl;
      //for(oClusterIter=listofClusters.begin();oClusterIter!=listofClusters.end(); oClusterIter++)
      //{
      //    tmpCluster= *oClusterIter;
      //    cout <<" Cluster length: "<< tmpCluster.size()<<endl;
      //    for(ostrawIter=tmpCluster.begin();ostrawIter!=tmpCluster.end(); ostrawIter++)
      //    {
      //       cout<<*ostrawIter<<endl;
      //      }
      //  }
      // }
      //evt.put(listofClusterspointer);
    }
  } // end of ::produce.

}


using mu2e::MakeStrawCluster;
DEFINE_ART_MODULE(MakeStrawCluster);
