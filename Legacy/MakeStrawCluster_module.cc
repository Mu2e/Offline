//
// $Id: MakeStrawCluster_module.cc,v 1.16 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
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
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawCluster.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
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
  void MakeStrawCluster::produce(art::Event  & evt)
  {
    if ( _diagLevel > 0 ) cout << "MakeStrawCluster: produce() begin" << endl;
    static int ncalls(0);
    ++ncalls;
    // A container to hold the output hits.
    unique_ptr<StrawClusterCollection> listofClusterspointer(new StrawClusterCollection);

    StrawCluster Cluster;
    StrawCluster tmpCluster;
    StrawClusterCollection&  listofClusters = *listofClusterspointer;

    StrawHitPtrVector ptrtoHits;

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
    StrawHitCollection const* hits = pdataHandle.product();
    for ( size_t i=0; i<hits->size(); ++i )
      {
        // Access data
        StrawHit hit = hits->at(i); 
        StrawIndex si = hit.strawIndex();
        Straw str = tracker.getStraw(si);         
        cout << " Hit Energy Deposition: "<< hit.energyDep()
             << " Hit time:    "<< hit.time()<<endl;
        //         StrawId sid = str.Id();
        bool used =false;
        for (size_t ii=0;ii<listofClusters.size();ii++)
          {         
            StrawCluster const& scluster =listofClusters.at(ii) ;
            StrawHitPtrVector const &  tmpptrtoHits = scluster.strawHits();
            for (size_t jj=0;jj<tmpptrtoHits.size();jj++)
              {
                StrawHit const& strawhit = *tmpptrtoHits[jj];
                if (strawhit==hit)
                  {                
                    used = true;
                    break;                                    }
              }
          }
        if ( !used )
          {
            ptrtoHits.push_back( StrawHitPtr(pdataHandle,i));
            const std::vector<StrawIndex> nearindex= str.nearestNeighboursByIndex();
            vector<StrawIndex>::const_iterator ncid;
            for(ncid=nearindex.begin(); ncid!=nearindex.end(); ncid++)
              {
                for ( size_t jj=0; jj<hits->size(); jj++ ) 
                  {
                    StrawHit nhit = hits->at(jj);
                    StrawIndex nsi = nhit.strawIndex();            
                    if (nsi==*ncid) ptrtoHits.push_back( StrawHitPtr(pdataHandle,jj));
                  } // end loop over all hits 
              } // end loop over neighbors 
            bool added=false; 
            if (ptrtoHits.size()>1) added = true;
            while (added)
              {        
                added = false;
                for(size_t kk=0;kk<ptrtoHits.size(); kk++)
                  {
                    StrawHit const& strawhit = *ptrtoHits[kk];
                    Straw straw = tracker.getStraw(strawhit.strawIndex());                      
                    const std::vector<StrawIndex> nnearindex= straw.nearestNeighboursByIndex();
                    vector<StrawIndex>::const_iterator nncid;
                    for(nncid=nnearindex.begin(); nncid!=nnearindex.end(); nncid++)
                      {
                        //
                        // first check if not already part of the cluster
                        bool usedincl=false;
                        for (size_t jjj=0;jjj<ptrtoHits.size();jjj++)
                          {
                            StrawHit const& strawhit = *ptrtoHits[jjj];
                            if (strawhit.strawIndex()==*nncid)
                              {
                                usedincl=true;
                                break;
                              }
                          }
                        if (!usedincl)
                          {
                            for ( size_t jj=0; jj<hits->size(); jj++ ) 
                              {
                                StrawHit const& nhit = hits->at(jj);
                                StrawIndex nsi = nhit.strawIndex();            
                                if (nsi==*nncid)
                                  {
                                    ptrtoHits.push_back( StrawHitPtr(pdataHandle,jj));
                                    added = true;
                                  }
                              } // end loop over all straws that fired
                          }  // end used in cluster     
                      }// end loop over neighbors
                  } // end loop over straws in cluster
              } // end while added
            if (ptrtoHits.size()>0) 
              {
                Cluster = StrawCluster(ptrtoHits);
                listofClusterspointer->push_back(Cluster);
                //listofptrtoHits->push_back(ptrtoHits);
              }
            ptrtoHits.clear();

          }// end if not used 
      } // end loop over all strawHits

    evt.put(std::move(listofClusterspointer));






    /*

    cout << "Number of StrawHits in Event:  "<< hits->size()<<endl;
    for ( size_t i=0; i<hits->size(); ++i ) {
    // Access data
    StrawHit        const&      hit(hits->at(i));
    //StrawIndex si = hit.strawIndex();
    //Straw str = tracker.getStraw(si);
    //StrawId sid = str.id();    
    ptrtoHits.push_back( StrawHitPtr(pdataHandle,i));

    }// loop over all straws that fired.
    cout << "size of index ptr vector:  " << ptrtoHits.size()<<endl;
    Cluster = StrawCluster(ptrtoHits);
    listofClusterspointer->push_back(Cluster);
    // Add the output hit collection to the event
    evt.put(std::move(listofClusterspointer));
    */
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
      //evt.put(std::move(listofClusterspointer));
    }
  } // end of ::produce.

}


using mu2e::MakeStrawCluster;
DEFINE_ART_MODULE(MakeStrawCluster);
