//
// $Id: MakeStrawCluster_module.cc,v 1.2 2011/05/17 22:22:46 wb Exp $
// $Author: wb $
// $Date: 2011/05/17 22:22:46 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <cmath>


// Framework includes.
//#include "art/Framework/Core/EDAnalyzer.h"
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

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawCluster.hh"
#include "ToyDP/inc/StrawClusterCollection.hh"
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
     auto_ptr<StrawClusterCollection>        listofClusterspointer(new StrawClusterCollection);
     
     StrawCluster Cluster;
     StrawCluster tmpCluster;
     StrawClusterCollection&  listofClusters = *listofClusterspointer;

     StrawCluster::const_iterator ostrawIter;
     StrawClusterCollection::const_iterator oClusterIter;    
     StrawCluster::const_iterator istrawIter;
     StrawClusterCollection::const_iterator iClusterIter;
     
     // Geometry info for the TTracker.
     // Get a reference to one of the L or T trackers.
     // Throw exception if not successful.
     const Tracker& tracker = getTrackerOrThrow();
     art::Handle<StrawHitCollection> pdataHandle;
     evt.getByLabel(_makerModuleLabel,pdataHandle);
     StrawHitCollection const* hits = pdataHandle.product();
     for ( size_t i=0; i<hits->size(); ++i ) {
       // Access data
       StrawHit        const&      hit(hits->at(i));
       StrawIndex si = hit.strawIndex();
       Straw str = tracker.getStraw(si);	 
       StrawId sid = str.Id();
       // first check if  straw already has been used
       bool used =false;
       for(oClusterIter=listofClusters.begin();oClusterIter!=listofClusters.end(); oClusterIter++)
	 {
	   tmpCluster= *oClusterIter;
	   for(ostrawIter=tmpCluster.begin();ostrawIter!=tmpCluster.end(); ostrawIter++)
	     {
	       if (sid == *ostrawIter)
		 {
		   used = true;
		   break;
		 }
	     }
	 }
       if ( !used )
	 {
	   Cluster.push_back(sid);
	   // get list of neighbors and check if they fired:
	   const std::vector<StrawId> nearid= str.nearestNeighboursById();
	   vector<StrawId>::const_iterator ncid;
	   for(ncid=nearid.begin(); ncid!=nearid.end(); ncid++)
	     {
	       for ( size_t jj=0; jj<hits->size(); ++jj ) 
		 {
		   StrawHit        const&      hit(hits->at(jj));
		   StrawIndex nsi = hit.strawIndex();	    
		   Straw nstr = tracker.getStraw(nsi);	 
		   StrawId nsid = nstr.Id();
		   if (nsid==*ncid)
		     {
		       bool nused =false;
		       for(oClusterIter=listofClusters.begin();oClusterIter!=listofClusters.end(); oClusterIter++)
			 {
			   tmpCluster= *oClusterIter;
			   for(ostrawIter=tmpCluster.begin();ostrawIter!=tmpCluster.end(); ostrawIter++)
			     {
			       if (nsid == *ostrawIter)
				 {
				   nused = true;
				   break;
				}
			     }
			 }
		       if ( !nused) Cluster.push_back(nsid);
		     } 
		 } // end loop over all hits 
	     } // end loop over neighbors 
	   bool added=false; 
	   if (Cluster.size()>1) added = true;
	   while (added)
	     {
	       added = false;
	       for(size_t kk=0;kk<Cluster.size(); kk++)
		 {
		   Straw straw = tracker.getStraw(Cluster[kk]);		      
		   const std::vector<StrawId> nnearid= straw.nearestNeighboursById();
		   vector<StrawId>::const_iterator nncid;
		   for(nncid=nnearid.begin(); nncid!=nnearid.end(); nncid++)
		     {
		       //
		       // first check if not already part of the cluster
		       //
		       vector<StrawId>::const_iterator sIter;
		       bool usedincl=false;
		      for(sIter=Cluster.begin();sIter!=Cluster.end(); sIter++)
			{
			  if (*sIter==*nncid) 
			    {
			      usedincl=true;
			      break;
			    }
			}
		      if (!usedincl)
			{
			  for ( size_t jj=0; jj<hits->size(); jj++ ) 
			    {
			      StrawHit        const&      hit(hits->at(jj));
			      StrawIndex nsi = hit.strawIndex();	    
			      Straw nstr = tracker.getStraw(nsi);	 
			      StrawId nsid = nstr.Id();
			      if (nsid==*nncid)
				{
				  bool nused =false;
				  for(iClusterIter=listofClusters.begin();iClusterIter!=listofClusters.end(); iClusterIter++)
				    {
				      tmpCluster= *iClusterIter;
				      for(istrawIter=tmpCluster.begin();istrawIter!=tmpCluster.end(); istrawIter++)
					{
					  if (nsid == *istrawIter)
					    {
					      nused = true;
					      break;
					    }
					}
				    }
				  if ( !nused)
				    {
				      Cluster.push_back(nsid);
				      added = true;
				    }
				}
			    } // end loop over all straws that fired
			}  // end used in cluster     
		    }// end loop over neighbors
		} // end loop over straws in cluster
	    } // end while added
	  if (Cluster.size()>0) 
	    {
	      listofClusters.push_back(Cluster);
	    }
	  Cluster.clear();
	}// loop over all straws that fired.
    }
    // Add the output hit collection to the event
    //evt.put(listofClusterspointer);
     
    if (_diagLevel>2){
      cout << " Nr of Hits:  "<< hits->size()<<endl;
      cout << " nr of clusters:  " <<listofClusters.size()<<endl;
      for(oClusterIter=listofClusters.begin();oClusterIter!=listofClusters.end(); oClusterIter++)
	{
	  tmpCluster= *oClusterIter;
	  cout <<" Cluster length: "<< tmpCluster.size()<<endl;
	  for(ostrawIter=tmpCluster.begin();ostrawIter!=tmpCluster.end(); ostrawIter++)
	    {
	      cout<<*ostrawIter<<endl;
	    }
	}
    }
    evt.put(listofClusterspointer);			
  } // end of ::analyze.
  
}


using mu2e::MakeStrawCluster;
DEFINE_ART_MODULE(MakeStrawCluster);
