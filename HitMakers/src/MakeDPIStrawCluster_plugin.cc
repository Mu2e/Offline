//
// $Id: MakeDPIStrawCluster_plugin.cc,v 1.1 2011/01/28 21:38:12 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/28 21:38:12 $
//
// Original author Hans Wenzel
// This modules create clusters of fired StrawHits in a panel 
// the clusters are stored in a DPIndexVectorCollection pointing 
// back to the strawhits used. 

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Provenance/interface/Provenance.h"

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
#include "ToyDP/inc/DPIndexVector.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
using namespace std;

namespace mu2e {
  
  //--------------------------------------------------------------------
  //
  // 
  class MakeDPIStrawCluster : public edm::EDProducer {
  public:
    explicit MakeDPIStrawCluster(edm::ParameterSet const& pset):
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel")),
      _messageCategory("StrawClusterMaker"){

      // Tell the framework what we make.
      produces<DPIndexVectorCollection>("DPIStrawCluster");
    }
    virtual ~MakeDPIStrawCluster() { }
    
    virtual void beginJob(edm::EventSetup const&);
    void produce( edm::Event& e, edm::EventSetup const&);
    
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
     
  void MakeDPIStrawCluster::beginJob(edm::EventSetup const& ){
    
    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;
    
    edm::Service<edm::TFileService> tfs;
  }
   void MakeDPIStrawCluster::produce(edm::Event& evt, edm::EventSetup const&)
   {
     if ( _diagLevel > 0 ) cout << "MakeDPIStrawCluster: produce() begin" << endl;
     static int ncalls(0);
     ++ncalls;
     // A container to hold the output.
     auto_ptr<DPIndexVectorCollection>   listofptrtoHits(new DPIndexVectorCollection);
     DPIndexVector ptrtoHits;
     
     // Geometry info for the TTracker.
     // Get a reference to one of the L or T trackers.
     // Throw exception if not successful.
     const Tracker& tracker = getTrackerOrThrow();
    // Ask the event to give us a handle to the requested StrawHits
     edm::Handle<StrawHitCollection> pdataHandle;
     evt.getByLabel(_makerModuleLabel,pdataHandle);
     edm::ProductID const& id(pdataHandle.id());
     StrawHitCollection const* hits = pdataHandle.product();
     //    sort(hits->begin(),hits->end());
     for ( size_t i=0; i<hits->size(); ++i )
       {
	 // Access data
	 StrawHit hit = hits->at(i); 
	 StrawIndex si = hit.strawIndex();
	 Straw str = tracker.getStraw(si);	 
	 //	 StrawId sid = str.Id();
	 bool used =false;
	 for (size_t ii=0;ii<listofptrtoHits->size();ii++)
	   {  
	     DPIndexVector tmpptrtoHits=listofptrtoHits->at(ii);
	     for (size_t jj=0;jj<tmpptrtoHits.size();jj++)
	       {
		 DPIndex const& junkie = tmpptrtoHits[jj];
		 StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
		 if (strawhit==hit)
		 {		
		   used = true;
		   break;		   		 }
	       }
	   }
	 if ( !used )
	   {
	     ptrtoHits.push_back(DPIndex(id,i));
	     const std::vector<StrawIndex> nearindex= str.nearestNeighboursByIndex();
	     vector<StrawIndex>::const_iterator ncid;
	     for(ncid=nearindex.begin(); ncid!=nearindex.end(); ncid++)
	       {
		 for ( size_t jj=0; jj<hits->size(); jj++ ) 
		   {
		     StrawHit nhit = hits->at(jj);
		     StrawIndex nsi = nhit.strawIndex();	    
		     if (nsi==*ncid) ptrtoHits.push_back(DPIndex(id,jj));
		   } // end loop over all hits 
	       } // end loop over neighbors 
	     bool added=false; 
	     if (ptrtoHits.size()>1) added = true;
	     while (added)
	       {	
		 added = false;
		 for(size_t kk=0;kk<ptrtoHits.size(); kk++)
		   {
		     DPIndex const& junkie = ptrtoHits[kk];
		     StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
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
			     DPIndex const& junkie = ptrtoHits[jjj];
			     StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
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
				     ptrtoHits.push_back(DPIndex(id,jj));
				     added = true;
				   }
			       } // end loop over all straws that fired
			   }  // end used in cluster     
		       }// end loop over neighbors
		   } // end loop over straws in cluster
	       } // end while added
	     if (ptrtoHits.size()>0) 
	       {
		 listofptrtoHits->push_back(ptrtoHits);
	       }
	     ptrtoHits.clear();

	   }// end if not used 
       } // end loop over all strawHits
     if ( _diagLevel > 1 )
       {      
	 cout << "Number of Clusters: "<< listofptrtoHits->size()<<"  Number of Hits:  "<<  hits->size()<<endl;
	 for (size_t ii=0;ii<listofptrtoHits->size();ii++)
	   {
	     DPIndexVector tmpptrtoHits=listofptrtoHits->at(ii);
	     cout << "Cluster Nr. "<<ii<<":  ";
	     for(size_t kk=0;kk<tmpptrtoHits.size(); kk++)
	       { 
	       DPIndex const& junkie = tmpptrtoHits[kk];
	       StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
	       cout <<"  "<<strawhit.strawIndex();
	       }
	     cout<<endl;
	   }
       }
     // Add the output  to the event
     evt.put(listofptrtoHits,"DPIStrawCluster");		
   } // end of ::produce.
  
}

using mu2e::MakeDPIStrawCluster;
DEFINE_FWK_MODULE(MakeDPIStrawCluster);
