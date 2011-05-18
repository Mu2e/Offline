//
// $Id: MakeDPIStrawCluster_module.cc,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
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
#include "ToyDP/inc/DPIndexVector.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class MakeDPIStrawCluster : public art::EDProducer {
  public:
    explicit MakeDPIStrawCluster(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _messageCategory("StrawClusterMaker"){

      // Tell the framework what we make.
      produces<DPIndexVectorCollection>("DPIStrawCluster");
    }
    virtual ~MakeDPIStrawCluster() { }

    virtual void beginJob();
    void produce( art::Event& e);

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

  void MakeDPIStrawCluster::beginJob(){

    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
         << endl;

    art::ServiceHandle<art::TFileService> tfs;
  }
   void MakeDPIStrawCluster::produce(art::Event& evt)
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
     art::Handle<StrawHitCollection> pdataHandle;
     evt.getByLabel(_makerModuleLabel,pdataHandle);
     art::ProductID const& id(pdataHandle.id());
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
DEFINE_ART_MODULE(MakeDPIStrawCluster);
