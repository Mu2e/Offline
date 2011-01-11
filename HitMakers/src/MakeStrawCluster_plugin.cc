//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: MakeStrawCluster_plugin.cc,v 1.5 2011/01/11 17:16:27 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/11 17:16:27 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <cmath>


// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
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

using namespace std;

namespace mu2e {
  
  
  class Vector
  {
  public:
    float x_, y_;  
    Vector(float f = 0.0f)
      : x_(f), y_(f) {}
    
    Vector(float x, float y)
      : x_(x), y_(y) {}
  };
  
  class LineSegment
  {
  public:
    Vector begin_;
    Vector end_;
    
    LineSegment(const Vector& begin, const Vector& end)
      : begin_(begin), end_(end) {}
    
    enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERSECTING, INTERSECTING };
    
    IntersectResult Intersect(const LineSegment& other_line, Vector& intersection)
    {
      float denom = ((other_line.end_.y_ - other_line.begin_.y_)*(end_.x_ - begin_.x_)) -
	((other_line.end_.x_ - other_line.begin_.x_)*(end_.y_ - begin_.y_));
      
      float nume_a = ((other_line.end_.x_ - other_line.begin_.x_)*(begin_.y_ - other_line.begin_.y_)) -
	((other_line.end_.y_ - other_line.begin_.y_)*(begin_.x_ - other_line.begin_.x_));
      
      float nume_b = ((end_.x_ - begin_.x_)*(begin_.y_ - other_line.begin_.y_)) -
	((end_.y_ - begin_.y_)*(begin_.x_ - other_line.begin_.x_));
      
      if(denom == 0.0f)
	{
	  if(nume_a == 0.0f && nume_b == 0.0f)
	    {
	      return COINCIDENT;
	    }
	  return PARALLEL;
	}
      
      float ua = nume_a / denom;
      float ub = nume_b / denom;
      
      if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
	{
	  // Get the intersection point.
	  intersection.x_ = begin_.x_ + ua*(end_.x_ - begin_.x_);
	  intersection.y_ = begin_.y_ + ua*(end_.y_ - begin_.y_);
	  
	  return INTERSECTING;
	}
      
      return NOT_INTERSECTING;
    }
  };
  //--------------------------------------------------------------------
  //
  // 
  class MakeStrawCluster : public edm::EDAnalyzer {
  public:
    explicit MakeStrawCluster(edm::ParameterSet const& pset):
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel"))
    {
    }
    virtual ~MakeStrawCluster() { }
    
    virtual void beginJob(edm::EventSetup const&);
    
    void analyze( edm::Event const& e, edm::EventSetup const&);
    
  private:
    
    // Diagnostics level.
    int _diagLevel;
    
    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    
    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;
    
    // Label of the module that made the hits.
    std::string _makerModuleLabel;
    
  };
  
  void MakeStrawCluster::beginJob(edm::EventSetup const& ){
    
    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;
    
    edm::Service<edm::TFileService> tfs;
  }
  
  void
  MakeStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&) {
     
    static int ncalls(0);
    ++ncalls;
    vector<StrawId> Cluster;
    vector<StrawId> tmpCluster;
    vector<vector<StrawId> > listofClusters;
    vector<StrawId>::const_iterator ostrawIter;
    vector<vector<StrawId> >::const_iterator oClusterIter;
    
    vector<StrawId>::const_iterator istrawIter;
    vector<vector<StrawId> >::const_iterator iClusterIter;
    
    // Geometry info for the TTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();
    edm::Handle<StrawHitCollection> pdataHandle;
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
  } // end of ::analyze.
  
}


using mu2e::MakeStrawCluster;
DEFINE_FWK_MODULE(MakeStrawCluster);
