//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: MakeStrawCluster_plugin.cc,v 1.3 2011/01/07 18:08:23 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/07 18:08:23 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <cmath>
#include <algorithm>
#include <utility>
#include <map>


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
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveTransients.hh"

using namespace std;

namespace mu2e {
class straw{

 public:
    Int_t   evt;
    Int_t   id;
    Int_t   lay;
    Int_t   did;
    Int_t   sec;
    Float_t hl;
    Float_t mpx;
    Float_t mpy;
    Float_t mpz;
    Float_t dirx;
    Float_t diry;
    Float_t dirz;
    Int_t PanelIndex()
    { 
    const int spdev   = 600;
    const int sppanel = 100;
    //const int splay   = 50;
      return id-(did*spdev)-(sec*sppanel);
    }
    bool operator>(const straw other) const {
      if (id > other.id) {
	return true;
      }
      else{
	return false;
      }
    }
   bool operator<(const straw other) const {
      if (id < other.id) {
	return true;
      }
      else{
	return false;
      }
   }
   bool operator==(const straw other) const {
      if (id == other.id) {
	return true;
      }
      else{
	return false;
      }
   }
    void Print()
    {
      cout<< "Straw:  " <<endl;
      cout<< "======  " <<endl;
      cout<< "Event:  "<<evt<<endl;
      cout<< "SID:  "<<id<<endl;
      cout<< "Layer:  "<<lay<<endl;
      cout<< "DID:    "<<did<<endl;
      cout<< "Sector: "<< sec<<endl;
      //    cout<< ":"<<hl<<endl;
      //cout<< "evt:"<<mpx<<endl;
      //cout<< "evt:"<< mpy<<endl;
      //cout<< "evt:"<<mpz<<endl;
      //cout<< "evt:"<< dirx<<endl;
      //cout<< "evt:"<<diry<<endl;
      //cout<< "evt:"<<dirz<<endl;     
    }
 };


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
    //    vector<StrawId>  ListofNeighbors(StrawId si);
    // vector<StrawId>  createStrawCluster(StrawId si);


    //    vector<vector<StrawId>> StrawCluster(

  };

  void MakeStrawCluster::beginJob(edm::EventSetup const& ){

    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;

    edm::Service<edm::TFileService> tfs;

;
  }

  void
  MakeStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;
    vector<StrawId> Cluster;
    vector<StrawId> tmpCluster;
    vector<StrawId>::const_iterator strawIter;
    vector<vector<StrawId> > listofClusters;
    vector<vector<StrawId> >::const_iterator ClusterIter;

    // Geometry info for the LTracker.
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
      //int sindex = si.asInt();
      Straw str = tracker.getStraw(si);	 
      StrawId sid = str.Id();
      cout << "Straw: " << sid<< endl;
      //      vector<StrawId> cluster = createStrawCluster(sid);
      const std::vector<StrawId> nearid= str.nearestNeighboursById();
      cout << "nr of neighbors:  " << nearid.size()<<endl;
      vector<StrawId>::const_iterator ncid;
      for(ncid=nearid.begin(); ncid!=nearid.end(); ncid++)
	{
	  cout << "cc StrawId: " <<*ncid << endl;
	  bool used =false;
	  for(ClusterIter=listofClusters.begin();ClusterIter!=listofClusters.end(); ClusterIter++)
	    {
	      tmpCluster= *ClusterIter;
	      for(strawIter=tmpCluster.begin();strawIter!=tmpCluster.end(); strawIter++)
		{
		  if (sid == *strawIter)
		    {
		      used = true;
		      break;
		    }
		}
	    }
	  if ( used==false ) Cluster.push_back(sid);
	  for ( size_t jj=0; jj<hits->size(); ++jj ) 
	    {
	      StrawHit        const&      hit(hits->at(jj));
	      StrawIndex nsi = hit.strawIndex();	    
	      Straw nstr = tracker.getStraw(nsi);	 
	      StrawId nsid = nstr.Id();
	      if (nsid==*ncid)
		{
		  cout<< " fired"<<endl;
		}
	    }
	  
	}
    }


  } // end of ::analyze.
/*
vector<StrawId>   MakeStrawCluster::createStrawCluster(StrawId sid)
 { 
   //   vector<StrawId> neighbors = ListofNeighbors(sid);
   // cout<< " number of neighbors: " << neighbors.size()<<endl;
   //vector<StrawId>::const_iterator ncii;
   //   for(ncii=neighbors.begin(); ncii!=neighbors.end(); ncii++)
   //  {
   //    cout << "Straw index: " <<*ncii << endl;
   //  }

   const Tracker& tracker = getTrackerOrThrow();
   edm::Handle<StrawHitCollection> pdataHandle;
   evt.getByLabel(_makerModuleLabel,pdataHandle);
   StrawHitCollection const* hits = pdataHandle.product();
   Straw str = tracker.getStraw(sid);	 
   //StrawId sid = str.Id();
   const std::vector<StrawId> nearid= str.nearestNeighboursById();
   cout << "nr of neighbors:  " << nearid.size()<<endl;
   vector<StrawId>::const_iterator ncid;
   for(ncid=nearid.begin(); ncid!=nearid.end(); ncid++)
     {
       cout << "cc StrawId: " <<*ncid << endl;
       for ( size_t jj=0; jj<hits->size(); ++jj ) {
	 StrawHit        const&      hit(hits->at(jj));
	 StrawIndex nsi = hit.strawIndex();	    
	 Straw nstr = tracker.getStraw(nsi);	 
	 StrawId nsid = nstr.Id();
	 if (nsid==*ncid)
	   {
	     cout<< " fired"<<endl;
	   }
       }
       
     }
   




    return nearid;    

}//end createStrawCluster
*/
//
// use for now in the end we probably want a straw to be able what it's neighbors are. 
//
/*
vector<StrawId>   MakeStrawCluster::ListofNeighbors(StrawId sid)
 {   
   LayerId lid      = sid.getLayerId();
   SectorId   secid = sid.getSectorId();
   int lay          = sid.getLayer();
   int sec          = sid.getSector(); 
   int did          = sid.getDevice(); 
   int ind          = sid.getStraw();
   StrawId si;
   const int splay  = 50;    // number of straws per panel layer
   vector<StrawId> neighbors;
   
   cout << "layer:  "<<lay<<"  Sector:  "<<sec<<"  Device:  "<< did<<" Straw:  "<<ind<<endl;
   if (lay==0)
     {
       if (ind-1>0)
	 {
	   si = StrawId(lid,ind-1);
	   neighbors.push_back(si);
	   si = StrawId(secid,lay+1,ind-1);
	   neighbors.push_back(si);
	 }
       if (ind+1<splay)
	 {
	   si = StrawId(lid,ind+1);
	   neighbors.push_back(si);
	 } 
	   si = StrawId(secid,lay+1,ind);
	   neighbors.push_back(si);
     }
   else
     {
       if (ind-1>0)
	 {
	   si = StrawId(lid,ind-1);
	   neighbors.push_back(si);
	 }
       if (ind+1<splay)
	 {
	   si = StrawId(lid,ind+1);
	   neighbors.push_back(si);
	   si = StrawId(secid,lay-1,ind+1);
	   neighbors.push_back(si);
	 }       
       si = StrawId(secid,lay-1,ind);
       neighbors.push_back(si);	
     }
    return neighbors;    

}//end of ListofNeighbors
*/ 
}


using mu2e::MakeStrawCluster;
DEFINE_FWK_MODULE(MakeStrawCluster);
