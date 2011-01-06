//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: MakeStrawCluster_plugin.cc,v 1.1 2011/01/06 00:11:48 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/06 00:11:48 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
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
  //const int spdev   = 600;
  //const int sppanel = 100;
  //const int splay   = 50;
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
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel")),
      _hHitTime(0),
      _hHitDeltaTime(0),
      _hHitEnergy(0),
      _hNHits(0),
      _hNHitsPerWire(0),
      _hDriftTime(0),
      _hDriftDistance(0),
      _hDistanceToMid(0),
      _hNG4Steps(0),
      _hT0(0),
      _hG4StepLength(0),
      _hG4StepEdep(0),
      _ntup(0)
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

    // Some diagnostic histograms.
    TH1F* _hHitTime;
    TH1F* _hHitDeltaTime;
    TH1F* _hHitEnergy;
    TH1F* _hNHits;
    TH1F* _hNHitsPerWire;
    TH1F* _hDriftTime;
    TH1F* _hDriftDistance;
    TH1F* _hDistanceToMid;
    TH1F* _hNG4Steps;
    TH1F* _hT0;
    TH1F* _hG4StepLength;
    TH1F* _hG4StepEdep;
    TNtuple* _ntup;
    vector<int>  ListofNeighbors(int straw);
    vector<StrawId>  nListofNeighbors(StrawId si);
    //void nListofNeighbors(StrawId si);
  };

  void MakeStrawCluster::beginJob(edm::EventSetup const& ){

    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;

    edm::Service<edm::TFileService> tfs;

    _hHitTime      = tfs->make<TH1F>( "hHitTime",      "Hit Time (ns)", 200, 0., 2000. );
    _hHitDeltaTime = tfs->make<TH1F>( "hHitDeltaTime", "Hit Delta Time (ns)", 80, -20.0, 20. );
    _hHitEnergy    = tfs->make<TH1F>( "hHitEnergy",    "Hit Energy (keV)", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",        "Number of straw hits", 500, 0., 500. );
    _hNHitsPerWire = tfs->make<TH1F>( "hNHitsPerWire", "Number of hits per straw", 10, 0., 10. );
    _hDriftTime    = tfs->make<TH1F>( "hDriftTime",    "Drift time, ns", 100, 0., 100. );
    _hDriftDistance= tfs->make<TH1F>( "hDriftDistance","Drift Distance, mm", 100, 0., 3. );
    _hDistanceToMid= tfs->make<TH1F>( "hDistanceToMid","Distance to wire center, mm", 160, -1600., 1600. );
    _hNG4Steps     = tfs->make<TH1F>( "hNG4Steps",     "Number of G4Steps per hit", 100, 0., 100. );
    _hT0           = tfs->make<TH1F>( "hT0",           "T0, ns", 100, -50., 50. );
    _hG4StepLength = tfs->make<TH1F>( "hG4StepLength", "Length of G4Steps, mm", 100, 0., 10. );
    _hG4StepEdep   = tfs->make<TH1F>( "hG4StepEdep",   "Energy deposition of G4Steps, keV", 100, 0., 10. );
    _ntup          = tfs->make<TNtuple>( "ntup", "Straw Hit ntuple", 
                      "evt:lay:did:sec:hl:mpx:mpy:mpz:dirx:diry:dirz:time:dtime:eDep:driftT:driftDistance:distanceToMid");
  }

  void
  MakeStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;


    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

  
    edm::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.

    edm::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs

    edm::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than 
    // use producer name directly ("g4run"). 

    edm::Handle<StepPointMCCollection> mchitsHandle;
    evt.getByLabel("g4run",_trackerStepPoints,mchitsHandle);
    StepPointMCCollection const* mchits = mchitsHandle.product();

    // Fill histograms

    _hNHits->Fill(hits->size());

    std::map<StrawIndex,int> nhperwire;

    for ( size_t i=0; i<hits->size(); ++i ) {

      // Access data
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));
      
      // Fill per-event histograms
      if( i==0 ) {
        _hT0->Fill(truth.t0());
      }

      // Use data from hits
      _hHitTime->Fill(hit.time());
      _hHitDeltaTime->Fill(hit.dt());
      _hHitEnergy->Fill(hit.energyDep()*1000.0);

      // Use MC truth data
      _hDriftTime->Fill(truth.driftTime());
      _hDriftDistance->Fill(truth.driftDistance());
      _hDistanceToMid->Fill(truth.distanceToMid());

      // Use data from G4 hits
      _hNG4Steps->Fill(mcptr.size());
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
        _hG4StepLength->Fill(mchit.stepLength());
        _hG4StepEdep->Fill(mchit.eDep()*1000.0);
      }
      StrawIndex si = hit.strawIndex();
      int sindex = si.asInt();
      Straw str = tracker.getStraw(si);	
      StrawId sid = str.Id();
      vector<StrawId> nneighbors = nListofNeighbors(sid);
      cout<< " number of neighbors: " << nneighbors.size()<<endl;
      vector<StrawId>::const_iterator ncii;
      for(ncii=nneighbors.begin(); ncii!=nneighbors.end(); ncii++)
	{
	  cout << "Straw index: " <<*ncii << endl;
	}
      LayerId lid = sid.getLayerId();
      vector<int> neighbors =  ListofNeighbors(sid.getStraw());
      vector<int>::const_iterator cii;
      for(cii=neighbors.begin(); cii!=neighbors.end(); cii++)
	{
	  cout << "Straw index: " <<sindex << "   Neighbor:  "<<*cii << endl;
	}
      DeviceId did = sid.getDeviceId();
      SectorId secid = sid.getSectorId();
      float nt[17];
      const CLHEP::Hep3Vector vec3junk = str.getMidPoint();
      const CLHEP::Hep3Vector vec3junk1 = str.getDirection();
      // Fill the ntuple:
      nt[0]  = evt.id().event();
      nt[1]  = lid.getLayer();
      nt[2]  = did;
      nt[3]  = secid.getSector();
      nt[4]  = str.getHalfLength();
      nt[5]  = vec3junk.getX();
      nt[6]  = vec3junk.getY();
      nt[7]  = vec3junk.getZ();
      nt[8]  = vec3junk1.getX();
      nt[9]  = vec3junk1.getY();
      nt[10] = vec3junk1.getZ();
      nt[11] = hit.time();
      nt[12] = hit.dt();
      nt[13] = hit.energyDep();
      nt[14] = truth.driftTime();
      nt[15] = truth.driftDistance();
      nt[16] = truth.distanceToMid();
      _ntup->Fill(nt);
      // Calculate number of hits per wire
      ++nhperwire[hit.strawIndex()];

    }

    for( std::map<StrawIndex,int>::iterator it=nhperwire.begin(); it!= nhperwire.end(); ++it ) {
      _hNHitsPerWire->Fill(it->second);
    }

  } // end of ::analyze.

 vector<int>   MakeStrawCluster::ListofNeighbors(int straw)
 {
   // some info about the ttracker:
   const int spdev   = 600;   // number of straws per device
   const int sppanel = 100;   // number of straws per panel 
   const int splay   = 50;    // number of straws per panel layer
   vector<int> neighbors;
   int did = straw/spdev;
   int sec = (straw-(did*spdev))/sppanel;
   int lay = (straw-(did*spdev)-(sec*sppanel))/splay;
   int ind =  straw-(did*spdev)-(sec*sppanel);

   if (lay==0) 
     {
       if ((ind-1)>-1)              neighbors.push_back(did*spdev+sec*sppanel+ind-1);
       if ((ind+1)<splay)           neighbors.push_back(did*spdev+sec*sppanel+ind+1);
       if ((ind+splay-1)>(splay-1)) neighbors.push_back(did*spdev+sec*sppanel+ind+splay-1);
       if ((ind+splay)<sppanel)     neighbors.push_back(did*spdev+sec*sppanel+ind+splay);
     }
   else 
     {
       if ((ind-splay)>-1)          neighbors.push_back(did*spdev+sec*sppanel+ind-splay);
       if ((ind-splay+1)<splay)     neighbors.push_back(did*spdev+sec*sppanel+ind-splay+1);
       if ((ind-1)>(splay-1))       neighbors.push_back(did*spdev+sec*sppanel+ind-1);
       if ((ind+1)<sppanel)         neighbors.push_back(did*spdev+sec*sppanel+ind+1);
     }

   return neighbors; 
}//end of ListofNeighbors 
vector<StrawId>   MakeStrawCluster::nListofNeighbors(StrawId sid)
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

   // some info about the ttracker:
     //   const int spdev   = 600;   // number of straws per device
     //const int sppanel = 100;   // number of straws per panel 
     //const int splay   = 50;    // number of straws per panel layer
   //   vector<StrawId> neighbors;
   //int did = straw/spdev;
   //int sec = (straw-(did*spdev))/sppanel;
   //int lay = (straw-(did*spdev)-(sec*sppanel))/splay;
   //int ind =  straw-(did*spdev)-(sec*sppanel);
   /* 
   if (lay==0) 
     {
       if ((ind-1)>-1)              neighbors.push_back(did*spdev+sec*sppanel+ind-1);
       if ((ind+1)<splay)           neighbors.push_back(did*spdev+sec*sppanel+ind+1);
       if ((ind+splay-1)>(splay-1)) neighbors.push_back(did*spdev+sec*sppanel+ind+splay-1);
       if ((ind+splay)<sppanel)     neighbors.push_back(did*spdev+sec*sppanel+ind+splay);
     }
   else 
     {
       if ((ind-splay)>-1)          neighbors.push_back(did*spdev+sec*sppanel+ind-splay);
       if ((ind-splay+1)<splay)     neighbors.push_back(did*spdev+sec*sppanel+ind-splay+1);
       if ((ind-1)>(splay-1))       neighbors.push_back(did*spdev+sec*sppanel+ind-1);
       if ((ind+1)<sppanel)         neighbors.push_back(did*spdev+sec*sppanel+ind+1);
     }

   return neighbors;
*/ 
}//end of nListofNeighbors 

}


using mu2e::MakeStrawCluster;
DEFINE_FWK_MODULE(MakeStrawCluster);
