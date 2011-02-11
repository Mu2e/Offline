//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadDPIStrawCluster_plugin.cc,v 1.5 2011/02/11 00:06:20 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/02/11 00:06:20 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <utility>

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
#include "TH2F.h"
#include "TF1.h"
#include "TNtuple.h"
#include <TStyle.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TArc.h"
#include "TMinuit.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "Mu2eG4/inc/ConvElecUtilities.hh" 
#include "ToyDP/inc/StatusG4.hh"

using namespace std;

namespace mu2e {
  enum PrintLevel { quiet  =-1,
		    normal = 0,
		    verbose= 1};
  
  TGraph *gr2;
  TGraphErrors *error;
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = error->GetN();
   f = 0;
   Double_t *x = error->GetX();
   Double_t *y = error->GetY();
   Double_t *ex = error->GetEX();
   Double_t *ey = error->GetEY();

   for (Int_t i=0;i<np;i++) {
      Double_t u = (x[i] - par[0]);
      Double_t v = (y[i] - par[1]);
      Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
      //f += (dr*dr)/(ex[i]*ex[i]+ey[i]*ey[i]);
      f += (dr*dr)/(0.25*ex[i]*ex[i]);
   }
}
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
 class pstraw{
   //
   // pseudo straw class 
   //
 public:
   Int_t   lay;
   Int_t   did;
   Int_t   sec;
   Float_t hl;
   Float_t mpx;
   Float_t mpy;
   Float_t mpz;
   Float_t dirx;
   Float_t diry;
   Float_t dirz; // should always be 0 
   /*
    bool operator>(const pstraw other) const {
      if (id > other.id) {
	return true;
      }
      else{
	return false;
      }
    }
   bool operator<(const pstraw other) const {
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
   */
    void Print()
    {
      cout<< "Straw:  " << endl;
      cout<< "======  " << endl;
      cout<< "Layer:  " << lay  <<endl;
      cout<< "DID:    " << did  <<endl;
      cout<< "Sector: " << sec  <<endl;
      cout<< "hl:     " << hl   <<endl;
      cout<< "mpx:    " << mpx  <<endl;
      cout<< "mpy:    " << mpy  <<endl;
      cout<< "mpz:    " << mpz  <<endl;
      cout<< "dirx:   " << dirx <<endl;
      cout<< "diry:   " << diry <<endl;
      cout<< "dirz:   " << dirz <<endl;     
    }
 };

  //--------------------------------------------------------------------
  //
  // 
  class ReadDPIStrawCluster : public edm::EDAnalyzer {
  public:
    explicit ReadDPIStrawCluster(edm::ParameterSet const& pset):
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
      _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
      _clmakerModuleLabel(pset.getParameter<std::string>("clmakerModuleLabel")),
      _hNInter(0),
      _hNClusters(0),
      _hNHits(0),
      _hNEleHits(0),
      _hNStraws(0),
      _R_rec(0),
      _Pt_in(0),
      _P_in(0),
      _Pz_in(0),
      _Pt_rec(0),
      _P_rec(0),
      _Pz_rec(0),
      _Pt_diff(0),
      _P_diff(0),
      _Pz_diff(0),
      _x0y0(0),
      _chi2(0)
    {
    }
    virtual ~ReadDPIStrawCluster() { }
    
    virtual void beginJob(edm::EventSetup const&);
    
    void analyze( edm::Event const& e, edm::EventSetup const&);
    void FitCircle(vector<double> X,vector<double> Y);
    void FitSinus( vector<double> R,vector<double> Z);

  private:
 
    // Diagnostics level.
    int _diagLevel;
    
    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
        
    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the Clusters.
    std::string _clmakerModuleLabel;

    // Some diagnostic histograms.
    TH1F* _hNInter;
    TH1F* _hNClusters;
    TH1F* _hNHits;
    TH1F* _hNEleHits;
    TH1F* _hNStraws;
    TH1F* _R_rec;
    TH1F* _Pt_in;
    TH1F* _P_in;
    TH1F* _Pz_in;
    TH1F* _Pt_rec;
    TH1F* _P_rec;
    TH1F* _Pz_rec;
    TH1F* _Pt_diff;
    TH1F* _P_diff;
    TH1F* _Pz_diff;
    TH2F* _x0y0;
    TH1F* _chi2;
    Double_t R_rec,x0,y0;
    Double_t Pt,Pz;
    CLHEP::Hep3Vector  X_in;  
    CLHEP::Hep3Vector  P_in; 
    Double_t Pt_inval;
    Double_t P_inval;
  };
  
  void ReadDPIStrawCluster::beginJob(edm::EventSetup const& ){
    
    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;
    
    edm::Service<edm::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );  
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 500, 0., 500. );
    _hNHits        = tfs->make<TH1F>( "hNHits","Number of straw Hits", 500, 0., 500. );
    _hNEleHits     = tfs->make<TH1F>( "hNEleHits","Number of straw Hits/conversion electron", 500, 0., 500. );
    _hNStraws      = tfs->make<TH1F>( "hNStraws",  "Number of straws/cluster", 5  , 0., 5. );
    _R_rec         = tfs->make<TH1F>( "R_rec",  "reconstructed track radius", 100, 0., 800. );
    _Pt_in         = tfs->make<TH1F>( "Pt_in",  "tansverse momentum Pt at first Hit", 100, 0., 160. );
    _P_in          = tfs->make<TH1F>( "P_in",  "momentum at first Hit", 100, 60., 140. );
    _Pz_in         = tfs->make<TH1F>( "Pz_in",  "longitudinal momentum at first Hit", 100, 0., 120. );
    _Pt_rec        = tfs->make<TH1F>( "Pt_rec",  "reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec         = tfs->make<TH1F>( "P_rec",  "reconstructed momentum", 100, 60., 140. );
    _Pz_rec        = tfs->make<TH1F>( "Pz_rec",  "reconstructed longitudinal momentum", 100, 0., 120. );
    _Pt_diff       = tfs->make<TH1F>( "Pt_diff",  "delta tansverse momentum Pt", 100, -20., 20. );
    _P_diff        = tfs->make<TH1F>( "P_diff",  "delta momentum", 100, -20., 20.);
    _Pz_diff       = tfs->make<TH1F>( "Pz_diff",  "delta longitudinal momentum", 100, -20., 20.);
    _x0y0          = tfs->make<TH2F>( "x0y0","x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2          = tfs->make<TH1F>( "chi2",  "chi2", 200, 0., 200. );
  }
  
  void ReadDPIStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&)
  {
    if ( _diagLevel > 2 ) cout << "ReadDPIStrawCluster: analyze() begin"<<endl;
    static int ncalls(0);
    ++ncalls;
    StrawId nsid;
    Straw str;
    StrawId sid;
    LayerId lid;
    DeviceId did;
    SectorId secid;

// calculate the average hit position of track at a plane
    double Xavg[36], Yavg[36], edep[36] ;
    int nhitdev[36];
    
    for (int i = 0; i < 36 ; i++) { 
      nhitdev[i] = 0 ;
      Xavg[i] = 0 ;
      Yavg[i] = 0 ;
      edep[i] = 0 ;
    } 
    
    // Geometry info for the TTracker.
    // Get a reference to one of the T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();
    multimap<int,pstraw> mpstraws;
    mpstraws.clear();

    //
    // Get the persistent data about pointers to StrawHits
    //
    edm::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_clmakerModuleLabel,"DPIStrawCluster",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    edm::Handle<StatusG4> g4StatusHandle;
    evt.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    evt.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    // Get handles to the generated and simulated particles.
    edm::Handle<ToyGenParticleCollection> genParticles;
    evt.getByType(genParticles);
    
    edm::Handle<SimParticleCollection> simParticles;
    evt.getByType(simParticles);
    
    // Handle to information about G4 physical volumes.
    edm::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByType(volumes);
    
    //Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    cout<<"_g4ModuleLabel:  "<<_g4ModuleLabel
	<<"  trackerStepPoints:  "<<_trackerStepPoints
	<<"  Hits:  " <<hits->size()
	<<endl;
    //Class constructor: you have to pass the event and two strings, one for the g4module, one for the tracker name.
    //The constructor checks that we have hits in the event
    ConvElecUtilities ConvElec(evt, _g4ModuleLabel, _trackerStepPoints); 
    
    //This gives you the number of hits associated to the generated conversion electron 
    int nCEHits = ConvElec.hasStepPointMC(); 
    
    //Here I check that there is one conversion electron
    if (nCEHits >=1) 
      {
	//This gives a reference to the first hit (temporally speaking) associated to the converse electron
	const StepPointMC & CEFirstHit = ConvElec.firstHit(); 
	
	//Now you can ask position, momentum and whatever to this hit, looking at ToyDP/inc/StepPointMC.hh
	cout << "Position of first Hit:     "<<CEFirstHit.position()   <<endl;   
	cout << "Momentum at first Hit:     "<<CEFirstHit.momentum()   <<endl;  
	cout << "Staw Index of first Hit:   "<<CEFirstHit.strawIndex() <<endl;
	X_in= CEFirstHit.position();  
	P_in= CEFirstHit.momentum();
	Pt_inval =  TMath::Sqrt(P_in.x()*P_in.x()+P_in.y()*P_in.y());
	P_inval  =  TMath::Sqrt(P_in.x()*P_in.x()+P_in.y()*P_in.y()+P_in.z()*P_in.z());
	_Pt_in->Fill( Pt_inval);
	_P_in->Fill(P_inval);
	_Pz_in->Fill(P_in.z());
	//The following line gives you the StrawIndex of the straw associated to the first hit of the conversion electron 
	StrawIndex Idx = ConvElec.earliestStrawIndex();
	
	//Using the StrawIndex above you can point to the straw, passing through the tracker
	const Straw& CEFirstStraw = tracker.getStraw(Idx);
	
	//You can ask the straw its Id or whatever, just looking at   ToyDP/inc/Straw.hh
	//cout << CEFirstStraw.Id() << endl; 
	
	//This method gives a vector of StrawIndex associated to all the straws hit by the conversion electron
	vector<StrawIndex> CEStrawsIdx = ConvElec.convElecStrawIdx(); 
	
	//You can loop through the conversion electron straws, passing the element of the vector (StrawIndex)
	// to the tracker and ask them whatever you prefer
	_hNEleHits->Fill( CEStrawsIdx.size());
	//for (vector<StrawIndex>::iterator i = CEStrawsIdx.begin(); 
	//     i != CEStrawsIdx.end(); ++i) { 
	// cout << tracker.getStraw(*i).Id() << endl; 
	//}
	vector<size_t> CEStepPoints =ConvElec.convElecHitsIdx();
	for (vector<size_t>::iterator i = CEStepPoints.begin(); 
	     i != CEStepPoints.end(); ++i) { 
	  const StepPointMC hit = (*hits)[*i];
	  cout << "Step point index:  "<<*i << endl;
	  const CLHEP::Hep3Vector& pos = hit.position();
	  const CLHEP::Hep3Vector& mom = hit.momentum();
	  cout << "Position of hits " << hit.position() << endl ;      
	  cout << "Momentum of hits " << hit.momentum() << endl ;  
	  // find out what device or plane we are in (layer goes from 0 to 35) 
      
	  const Straw&      straw = tracker.getStraw( hit.strawIndex() );
	  StrawId nsid = straw.Id() ;	  
	  nhitdev[nsid.getDevice()]++;
	  Xavg[nsid.getDevice()] = Xavg[nsid.getDevice()] + pos.x() ;
	  Yavg[nsid.getDevice()] = Yavg[nsid.getDevice()] + pos.y() ;
	  edep[nsid.getDevice()] = edep[nsid.getDevice()] + hit.eDep() ;  
	}
 
      }
    else
      {
	cout << "No hits associated to Conversion Electron" << endl; 
      } 
    //  calculate the average hit position of the MC hit within each layer
    for (int i = 0; i < 36 ; i++) { 
      if (nhitdev[i] != 0) {
	Xavg[i] = Xavg[i]/nhitdev[i]; 
	Yavg[i] = Yavg[i]/nhitdev[i]; 
	//	 cout << i << " " << nhitdev[i] << " " << Xavg[i] << " " << Yavg[i] << endl ; 
	//_eplane -> Fill(edep[i]) ;  // energy deposited in plane
	//_nhitpl -> Fill(nhitdev[i]) ; // number of hits in plane
      }
    } 
    
    vector<double> X;
    vector<double> Y;
    vector<double> R;
    vector<double> Z;
    X.clear();
    Y.clear();
    R.clear();
    Z.clear();
    Int_t totalHits=0;
    _hNClusters->Fill(mcptrHandle->size());
    CLHEP::Hep3Vector dvec;
    for ( size_t i=0; i< mcptrHandle->size(); ++i ) {
      double hlen=9999999.;
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));
      CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
      _hNStraws->Fill(mcptr.size());
      totalHits = totalHits+mcptr.size();
      for( size_t j=0; j<mcptr.size(); j++ ) {
	DPIndex const& junkie = mcptr[j];
       	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
	str = tracker.getStraw(strawhit.strawIndex());
	sid = str.Id();
	lid = sid.getLayerId();
	did = sid.getDeviceId();
	secid = sid.getSectorId();
	const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
	const CLHEP::Hep3Vector dirvec = str.getDirection();
	dvec = CLHEP::Hep3Vector(dirvec.getX(),dirvec.getY(),dirvec.getZ());
	pvec = pvec + mpvec;
	if (str.getHalfLength()<hlen)
	  {
	    hlen=str.getHalfLength();
	  }
      }
      double a = 1./double(mcptr.size());
      pvec = pvec*a;
      pstraw pstr;
      pstr.lay=lid.getLayer();
      pstr.did=did;
      pstr.sec=secid.getSector();
      pstr.hl=hlen;
      pstr.mpx=pvec.getX();
      pstr.mpy=pvec.getY();
      pstr.mpz=pvec.getZ();
      pstr.dirx=dvec.getX();
      pstr.diry=dvec.getY();
      pstr.dirz=dvec.getZ();
      mpstraws.insert(pair<int,pstraw>(did,pstr));
    }
    _hNHits->Fill(totalHits);
    //cout << " size of pseudo straw map: " <<mpstraws.size()<<endl; 

    Int_t nint = 0;
    for (int i = 0;i<36;i++)
      {
	if (mpstraws.count(i)>1) 
	  {
	    pair<multimap<int,pstraw>::iterator, multimap<int,pstraw>::iterator> ppp1;
	    ppp1 = mpstraws.equal_range(i);
	    multimap<int,pstraw>::iterator first1 = ppp1.first;
	    multimap<int,pstraw>::iterator first2 = ppp1.first;
	    multimap<int,pstraw>::iterator last1 = ppp1.second;
	    last1--;
	    multimap<int,pstraw>::iterator last2 = ppp1.second;
	    for (first1;first1 != last1;++first1)
	      {
		first2=first1;
		first2++;
		for (first2;first2 != last2;++first2)
		  {
		    pstraw junk  = (*first1).second;
		    pstraw pjunk = (*first2).second;
		    const Vector p0= Vector(junk.mpx-junk.hl*junk.dirx,junk.mpy-junk.hl*junk.diry);
		    const Vector p1= Vector(junk.mpx+junk.hl*junk.dirx,junk.mpy+junk.hl*junk.diry);
		    const Vector p2= Vector(pjunk.mpx-pjunk.hl*pjunk.dirx,pjunk.mpy-pjunk.hl*pjunk.diry); 
		    const Vector p3= Vector(pjunk.mpx+pjunk.hl*pjunk.dirx,pjunk.mpy+pjunk.hl*pjunk.diry);
		    LineSegment linesegment0(p0, p1);
		    LineSegment linesegment1(p2, p3);
		    Vector intersection;
		    switch(linesegment0.Intersect(linesegment1, intersection))
		      {
		      case LineSegment::PARALLEL:
			//std::cout << "The lines are parallel\n\n";
			break;
		      case LineSegment::COINCIDENT:
			//std::cout << "The lines are coincident\n\n";
			break;
		      case LineSegment::NOT_INTERSECTING:
			//std::cout << "The lines do not intersect\n\n";
			break;
		      case LineSegment::INTERSECTING:
			X.push_back(intersection.x_);
			Y.push_back(intersection.y_);
			Z.push_back(0.5*(junk.mpz+pjunk.mpz));
			R.push_back( TMath::Sqrt(intersection.x_*intersection.x_ + intersection.y_*intersection.y_));

			nint ++;
			//cout<<nint<<endl;
			break;
		      }  // end switch 
		  } // end for first2
	      }// end for first1
	  }// end count >1
      }   ///endloop over all devices
    _hNInter->Fill(X.size());
    if (X.size()>2)
      {
	FitCircle(X, Y);
	Double_t Bmagnet=10.;   // 10 KGauss magnetic field
	Double_t Const=1.49898e-4;
	Pt =1000.*R_rec*0.1 * 2. * Bmagnet* Const;
	_Pt_rec->Fill(Pt);
	_Pt_diff->Fill(Pt-Pt_inval);
	FitSinus(R, Z);
	_Pz_rec->Fill(Pz);
	_Pz_diff->Fill(Pz-P_in.z());
	Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
	_P_rec->Fill(Ptot);
	_P_diff->Fill(Ptot-P_inval);

      }
  } // end of ::analyze.

  void ReadDPIStrawCluster::FitSinus(    vector<double> R,     vector<double> Z)
  {
    Int_t n = R.size();
    Double_t z[n];
    Double_t r[n];
    for ( size_t i=0; i<R.size(); ++i ) {
      z[i]=Z[i];
      r[i]=R[i];
    }
    TVirtualFitter::SetDefaultFitter("Minuit");
    gMinuit->SetPrintLevel(quiet);
    gr2 = new TGraph(n,z,r);
    TF1 *f2 = new TF1("f2", "[0]+[1]*sin([2]*x-[3])", -2000., 2000.);
    Double_t offset =  TMath::Sqrt(x0*x0+y0*y0);
    Double_t radius= R_rec;
    gMinuit->SetPrintLevel(quiet);
    f2->SetParameters(offset,radius,0.005,0.1);
    f2->FixParameter(0,offset);
    f2->FixParameter(1,radius);
    gr2->Fit(f2);
    Double_t p2 = f2->GetParameter(2);
    Pz = 10./(33.36*p2);
    
  }
  void ReadDPIStrawCluster::FitCircle(    vector<double> X,vector<double> Y)
  {
    Int_t n = X.size();
    Double_t x[n];
    Double_t y[n];
    Double_t ex[n] ; Double_t ey[n] ;
    for ( size_t i=0; i<X.size(); ++i ) {
      x[i]=X[i];
      y[i]=Y[i];
      ex[i] = 5.0 ; 
      ey[i] = 5.0 ;
    }
    error = new TGraphErrors(n,x,y,ex,ey);
    TMinuit *gmMinuit = new TMinuit(3); 
    gmMinuit->SetPrintLevel(quiet);
    gmMinuit->SetFCN(myfcn);
    const int dim(3);
    const char par_name[dim][20]={"x0","y0","R"};
    static Double_t step[dim] = {0.001,0.001,0.001};
    Double_t sfpar[dim]={0.0,0.0,175.};
    Double_t errsfpar[dim]={0.0,0.0,0.0};
    int ierflg = 0;
    for (int ii = 0; ii<dim; ii++) {    
      gmMinuit->mnparm(ii,par_name[ii],sfpar[ii], step[ii], 0,0,ierflg);
    }
    int result=gmMinuit->Migrad();
    bool converged = gmMinuit->fCstatu.Contains("CONVERGED");
    if (!converged) 
      {
	cout <<"-----------fit didn't converge---------------------------" <<endl;
	return;
      }
    for (int i = 0;i<3;i++) {
     gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
    } 
    x0    = sfpar[0];
    y0    = sfpar[1];
    R_rec = sfpar[2];
    Double_t chi2, edm, errdef; 
    Int_t nvpar, nparx,istat;
    gmMinuit->mnstat(chi2,edm,errdef,nvpar,nparx,istat);
    _x0y0->Fill(x0,y0);
    _R_rec->Fill(R_rec);
    _chi2 -> Fill(chi2) ;
  }


}


using mu2e::ReadDPIStrawCluster;
DEFINE_FWK_MODULE(ReadDPIStrawCluster);
