//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadDPIStrawCluster_plugin.cc,v 1.1 2011/01/28 21:38:12 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/28 21:38:12 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
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
#include "TH2F.h"
#include "TNtuple.h"
#include <TStyle.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>
#include <iostream>
#include <map>
#include <utility>
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TArc.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawClusterCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"

using namespace std;

namespace mu2e {
  TGraph *gr;
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = gr->GetN();
   f = 0;
   Double_t *x = gr->GetX();
   Double_t *y = gr->GetY();
   for (Int_t i=0;i<np;i++) {
      Double_t u = x[i] - par[0];
      Double_t v = y[i] - par[1];
      Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
      f += dr*dr;
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
      _clmakerModuleLabel(pset.getParameter<std::string>("clmakerModuleLabel")),
      _hNInter(0),
      _hNClusters(0),
      _hNStraws(0),
      _R_rec(0),
      _x0y0(0)
    {
    }
    virtual ~ReadDPIStrawCluster() { }
    
    virtual void beginJob(edm::EventSetup const&);
    
    void analyze( edm::Event const& e, edm::EventSetup const&);
    void FitCircle(    vector<double> X,vector<double> Y);
  private:
    
    // Diagnostics level.
    int _diagLevel;
    
    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    
    // Label of the module that made the Clusters.
    std::string _clmakerModuleLabel;
    // Some diagnostic histograms.
    TH1F* _hNInter;
    TH1F* _hNClusters;
    TH1F* _hNStraws;
    TH1F* _R_rec;
    TH2F* _x0y0;
  };
  
  void ReadDPIStrawCluster::beginJob(edm::EventSetup const& ){
    
    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;
    
    edm::Service<edm::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );  
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 500, 0., 500. );
    _hNStraws      = tfs->make<TH1F>( "hNStraws",  "Number of straws/cluster", 5  , 0., 5. );
    _R_rec         = tfs->make<TH1F>( "R_rec",  "reconstructed track radius", 100, 250., 350. );
    _x0y0          = tfs->make<TH2F>( "x0y0","x0 of circle vs y0 of circle ", 500,-650.,650.,500.,-650.,650.);

  }
  
  void ReadDPIStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&)
  {
    if ( _diagLevel > 2 ) cout << "ReadDPIStrawCluster: analyze() begin"<<endl;
    static int ncalls(0);
    ++ncalls;
    
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
    cout << "Length of pointers:  " << mcptrHandle->size()<<endl;
    _hNClusters->Fill(mcptrHandle->size());
    for ( size_t i=0; i< mcptrHandle->size(); ++i ) {
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));
       _hNStraws->Fill(mcptr.size());
       // cout << " Nr of Hits " << mcptr.size()<<endl;
      for( size_t j=0; j<mcptr.size(); ++j ) {
	DPIndex const& junkie = mcptr[j];
       	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
	//resolveDPIndex<StrawHitCollection>(evt,junkie);
	//cout << strawhit.strawIndex()<<endl; 	
	//   _hG4StepLength->Fill(mchit.stepLength());
	// _hG4StepEdep->Fill(mchit.eDep()*1000.0);
      }
    }

  } // end of ::analyze.
  void ReadDPIStrawCluster::FitCircle(    vector<double> X,vector<double> Y)
  {
    Int_t n = X.size();
    Double_t x[n];
    Double_t y[n];
    for ( size_t i=0; i<X.size(); ++i ) {
      x[i]=X[i];
      y[i]=Y[i];
    }
    gr = new TGraph(n,x,y);
    //   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",700,700);
    //c1->SetGrid();
    // c1->DrawFrame(-700.,-700.,700.,700.);

    //TEllipse el1(0.,0.,380.,380.);
    //TEllipse *ellipse = new TEllipse(0,0,680,680,0,360,0);
    //Int_t ci;   // for color index setting
    //Int_t ci = TColor::GetColor("#ffffcc");
    //ellipse->SetFillColor(ci);
    //ellipse->Draw();
    
    //el1.Draw("SAME");
    //gr->SetLineColor(2);
    //gr->SetLineWidth(4);
    //gr->SetMarkerColor(4);
    //gr->SetMarkerStyle(21);
    //gr->SetTitle("x-y projection");
    //gr->GetXaxis()->SetTitle("x [mm]");
    //gr->GetYaxis()->SetTitle("y [mm]");
    //gr->Draw("P");
    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    fitter->SetFCN(myfcn);
    fitter->SetParameter(0, "x0",   0, 0.1, 0,0);
    fitter->SetParameter(1, "y0",   0, 0.1, 0,0);
    fitter->SetParameter(2, "R",    175., 0.1, 0,0);
    //	      fitter->SetParameter(3, "omega",    175., 0.1, 0,0);
    //	      fitter->SetParameter(4, "phase",    175., 0.1, 0,0);
    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
    cout << "x0:   " << fitter->GetParameter(0)
	 << " y0:  " << fitter->GetParameter(1)
	 << " r:   " << fitter->GetParameter(2)
    <<endl;
    _x0y0->Fill(fitter->GetParameter(0),fitter->GetParameter(1));
    _R_rec->Fill(fitter->GetParameter(2));
    //Draw the circle on top of the points
    //TArc *arc = new TArc(fitter->GetParameter(0),
    //			 fitter->GetParameter(1),fitter->GetParameter(2));
  //arc->SetLineColor(kRed);
  //arc->SetLineWidth(4);
  // arc->Draw();
  //int age;
  //cin >> age;
  }
}


using mu2e::ReadDPIStrawCluster;
DEFINE_FWK_MODULE(ReadDPIStrawCluster);
