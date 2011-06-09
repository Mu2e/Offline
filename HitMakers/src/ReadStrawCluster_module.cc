//
// Plugin to test that I can read back the persistent data about straw hits.
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadStrawCluster_module.cc,v 1.16 2011/06/09 21:23:03 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/09 21:23:03 $
//
// Original author Hans Wenzel
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LineSegmentPCA.hh"
#include "Mu2eUtilities/inc/StrawClusterUtilities.hh"
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

  //--------------------------------------------------------------------
  //
  //
  class ReadStrawCluster : public art::EDAnalyzer {
  public:
    explicit ReadStrawCluster(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _clmakerModuleLabel(pset.get<std::string>("clmakerModuleLabel")),
      _hNInter(0),
      _hNClusters(0),
      _hNStraws(0),
      _R_rec(0),
      _x0y0(0)
    {
    }
    virtual ~ReadStrawCluster() { }

    virtual void beginJob();

    void analyze( art::Event const& e);
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

  void ReadStrawCluster::beginJob(){

    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
         << endl;

    art::ServiceHandle<art::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 500, 0., 500. );
    _hNStraws      = tfs->make<TH1F>( "hNStraws",  "Number of straws/cluster", 5  , 0., 5. );
    _R_rec         = tfs->make<TH1F>( "R_rec",  "reconstructed track radius", 100, 250., 350. );
    _x0y0          = tfs->make<TH2F>( "x0y0","x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);

  }

  void ReadStrawCluster::analyze(art::Event const& evt) {
    cout << "ReadStrawCluster: analyze() begin"<<endl;
    if ( _diagLevel > 2 ) cout << "ReadStrawCluster: analyze() begin"<<endl;
    static int ncalls(0);
    StrawClusterUtilities sutils=StrawClusterUtilities() ;

	
    multimap<int,StrawCluster> mpcluster;
    mpcluster.clear();
    art::Handle<StrawClusterCollection> pdataHandle;
    evt.getByLabel(_clmakerModuleLabel,pdataHandle);
    StrawClusterCollection const* clusters = pdataHandle.product();
    cout << "Nr of clusters:   " << clusters->size()<<endl;   
    for ( size_t cluster=0; cluster<clusters->size(); ++cluster) // Loop over StrawClusters
      {
	StrawCluster const& scluster = clusters->at(cluster);	
	cout << "averageT:  "<<sutils.averageT(scluster,evt)
	     <<"  did: "<< sutils.did(scluster,evt)
	     <<"  Station: "<<sutils.Station(scluster,evt)
	     <<"  sector: "<< sutils.secid(scluster,evt)
	     <<" X:  "<<sutils.X(scluster,evt)
	     <<" dirX:  "<<sutils.dirX(scluster,evt)
	     <<endl;
	mpcluster.insert(pair<int,StrawCluster>(sutils.did(scluster,evt),scluster));
      }
    cout <<"map size: " << mpcluster.size()<<endl;

    Int_t nint = 0;
    for (int i = 0;i<36;i++)
      {
	cout << " did:  "<<i << "  Count:  "<< mpcluster.count(i)<<endl;
	if (mpcluster.count(i)>1) 
	  {
	    pair<multimap<int,StrawCluster>::iterator, multimap<int,StrawCluster>::iterator> ppp1;
	    ppp1 = mpcluster.equal_range(i);
	    multimap<int,StrawCluster>::iterator first1 = ppp1.first;
	    multimap<int,StrawCluster>::iterator first2 = ppp1.first;
	    multimap<int,StrawCluster>::iterator last1 = ppp1.second;
	    last1--;
	    multimap<int,StrawCluster>::iterator last2 = ppp1.second;
	    for (first1;first1 != last1;++first1)
	      {
		first2=first1;
		first2++;
		for (first2;first2 != last2;++first2)
		  {
		  StrawCluster junk  = (*first1).second;
		  StrawCluster pjunk = (*first2).second;
		  LineSegmentPCA linesegment0 = sutils.linesegment(junk,evt);
		  LineSegmentPCA linesegment1 = sutils.linesegment(pjunk,evt);
		  CLHEP::Hep2Vector intersection;
		  switch(linesegment0.Intersect(linesegment1,intersection))
		    {
		    case LineSegmentPCA::PARALLEL:
		      //std::cout << "The lines are parallel\n\n";
		      break;
		    case LineSegmentPCA::COINCIDENT:
		      //std::cout << "The lines are coincident\n\n";
		      break;
		    case LineSegmentPCA::NOT_INTERSECTING:
		      // std::cout << "The lines do not intersect\n\n";
		      break;
		    case LineSegmentPCA::INTERSECTING:
		      std::cout << "The lines intersect at (" << intersection.x() << ", " << intersection.y() << ")\n\n";
		      //X.push_back(intersection.x_);
		      //Y.push_back(intersection.y_);
		      //Z.push_back(0.5*(junk.mpz+pjunk.mpz));
		      //R.push_back(sqrt(intersection.x_*intersection.x_ + intersection.y_+intersection.y_));
		      nint ++;
		      break;
		    }  // end switch 
		} // end for first2
	    }// end for first1
	}// end count >1
	  // cout << "Number of elements with key: "<<i<<"  " << m.count(i) << endl;
	  //pair<multimap< int,straw>::iterator, multimap<int,straw>::iterator> ppp;

      }   ///endloop over all devices
    cout<<"nint:  "<<nint<<endl;


















    ++ncalls;  
  } // end of ::analyze.
  void ReadStrawCluster::FitCircle(    vector<double> X,vector<double> Y)
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
    //        fitter->SetParameter(3, "omega",    175., 0.1, 0,0);
    //        fitter->SetParameter(4, "phase",    175., 0.1, 0,0);
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
    //                   fitter->GetParameter(1),fitter->GetParameter(2));
  //arc->SetLineColor(kRed);
  //arc->SetLineWidth(4);
  // arc->Draw();
  //int age;
  //cin >> age;
  }
}


using mu2e::ReadStrawCluster;
DEFINE_ART_MODULE(ReadStrawCluster);
