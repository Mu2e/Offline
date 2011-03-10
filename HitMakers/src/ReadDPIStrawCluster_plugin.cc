//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadDPIStrawCluster_plugin.cc,v 1.12 2011/03/10 02:07:37 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/03/10 02:07:37 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <utility>

#include "CLHEP/Vector/TwoVector.h"

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
#include "TStyle.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TArc.h"
#include "TMinuit.h"
#include "TGeoHelix.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "Mu2eG4/inc/ConvElecUtilities.hh" 
#include "ToyDP/inc/StatusG4.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "ToyDP/inc/ToyGenParticle.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/GenId.hh"

using namespace std;

namespace mu2e {;
  enum PrintLevel { quiet  =-1,
		    normal = 0,
		    verbose= 1};
  enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERSECTING, INTERSECTING };
  Double_t Radius;  
  Double_t curv;
  Double_t zstep;
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

  
  class LineSegment
  {
  private:
    CLHEP::Hep2Vector begin_;
    CLHEP::Hep2Vector end_;
  public:
    LineSegment(const CLHEP::Hep2Vector& begin, const CLHEP::Hep2Vector& end)
      : begin_(begin), end_(end) {}
    enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERSECTING, INTERSECTING };
    
    IntersectResult Intersect(const LineSegment& other_line, CLHEP::Hep2Vector& intersection)
    {
      float denom = ((other_line.end_.y() - other_line.begin_.y())*(end_.x() - begin_.x())) -
	((other_line.end_.x() - other_line.begin_.x())*(end_.y() - begin_.y()));
      
      float nume_a = ((other_line.end_.x() - other_line.begin_.x())*(begin_.y() - other_line.begin_.y())) -
	((other_line.end_.y() - other_line.begin_.y())*(begin_.x() - other_line.begin_.x()));
      
      float nume_b = ((end_.x() - begin_.x())*(begin_.y() - other_line.begin_.y())) -
	((end_.y() - begin_.y())*(begin_.x() - other_line.begin_.x()));
      
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
	  intersection =CLHEP::Hep2Vector(begin_.x() + ua*(end_.x() - begin_.x()),begin_.y() + ua*(end_.y() - begin_.y()));
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
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel")),
      _clmakerModuleLabel(pset.getParameter<std::string>("clmakerModuleLabel")),
      _hNInter(0),
      _hNClusters(0),
      _hNHits(0),
      _hNEleHits(0),
      _hNStraws(0),
      _Polar(0),
      _R_rec(0),
      _Pt_in_si(0),
      _Pz_in_si(0),
      _P_in_si(0),
      _Pt_out_si(0),
      _Pz_out_si(0),
      _P_out_si(0),
      _Pt_rec(0),
      _P_rec(0),
      _Pz_rec(0),
      _Pt_diff(0),
      _P_diff(0),
      _Pz_diff(0),
      _x0y0(0),
      _chi2(0),
      _Xdiff(0),
      _Ydiff(0),
      _Rdiff(0),
      _Phidiff(0),
      _R_rec_s(0),
      _Pt_rec_s(0),
      _P_rec_s(0),
      _Pz_rec_s(0),
      _Pt_diff_s(0),
      _P_diff_s(0),
      _Pz_diff_s(0),
      _x0y0_s(0),
      _chi2_s(0),
      _Xdiff_s(0),
      _Ydiff_s(0),
      _Rdiff_s(0),
      _Phidiff_s(0),
      _R_rec_c(0),
      _Pt_rec_c(0),
      _P_rec_c(0),
      _Pz_rec_c(0),
      _Pt_diff_c(0),
      _P_diff_c(0),
      _Pz_diff_c(0),
      _x0y0_c(0),
      _chi2_c(0),
      _Xdiff_c(0),
      _Ydiff_c(0),
      _Rdiff_c(0),
      _Phidiff_c(0),
      _eplane(0)
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
    // Label of the module that made the StrawHits.
    std::string _makerModuleLabel;
    // Label of the module that made the Clusters.
    std::string _clmakerModuleLabel;
    
    // Some diagnostic histograms.
    TH1F* _hNInter;
    TH1F* _hNClusters;
    TH1F* _hNHits;
    TH1F* _hNEleHits;
    TH1F* _hNStraws;
    TH1F*_Polar;
    TH1F* _R_rec;
    //
    TH1F* _Pt_in_si;
    TH1F* _Pz_in_si;
    TH1F* _P_in_si;
    //
    TH1F* _Pt_out_si;
    TH1F* _Pz_out_si;
    TH1F* _P_out_si;
    //
    TH1F* _Pt_rec;
    TH1F* _P_rec;
    TH1F* _Pz_rec;
    TH1F* _Pt_diff;
    TH1F* _P_diff;
    TH1F* _Pz_diff;
    TH2F* _x0y0;
    TH1F* _chi2;
    TH1F* _Xdiff;
    TH1F* _Ydiff;
    TH1F* _Rdiff;
    TH1F* _Phidiff;
    //
    TH1F* _R_rec_s;
    TH1F* _Pt_rec_s;
    TH1F* _P_rec_s;
    TH1F* _Pz_rec_s;
    TH1F* _Pt_diff_s;
    TH1F* _P_diff_s;
    TH1F* _Pz_diff_s;
    TH2F* _x0y0_s;
    TH1F* _chi2_s;
    TH1F* _Xdiff_s;
    TH1F* _Ydiff_s;
    TH1F* _Rdiff_s;
    TH1F* _Phidiff_s;
    //
    TH1F* _R_rec_c;
    TH1F* _Pt_rec_c;
    TH1F* _P_rec_c;
    TH1F* _Pz_rec_c;
    TH1F* _Pt_diff_c;
    TH1F* _P_diff_c;
    TH1F* _Pz_diff_c;
    TH2F* _x0y0_c;
    TH1F* _chi2_c;
    TH1F* _Xdiff_c;
    TH1F* _Ydiff_c;
    TH1F* _Rdiff_c;
    TH1F* _Phidiff_c;
    //
    TH1F* _eplane;
    Double_t R_rec,x0,y0,chi2;
    Double_t Pt,Pz;
    CLHEP::Hep3Vector  X_in;  
    CLHEP::Hep3Vector  P_in_si;
    CLHEP::Hep3Vector  P_out_si;
    //
    Double_t Pt_inval_si;
    Double_t P_inval_si;
    Double_t Pt_outval_si;
    Double_t P_outval_si;
    //
    CLHEP::Hep3Vector B;
   };
  
  void ReadDPIStrawCluster::beginJob(edm::EventSetup const& ){
    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint<<endl; 

    edm::Service<edm::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );  
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 500, 0., 500. );
    _hNHits        = tfs->make<TH1F>( "hNHits",    "Number of straw Hits", 500, 0., 500. );
    _hNEleHits     = tfs->make<TH1F>( "hNEleHits", "Number of straw Hits/conversion electron", 500, 0., 500. );
    _hNStraws      = tfs->make<TH1F>( "hNStraws",  "Number of straws/cluster", 5  , 0., 5. );
    _Polar         = tfs->make<TH1F>( "Polar",     "polar angle at production", 100, -1, 1. );
    _R_rec         = tfs->make<TH1F>( "R_rec",     "reconstructed track radius", 100, 0., 800. );
    //
    _Pt_in_si      = tfs->make<TH1F>( "Pt_in_si",     "tansverse momentum Pt_si at first Hit", 100, 0., 120. );
    _Pz_in_si      = tfs->make<TH1F>( "Pz_in_si",     "longitudinal momentum Pz_si at first Hit", 100, 0., 100. );
    _P_in_si       = tfs->make<TH1F>( "P_in_si",      "momentum at first Hit P_si", 100, 0., 110. );
    //
    _Pt_out_si     = tfs->make<TH1F>( "Pt_out_si",     "tansverse momentum Pt_out_si at last Hit", 100, 0., 120. );
    _Pz_out_si     = tfs->make<TH1F>( "Pz_out_si",     "longitudinal momentum Pz_out_si at last Hit", 100, 0., 100. );
    _P_out_si      = tfs->make<TH1F>( "P_out_si",      "momentum P_out_si at last Hit", 100, 0., 110. );
    //
    _Pt_rec        = tfs->make<TH1F>( "Pt_rec",    "reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec         = tfs->make<TH1F>( "P_rec",     "reconstructed momentum", 100, 60., 140. );
    _Pz_rec        = tfs->make<TH1F>( "Pz_rec",    "reconstructed longitudinal momentum", 100, 0., 120. );
    _Pt_diff       = tfs->make<TH1F>( "Pt_diff",   "delta tansverse momentum Pt", 100, -20., 20. );
    _P_diff        = tfs->make<TH1F>( "P_diff",    "delta momentum", 100, -20., 20.);
    _Pz_diff       = tfs->make<TH1F>( "Pz_diff",   "delta longitudinal momentum", 100, -20., 20.);
    _x0y0          = tfs->make<TH2F>( "x0y0",      "x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2          = tfs->make<TH1F>( "chi2",      "chi2 of 2d circle fit", 200, 0., 200. );
    _Xdiff         = tfs->make<TH1F>( "Xdiff",     "X(stepMC)  - X(stereo hit) mm ", 100, -20., 20.);
    _Ydiff         = tfs->make<TH1F>( "Ydiff",     "Y(stepMC)  - Y(stereo hit)  mm ", 100, -20., 20.);
    _Rdiff         = tfs->make<TH1F>( "Rdiff",     "R(stepMC)  - R(stereo hit) mm ", 100, -100., 20.);
    _Phidiff       = tfs->make<TH1F>( "Phidiff",   "Phi(stepMC)- Phi(stereo hit) mm ", 100, -0.6, 0.6);
    //
    _R_rec_s         = tfs->make<TH1F>( "R_rec_s",     "s reconstructed track radius", 100, 0., 800. );
    _Pt_rec_s        = tfs->make<TH1F>( "Pt_rec_s",    "s reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec_s         = tfs->make<TH1F>( "P_rec_s",     "s reconstructed momentum", 100, 60., 140. );
    _Pz_rec_s        = tfs->make<TH1F>( "Pz_rec_s",    "s reconstructed longitudinal momentum", 100, 0., 120. );
    _Pt_diff_s       = tfs->make<TH1F>( "Pt_diff_s",   "s delta tansverse momentum Pt", 100, -20., 20. );
    _P_diff_s        = tfs->make<TH1F>( "P_diff_s",    "s delta momentum", 100, -20., 20.);
    _Pz_diff_s       = tfs->make<TH1F>( "Pz_diff_s",   "s delta longitudinal momentum", 100, -20., 20.);
    _x0y0_s          = tfs->make<TH2F>( "x0y0_s",      "s x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2_s          = tfs->make<TH1F>( "chi2_s",      "s chi2 of 2d circle fit", 200, 0., 200. );
    _Xdiff_s         = tfs->make<TH1F>( "Xdiff_s",     "s X(stepMC)  - X(straw hit) mm ", 100, -20., 20.);
    _Ydiff_s         = tfs->make<TH1F>( "Ydiff_s",     "s Y(stepMC)  - Y(straw hit)  mm ", 100, -20., 20.);
    _Rdiff_s         = tfs->make<TH1F>( "Rdiff_s",     "s R(stepMC)  - R(straw hit) mm ", 100, -10., 10.);
    _Phidiff_s       = tfs->make<TH1F>( "Phidiff_s",   "s Phi(stepMC)- Phi(straw hit) mm ", 100, -0.05, 0.05);
    //
    //
    _R_rec_c         = tfs->make<TH1F>( "R_rec_c",     "s reconstructed track radius", 100, 0., 800. );
    _Pt_rec_c        = tfs->make<TH1F>( "Pt_rec_c",    "s reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec_c         = tfs->make<TH1F>( "P_rec_c",     "s reconstructed momentum", 100, 60., 140. );
    _Pz_rec_c        = tfs->make<TH1F>( "Pz_rec_c",    "s reconstructed longitudinal momentum", 100, 0., 120. );
    _Pt_diff_c       = tfs->make<TH1F>( "Pt_diff_c",   "s delta tansverse momentum Pt", 100, -20., 20. );
    _P_diff_c        = tfs->make<TH1F>( "P_diff_c",    "s delta momentum", 100, -20., 20.);
    _Pz_diff_c       = tfs->make<TH1F>( "Pz_diff_c",   "s delta longitudinal momentum", 100, -20., 20.);
    _x0y0_c          = tfs->make<TH2F>( "x0y0_c",      "s x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2_c          = tfs->make<TH1F>( "chi2_c",      "s chi2 of 2d circle fit", 200, 0., 200. );
    _Xdiff_c         = tfs->make<TH1F>( "Xdiff_c",     "s X(stepMC)  - X(stereo hit) mm ", 100, -20., 20.);
    _Ydiff_c         = tfs->make<TH1F>( "Ydiff_c",     "s Y(stepMC)  - Y(stereo hit)  mm ", 100, -20., 20.);
    _Rdiff_c         = tfs->make<TH1F>( "Rdiff_c",     "s R(stepMC)  - R(stereo hit) mm ", 100, -50., 50.);
    _Phidiff_c       = tfs->make<TH1F>( "Phidiff_c",   "s Phi(stepMC)- Phi(stereo hit) mm ", 100, -0.6, 0.6);
    //
    _eplane        = tfs->make<TH1F>( "eplane",    "edep in plane (MeV)", 200, 0., .1 );
  }
  
  void ReadDPIStrawCluster::analyze(edm::Event const& evt, edm::EventSetup const&)
  {
    if ( _diagLevel > 2 ) cout << "ReadDPIStrawCluster: analyze() begin"<<endl;
    // Geometry info for the TTracker.
    // Get a reference to one of the T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();
    //
    // Get the magnetic field in the tracker: 
    // Note there are some issues with the magnetic field units 
    // that's why the two methods yield different results.
    // 
    GeomHandle<BFieldManager> bfMgr;
    B= bfMgr->getDSUniformValue(); 
    cout << " B-field:  " <<B<<endl;
    //
    // Position of the center of the tracker in mu2e coordinates.
    //
    CLHEP::Hep2Vector tt =CLHEP::Hep2Vector( 0.0, 10200.);
    CLHEP::Hep3Vector point =CLHEP::Hep3Vector( -3904, 0.0, 10200.);
    CLHEP::Hep3Vector bf = bfMgr->getBField(point); 
    cout << " B-field:  " <<bf<<endl;
    static int ncalls(0);
    ++ncalls;
    StrawId nsid;
    Straw str;
    StrawId sid;
    LayerId lid;
    DeviceId did;
    SectorId secid;
    vector<double> X;       // x of cluster intersections
    vector<double> Y;       // y of cluster intersections
    vector<double> Z;       // z of cluster intersections
    vector<double> X_res;   // x residuals
    vector<double> Y_res;   // y residual 
    vector<double> R_res;   // R residuals
    vector<double> Phi_res; // phi residual 
    vector<double> R;       // radius of cluster intersections
    vector<double> Phi;     // angle of cluster intersections
    
    vector<CLHEP::Hep3Vector> Points3d; 
    //
    vector<double> X_straw;       // x of dt straw
    vector<double> Y_straw; 
    vector<double> Z_straw; 
    vector<double> X_res_straw;   // x residuals
    vector<double> Y_res_straw;   // y residual 
    vector<double> R_res_straw;   // R residuals
    vector<double> Phi_res_straw; // phi residual 
    vector<double> R_straw;       // radius
    vector<double> Phi_straw;     // angle 
    vector<CLHEP::Hep3Vector> Points3d_straw; // x,y measurement packed in

    vector<double> X_cluster;       // x of dt straw
    vector<double> Y_cluster; 
    vector<double> Z_cluster; 
    vector<double> X_res_cluster;   // x residuals
    vector<double> Y_res_cluster;   // y residual 
    vector<double> R_res_cluster;   // R residuals
    vector<double> Phi_res_cluster; // phi residual 
    vector<double> R_cluster;       // radius
    vector<double> Phi_cluster;     // angle 
    vector<CLHEP::Hep3Vector> Points3d_cluster; // x,y measurement packed in
    double  edep[36] ;
    int nhitdev[36];
    CLHEP::Hep3Vector  MCPoint[36];
   
    X.clear();
    Y.clear();
    Z.clear();
    X_res.clear();
    Y_res.clear();
    R_res.clear();
    Phi_res.clear();
    R.clear();
    Phi.clear();
    Points3d.clear();
    
    X_straw.clear();
    Y_straw.clear();
    Z_straw.clear();
    X_res_straw.clear();
    Y_res_straw.clear();
    R_res_straw.clear();
    Phi_res_straw.clear();
    R_straw.clear();
    Phi_straw.clear();
    Points3d_straw.clear();

    X_cluster.clear();
    Y_cluster.clear();
    Z_cluster.clear();
    X_res_cluster.clear();
    Y_res_cluster.clear();
    R_res_cluster.clear();
    Phi_res_cluster.clear();
    R_cluster.clear();
    Phi_cluster.clear();
    Points3d_cluster.clear();
    multimap<int,pstraw> mpstraws;
    mpstraws.clear();

    //
    // Get the persistent data about pointers to StrawHits
    //
    edm::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_clmakerModuleLabel,"DPIStrawCluster",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

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
    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits sims( evt,
                               _g4ModuleLabel, 
                               _makerModuleLabel,
                               "tracker",
                               0.001,
                               5 );

    typedef SimParticlesWithHits::map_type map_type;
    for ( map_type::const_iterator i=sims.begin();
          i != sims.end(); ++i ){
      // All information about this SimParticle
      SimParticleInfo const& simInfo = i->second;
      // Information about StrawHits that belong on this SimParticle.
      vector<StrawHitMCInfo> const& infos = simInfo.strawHitInfos();
      if (simInfo.simParticle().generatorIndex()>=0)
	{
	  const ToyGenParticle genpar  =genParticles->at(simInfo.simParticle().generatorIndex());
	  
	  //cout<< genpar.generatorId()<<endl;
	  if (genpar.generatorId()== GenId::conversionGun)
	    {
	      cout << "SimParticle associated to conversion electron: "
		   << " Event: " << evt.id().event()
		   << " Track: " << i->first 
		   << " PdgId: " << simInfo.simParticle().pdgId() 
		   << " |p|: "   << simInfo.simParticle().startMomentum().vect().mag()
		   << " Hits: "  << infos.size()
		   << " CC:   "  << simInfo.simParticle().creationCode() 
		   << " GI:   "  << simInfo.simParticle().generatorIndex() 
		   << endl;
	      cout << "Polar: "<<simInfo.simParticle().startMomentum().vect().getTheta ()<<endl;
	      StepPointMC const& fstep =simInfo.firstStepPointMCinTracker();
	      StepPointMC const& lstep =simInfo.lastStepPointMCinTracker();
	      P_in_si= fstep.momentum();     // momentum as the track enters the tracker
	      Pt_inval_si =  P_in_si.rho();
	      P_inval_si  =  P_in_si.mag();
	      _Pt_in_si->Fill(Pt_inval_si);
	      _P_in_si ->Fill(P_inval_si);
	      _Pz_in_si->Fill(P_in_si.z());
	      P_out_si= lstep.momentum();   // momentum as the track leaves the tracker
	      Pt_outval_si =  P_out_si.rho();
	      P_outval_si  =  P_out_si.mag();
	      _Pt_out_si->Fill(Pt_outval_si);
	      _P_out_si ->Fill(P_outval_si);
	      _Pz_out_si->Fill(P_out_si.z());
	      // Loop over all StrawsHits to which this SimParticle contributed.
	      double _timetodist=149.8962;
	      _hNEleHits->Fill(infos.size());
	      // calculate the average hit position of track at a plane
	      
	      for (int i = 0; i < 36 ; i++) { 
		nhitdev[i] = 0 ;
		edep[i] = 0 ;
		MCPoint[i] = CLHEP::Hep3Vector(0.,0.,0.);
	      }
	      
	      for ( size_t j=0; j<infos.size(); ++j) // Loop over associated Hits
		{
		  StrawHitMCInfo const& info = infos.at(j);
		  StrawHit const& hit        = info.hit();
		  Straw const& str           = tracker.getStraw(hit.strawIndex());
		  sid = str.Id();
		  did = sid.getDeviceId();
		  std::vector<StepPointMC const *> const& steps = info.steps();
		  for ( size_t k=0; k<steps.size(); ++k){
		    StepPointMC const& step = *(steps[k]);
		    MCPoint[did] =  MCPoint[did]+step.position();
		    edep[did] =  edep[did]+step.totalEDep();
		    nhitdev[did]++;
		  }
		}
	      for (int i = 0; i < 36 ; i++) {
		if (nhitdev[i]>0)
		  {
		    MCPoint[i] = MCPoint[i]/float(nhitdev[i]);
		    //cout<<"Edep:  " <<  edep[i]<<"  MCPoint:  " << MCPoint[i]<<endl;
		  }      
	      }
	      for ( size_t j=0; j<infos.size(); ++j) // Loop over associated Hits
		{
		  StrawHitMCInfo const& info = infos.at(j);
		  StrawHit const& hit        = info.hit();
		  Straw const& str           = tracker.getStraw(hit.strawIndex());
		  const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
		  const CLHEP::Hep3Vector dirvec = str.getDirection();
		  double dt =hit.dt();
		  double disttomid =dt* _timetodist;
		  CLHEP::Hep3Vector hitpos = mpvec+disttomid*dirvec;
		  X_straw.push_back(hitpos.getX());
		  Y_straw.push_back(hitpos.getY());
		  Z_straw.push_back(hitpos.getZ());
		  Double_t Rhit = hitpos.rho();
		  R_straw.push_back(Rhit);
		  CLHEP::Hep3Vector smcpos = CLHEP::Hep3Vector( 0.0, 0.0, 0.0);
		  //		  cout << "steps.size():  " << steps.size()<<endl;
		  std::vector<StepPointMC const *> const& steps = info.steps();
		  for ( size_t k=0; k<steps.size(); ++k){
		    StepPointMC const& step = *(steps[k]);
		    smcpos= smcpos+step.position();
		  }
		  smcpos=smcpos/steps.size();
		  Points3d_straw.push_back(smcpos);
		  Double_t Rmc  = smcpos.rho();
		  R_res_straw.push_back(Rmc-Rhit);
		  _Rdiff_s-> Fill(Rmc-Rhit);
		  _Phidiff_s->Fill(smcpos.phi()-hitpos.phi());
		  X_res_straw.push_back(smcpos.getX()-hitpos.getX());
		  Y_res_straw.push_back(smcpos.getY()-hitpos.getY());
		  _Xdiff_s -> Fill(smcpos.getX()-hitpos.getX());
		  _Ydiff_s -> Fill(smcpos.getY()-hitpos.getY());		
		}
	      
	      FitCircle(X_straw, Y_straw);
	      _x0y0_s->Fill(x0,y0);
	      _R_rec_s->Fill(R_rec);
	      _chi2_s -> Fill(chi2) ;
	      
	      Double_t Bmagnet=10.;   // 10 KGauss magnetic field (hard wired should get from framework) 
	      Double_t Const=1.49898e-4;
	      cout << "Pt_inval:  "<< Pt_inval_si
		   << "  Radius:  "<< Radius
		   << "  curv:    "<< curv
		   << "  R_rec:  " << R_rec <<endl;
	      Pt =1000.*R_rec*0.1 * 2. * Bmagnet* Const;
	      _Pt_rec_s->Fill(Pt);
	      _Pt_diff_s->Fill(Pt-Pt_inval_si);
	      FitSinus(R_straw, Z_straw);
	      _Pz_rec_s->Fill(Pz);
	      _Pz_diff_s->Fill(Pz-P_in_si.z());
	      Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
	      _P_rec_s->Fill(Ptot);
	      _P_diff_s->Fill(Ptot-P_inval_si);
	    }
	}
    }
    Int_t totalHits=0;
    _hNClusters->Fill(mcptrHandle->size());
    CLHEP::Hep3Vector dvec;

    for ( size_t i=0; i< mcptrHandle->size(); ++i ) {
      double hlen=9999999.;
      DPIndexVector   const&    mcptr(hits_mcptr->at(i));
      CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
      _hNStraws->Fill(mcptr.size());
      totalHits = totalHits+mcptr.size();
      Double_t totalEnergy = 0.0;
      for( size_t j=0; j<mcptr.size(); j++ ) {
	DPIndex const& junkie = mcptr[j];
       	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
	Double_t Energy = strawhit.energyDep();
	//Double_t Time   = strawhit.time();
	//Double_t deltaT = strawhit.dt();
	totalEnergy=totalEnergy+Energy;
	str = tracker.getStraw(strawhit.strawIndex());
	sid = str.Id();
	lid = sid.getLayerId();
	did = sid.getDeviceId();
	secid = sid.getSectorId();
	const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
	const CLHEP::Hep3Vector dirvec = str.getDirection();
	dvec = CLHEP::Hep3Vector(dirvec.getX(),dirvec.getY(),dirvec.getZ());
	//pvec = pvec + Energy * mpvec;     // weight straw by energy deposition
	pvec = pvec + mpvec;
	if (str.getHalfLength()<hlen)
	  {
	    hlen=str.getHalfLength();
	  }
      }
      double a = 1./double(mcptr.size());
      pvec = pvec*a;
      //pvec = pvec/totalEnergy;             // mean weighted by energy deposition
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
	    multimap<int,pstraw>::iterator first11 = ppp1.first;
	    multimap<int,pstraw>::iterator first22 = ppp1.first;
	    multimap<int,pstraw>::iterator last1 = ppp1.second;
	    last1--;
	    multimap<int,pstraw>::iterator last2 = ppp1.second;
	    for ( multimap<int,pstraw>::iterator first1=first11;first1 != last1;first1++)
	      {
		first22=first1;
		first22++;
		for ( multimap<int,pstraw>::iterator first2=first22;first2 != last2;++first2)
		  {
		    pstraw junk  = (*first1).second;
		    pstraw pjunk = (*first2).second;
		    const CLHEP::Hep2Vector p0 = 
		      CLHEP::Hep2Vector(junk.mpx-junk.hl*junk.dirx,junk.mpy-junk.hl*junk.diry);
		    const CLHEP::Hep2Vector p1 = 
		      CLHEP::Hep2Vector(junk.mpx+junk.hl*junk.dirx,junk.mpy+junk.hl*junk.diry);
		    const CLHEP::Hep2Vector p2 =
		      CLHEP::Hep2Vector(pjunk.mpx-pjunk.hl*pjunk.dirx,pjunk.mpy-pjunk.hl*pjunk.diry); 
		    const CLHEP::Hep2Vector p3 = 
		      CLHEP::Hep2Vector(pjunk.mpx+pjunk.hl*pjunk.dirx,pjunk.mpy+pjunk.hl*pjunk.diry); 
		    LineSegment linesegment0(p0, p1);
		    LineSegment linesegment1(p2, p3);
		    CLHEP::Hep2Vector intersection;
		    switch(linesegment0.Intersect(linesegment1, intersection))
		      {
		      case PARALLEL:
			//std::cout << "The lines are parallel\n\n";
			break;
		      case COINCIDENT:
			//std::cout << "The lines are coincident\n\n";
			break;
		      case NOT_INTERSECTING:
			//std::cout << "The lines do not intersect\n\n";
			break;
		      case INTERSECTING:
			X.push_back(intersection.x());
			Y.push_back(intersection.y());
			Z.push_back(0.5*(junk.mpz+pjunk.mpz));
			CLHEP::Hep3Vector point3d  =
			  CLHEP::Hep3Vector(intersection.x(),intersection.y(),0.5*(junk.mpz+pjunk.mpz));
			Double_t Rhit=point3d.rho();
			Points3d.push_back(point3d);
			Double_t Rmc  = MCPoint[i].rho();
			R.push_back(Rhit);
			R_res.push_back(Rmc-Rhit);
			_Rdiff-> Fill(Rmc-Rhit);
			_Phidiff->Fill(MCPoint[i].phi()-point3d.phi());
			X_res.push_back(MCPoint[i].getX()-intersection.x());
			Y_res.push_back(MCPoint[i].getY()-intersection.y());
			_Xdiff -> Fill(MCPoint[i].getX()-intersection.x());
			_Ydiff -> Fill(MCPoint[i].getX()-intersection.y());
			nint ++;
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
	_x0y0->Fill(x0,y0);
	_R_rec->Fill(R_rec);
	_chi2 -> Fill(chi2) ;

	Double_t Bmagnet=10.;   // 10 KGauss magnetic field (hard wired should get from framework) 
	Double_t Const=1.49898e-4;
	cout << "Pt_inval:  "<< Pt_inval_si
	     << "  Radius:  "<< Radius
	     << "  curv:    "<< curv
             << "  zstep:   "<< zstep 
             << "  R_rec:  " << R_rec <<endl;
	Pt =1000.*R_rec*0.1 * 2. * Bmagnet* Const;
	_Pt_rec->Fill(Pt);
	_Pt_diff->Fill(Pt-Pt_inval_si);
	FitSinus(R, Z);
	_Pz_rec->Fill(Pz);
	_Pz_diff->Fill(Pz-P_in_si.z());
	Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
	_P_rec->Fill(Ptot);
	_P_diff->Fill(Ptot-P_inval_si);

      }

    int nclusters=0;
    double _timetodist=149.8962;
    for ( size_t i=0; i< mcptrHandle->size(); ++i )// Loop over Clusters
      {
	DPIndexVector   const&    mcptr(hits_mcptr->at(i));
	CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
	CLHEP::Hep3Vector clusterpos =  CLHEP::Hep3Vector(0.,0.,0.);
	//totalHits = totalHits+mcptr.size();
	Double_t totalEnergy = 0.0;
	for( size_t j=0; j<mcptr.size(); j++ ) // Loop over Straws in Cluster
	  {
	    DPIndex const& junkie = mcptr[j];
	    StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
	    Double_t Energy = strawhit.energyDep();
	    //Double_t Time   = strawhit.time();
	    Double_t deltaT = strawhit.dt();
	    StrawIndex si   = strawhit.strawIndex();
	    Straw str       = tracker.getStraw(si);	 
	    const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
	    const CLHEP::Hep3Vector dirvec = str.getDirection();
	    double disttomid = deltaT* _timetodist;   // convert delta T into delta x along the wire
	    CLHEP::Hep3Vector hitpos = mpvec+disttomid*dirvec;
	    clusterpos=clusterpos+hitpos;
	    totalEnergy=totalEnergy+Energy;
	  } // end loop over straws in cluster
	double a = 1./double(mcptr.size());
	clusterpos=clusterpos*a;
	X_cluster.push_back(clusterpos.getX());
	Y_cluster.push_back(clusterpos.getY());
	Z_cluster.push_back(clusterpos.getZ());
	Double_t Rhit =clusterpos.rho();
	R_cluster.push_back(Rhit);
	Points3d_cluster.push_back(clusterpos);
	Double_t Rmc  = MCPoint[i].rho();
	R_res_cluster.push_back(Rmc-Rhit);
	_Rdiff_c-> Fill(Rmc-Rhit);
	_Phidiff_c->Fill(MCPoint[i].phi()- clusterpos.phi());
	X_res_cluster.push_back(MCPoint[i].getX()-clusterpos.getX());
	Y_res_cluster.push_back(MCPoint[i].getY()-clusterpos.getY());
	_Xdiff_c -> Fill(MCPoint[i].getX()-clusterpos.getX());
	_Ydiff_c -> Fill(MCPoint[i].getY()-clusterpos.getY());
	nclusters++;

      } //  end Loop over Clusters      
    if (X_cluster.size()>2)
      {
	FitCircle(X_cluster, Y_cluster);
	_x0y0_c->Fill(x0,y0);
	_R_rec_c->Fill(R_rec);
	_chi2_c -> Fill(chi2) ;
	
	Double_t Bmagnet=10.;   // 10 KGauss magnetic field (hard wired should get from framework) 
	Double_t Const=1.49898e-4;
	cout << "Pt_inval:  "<< Pt_inval_si
	     << "  Radius:  "<< Radius
	     << "  curv:    "<< curv
             << "  zstep:   "<< zstep 
             << "  R_rec:  " << R_rec <<endl;
	Pt =1000.*R_rec*0.1 * 2. * Bmagnet* Const;
	_Pt_rec_c->Fill(Pt);
	_Pt_diff_c->Fill(Pt-Pt_inval_si);
	FitSinus(R_cluster, Z_cluster);
	_Pz_rec_c->Fill(Pz);
	_Pz_diff_c->Fill(Pz-P_in_si.z());
	Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
	_P_rec_c->Fill(Ptot);
	_P_diff_c->Fill(Ptot-P_inval_si);	
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
    cout << " Result: "<< result <<endl;
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
    Double_t  edm, errdef; 
    Int_t nvpar, nparx,istat;
    gmMinuit->mnstat(chi2,edm,errdef,nvpar,nparx,istat);

  }
}


using mu2e::ReadDPIStrawCluster;
DEFINE_FWK_MODULE(ReadDPIStrawCluster);
