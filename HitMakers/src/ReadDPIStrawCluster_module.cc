//=============================================================================
//
// Plugin to test that I can read back the persistent data about straw hits.
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
// Uses three methods to reconstruct the track helix of a conversion electron:
// * use single straws associated to Simparticle., use deltat to form 3D
//   reconstructed points
// * use intersection of Pseudo straws and use stereo information to form
//   3D reconstructed points
// * form pseudo straws from clusters in a panel. Use the combined information
//   to form 3D reconstructed points
// For all three cases estimate Pt,Pz of the conversion electron by performing
// a simple circle/sinus fit.
//
// $Id: ReadDPIStrawCluster_module.cc,v 1.6 2011/05/19 23:51:50 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/19 23:51:50 $
//
// Original author: Hans Wenzel
//
//=============================================================================
// C++ includes
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "CLHEP/Vector/TwoVector.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"

// Mu2e includes.
#include "BFieldGeom/inc/BFieldManager.hh"
#include "CDFTrajectory/inc/Helix.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "Mu2eUtilities/inc/LineSegmentPCA.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/GenId.hh"
#include "ToyDP/inc/GenParticle.hh"
#include "ToyDP/inc/GenParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/StatusG4.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "TrackerGeom/inc/Tracker.hh"
using namespace std;

namespace mu2e {;
  enum PrintLevel { quiet  =-1,
                    normal = 0,
                    verbose= 1};
  Double_t Radius;
  Double_t curv;
  Double_t zstep;
  static bool magset(false);
  static Double_t Bmagnet;
  static Double_t Const(1.49898e-4);
  TGraphErrors *error;
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = error->GetN();
   f = 0;
   Double_t *x = error->GetX();
   Double_t *y = error->GetY();
   Double_t *ex = error->GetEX();
   //Double_t *ey = error->GetEY();

   for (Int_t i=0;i<np;i++) {
      Double_t u = (x[i] - par[0]);
      Double_t v = (y[i] - par[1]);
      Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
      //f += (dr*dr)/(ex[i]*ex[i]+ey[i]*ey[i]);
      f += (dr*dr)/(0.25*ex[i]*ex[i]);
   }
}
void myfcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = error->GetN();
   f = 0;
   Double_t *x = error->GetX();
   Double_t *y = error->GetY();
   //Double_t *ex = error->GetEX();
   Double_t *ey = error->GetEY();
   for (Int_t i=0;i<np;i++) {
     double_t dr = y[i] - (par[0]+par[1]*TMath::Sin(par[2]*x[i]-par[3]));
     f += (dr*dr)/(ey[i]*ey[i]);
   }
}
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
  class ReadDPIStrawCluster : public art::EDAnalyzer {
  public:
    explicit ReadDPIStrawCluster(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _clmakerModuleLabel(pset.get<std::string>("clmakerModuleLabel")),
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
      _ntup(0)
    {
    }
    virtual ~ReadDPIStrawCluster() { }

    virtual void beginJob();

    void analyze( art::Event const& e);
    bool FitCircle(vector<double> X,vector<double> Y);
    bool FitSinus2( vector<double> R,vector<double> Z);
  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    // Module label of the geerator module.
    std::string _generatorModuleLabel;
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
    TNtuple* _ntup;
    //
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


  void ReadDPIStrawCluster::beginJob(){
    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint<<endl;

    art::ServiceHandle<art::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",    "Number of straw Hits", 100, 0., 100. );
    _hNEleHits     = tfs->make<TH1F>( "hNEleHits", "Number of straw Hits/conversion electron", 100, 0., 100. );
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
    _Rdiff_s         = tfs->make<TH1F>( "Rdiff_s",     "s R(stepMC)  - R(straw hit) mm ", 100, -15., 15.);
    _Phidiff_s       = tfs->make<TH1F>( "Phidiff_s",   "s Phi(stepMC)- Phi(straw hit) mm ", 100, -0.5, 0.5);
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
    _Phidiff_c       = tfs->make<TH1F>( "Phidiff_c",   "s Phi(stepMC)- Phi(stereo hit) mm ", 100, -0.5, 0.5);

    _ntup          = tfs->make<TNtuple>( "ntup", "Pattern Recognition Ntuple",
                      "evt:Pgenx:Pgeny:Pgenz:Pinx:Piny:Pinz:Poutx:Pouty:Poutz:Nint:Rrec:Ptrec:Precz:Nstraws:Rrec_s:Ptrec_s:Prec_sz:NClusters:Rrec_c:Ptrec_c:Prec_cz");
    //
  }

  void ReadDPIStrawCluster::analyze(art::Event const& evt)
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
    if (!magset)
      {
        GeomHandle<BFieldManager> bfMgr;
        B= bfMgr->getDSUniformValue();
        cout << " B-field (getDSUniformValue()):  " <<B<<endl;
        Bmagnet=B.getZ();
        magset=true;

        const HepGeom::Point3D<double> position = HepGeom::Point3D<double>(0.,0.,0.);
        const HepGeom::Vector3D<double> direction=HepGeom::Vector3D<double>(1.,1.,1.);
        Helix helix = Helix(direction,position,-1,Bmagnet);
        cout << "Curvature:           " << helix.getCurvature()<<endl;
        cout << "Helicity:            " << helix.getHelicity() << endl;
        cout << "cotangent of theta:  " << helix.getCotTheta() << endl;
        cout << "phi0:                " << helix.getPhi0() << endl;
        cout << "d0:                  " << helix.getD0() << endl;
        cout << "Z0:                  " << helix.getZ0() << endl;
        cout << "Position after s=3.: " << helix.getPosition(3.) <<endl;
        cout << "Direction after s=3: " << helix.getDirection(3.) <<endl;
      }

    static int ncalls(0);
    ++ncalls;
    //    static Double_t Const(1.49898e-4);
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
    float nt[22];
    enum ntpos{EVT,
               PGENX    ,PGENY ,PGENZ  ,
               PINX     ,PINY  ,PINZ   ,
               POUTX    ,POUTY ,POUTZ  ,
               NINT     ,RREC  ,PTREC ,PRECZ,
               NSTRAWS  ,RREC_S,PTREC_S,PREC_SZ,
               NCLUSTERS,RREC_C,PTREC_C,PREC_CZ};
    for (int i = 0;i<22;i++)
      {
        nt[i]=-9999.;
      }
    //
    // Get the persistent data about pointers to StrawHits
    //
    art::Handle<DPIndexVectorCollection> mcptrHandle;
    evt.getByLabel(_clmakerModuleLabel,"DPIStrawCluster",mcptrHandle);
    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    evt.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);

    //Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits sims( evt,
                               _g4ModuleLabel,
                               _makerModuleLabel,
                               "tracker",
                               0.001,
                               5 );
    if (sims.size()<1) return;  // no sim particles found
    nt[EVT]  = evt.id().event();
    bool foundcele=false;
    typedef SimParticlesWithHits::map_type map_type;
    for ( map_type::const_iterator i=sims.begin();
          i != sims.end(); ++i )                      // loop over simparticles
      {
        // All information about this SimParticle
        SimParticleInfo const& simInfo = i->second;
        // Information about StrawHits that belong on this SimParticle.
        vector<StrawHitMCInfo> const& infos = simInfo.strawHitInfos();
        if (simInfo.simParticle().generatorIndex()>=0)
          {
            const GenParticle genpar  =genParticles->at(simInfo.simParticle().generatorIndex());
            if (genpar.generatorId()== GenId::conversionGun)
              {
                nt[PGENX]=genpar.momentum().getX();
                nt[PGENY]=genpar.momentum().getY() ;
                nt[PGENZ]=genpar.momentum().getZ() ;

                foundcele=true;
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
                nt[PINX] = P_in_si.getX();
                nt[PINY] = P_in_si.getY();
                nt[PINZ] = P_in_si.getZ();
                nt[POUTX]= P_out_si.getX();
                nt[POUTY]= P_out_si.getY();
                nt[POUTZ]= P_out_si.getZ() ;
                // Loop over all StrawsHits to which this SimParticle contributed.
                double _timetodist=149.8962;
                _hNEleHits->Fill(infos.size());
                // calculate the average hit position of track at a plane

                for (int idev = 0; idev < 36 ; idev++) {
                  nhitdev[idev] = 0 ;
                  edep[idev] = 0.0 ;
                  MCPoint[idev] = CLHEP::Hep3Vector(0.,0.,0.);
                }

                for ( size_t jdev=0; jdev<infos.size(); ++jdev) // Loop over associated Hits
                  {
                    StrawHitMCInfo const& info = infos.at(jdev);
                    StrawHit const& hit        = info.hit();
                    Straw const& str           = tracker.getStraw(hit.strawIndex());
                    sid = str.Id();
                    did = sid.getDeviceId();
                    std::vector<StepPointMC const *> const& steps = info.steps();
                    for ( size_t ks=0; ks<steps.size(); ++ks){
                      StepPointMC const& step = *(steps[ks]);
                      if (step.momentum().mag()>5)
                        {
                          MCPoint[did] =  MCPoint[did]+step.position();

                          edep[did] =  edep[did]+step.totalEDep();
                          nhitdev[did]++;
                        }

                    }
                  }                                             // end loop over associated Hits
                for (int idev = 0; idev < 36 ; idev++) {
                  if (nhitdev[idev]>0)
                    {
                      double a = 1.0/double(nhitdev[idev]);
                      MCPoint[idev] = MCPoint[idev]*a;
                    }
                }
                for ( size_t jhit=0; jhit<infos.size(); ++jhit) // Loop over associated Hits
                  {
                    StrawHitMCInfo const& info = infos.at(jhit);
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
                  }                        // end loop over hits
                nt[NSTRAWS]= X_straw.size();
                if ( X_straw.size()>4)
                  {
                    if (FitCircle(X_straw, Y_straw))
                      {
                        _x0y0_s->Fill(x0,y0);
                        _R_rec_s->Fill(R_rec);
                        _chi2_s -> Fill(chi2) ;
                        Pt =1000.*R_rec  * 2. * Bmagnet* Const;
                        _Pt_rec_s->Fill(Pt);
                        _Pt_diff_s->Fill(Pt-Pt_inval_si);
                        nt[RREC_S] = R_rec;
                        nt[PTREC_S]= Pt;
                        if (FitSinus2(Z_straw, R_straw))
                          {
                            _Pz_rec_s->Fill(Pz);
                            _Pz_diff_s->Fill(Pz-P_in_si.z());
                            Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                            _P_rec_s->Fill(Ptot);
                            _P_diff_s->Fill(Ptot-P_inval_si);
                            nt[PREC_SZ]= Pz;
                          }
                      }
                  }
              }
        }
      }           // end loop over simparticles
    if (!foundcele) return;       // no conversion electron found
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
                    LineSegmentPCA linesegment0(p0, p1);
                    LineSegmentPCA linesegment1(p2, p3);
                    CLHEP::Hep2Vector intersection;
                    switch(linesegment0.Intersect(linesegment1, intersection))
                      {
                      case LineSegmentPCA::PARALLEL:
                        //std::cout << "The lines are parallel\n\n";
                        break;
                      case LineSegmentPCA::COINCIDENT:
                        //std::cout << "The lines are coincident\n\n";
                        break;
                      case LineSegmentPCA::NOT_INTERSECTING:
                        //std::cout << "The lines do not intersect\n\n";
                        break;
                      case LineSegmentPCA::INTERSECTING:
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
    nt[NINT] = X.size();
    if (X.size()>4)
      {
        if (FitCircle(X, Y))
          {
            _x0y0->Fill(x0,y0);
            _R_rec->Fill(R_rec);
            _chi2 -> Fill(chi2) ;
            Pt =1000.*R_rec * 2. * Bmagnet* Const;
            _Pt_rec->Fill(Pt);
            _Pt_diff->Fill(Pt-Pt_inval_si);
            nt[RREC] = R_rec;
            nt[PTREC]= Pt;
            if (FitSinus2(Z, R))
              {
                _Pz_rec->Fill(Pz);
                _Pz_diff->Fill(Pz-P_in_si.z());
                Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                _P_rec->Fill(Ptot);
                _P_diff->Fill(Ptot-P_inval_si);
                nt[PRECZ]= Pz;
              }
          }
      }

    int nclusters=0;
    double _timetodist=149.8962;
    for ( size_t i=0; i< mcptrHandle->size(); ++i )// Loop over Clusters
      {
        DPIndexVector   const&    mcptr(hits_mcptr->at(i));
        CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
        CLHEP::Hep3Vector clusterpos =  CLHEP::Hep3Vector(0.,0.,0.);
        //totalHits = totalHits+mcptr.size();
        double totalEnergy = 0.0;
        for( size_t j=0; j<mcptr.size(); j++ ) // Loop over Straws in Cluster
          {
            DPIndex const& junkie = mcptr[j];
            StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(evt,junkie);
            double Energy = strawhit.energyDep();
            //double Time   = strawhit.time();
            double deltaT = strawhit.dt();
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
        double Rhit =clusterpos.rho();
        R_cluster.push_back(Rhit);
        Points3d_cluster.push_back(clusterpos);
        double Rmc  = MCPoint[i].rho();
        R_res_cluster.push_back(Rmc-Rhit);
        _Rdiff_c-> Fill(Rmc-Rhit);
        _Phidiff_c->Fill(MCPoint[i].phi()- clusterpos.phi());
        X_res_cluster.push_back(MCPoint[i].getX()-clusterpos.getX());
        Y_res_cluster.push_back(MCPoint[i].getY()-clusterpos.getY());
        _Xdiff_c -> Fill(MCPoint[i].getX()-clusterpos.getX());
        _Ydiff_c -> Fill(MCPoint[i].getY()-clusterpos.getY());
        nclusters++;

      } //  end Loop over Clusters
    nt[NCLUSTERS]= X_cluster.size();
    if (X_cluster.size()>4)
      {
        if (FitCircle(X_cluster, Y_cluster))
          {
            _x0y0_c->Fill(x0,y0);
            _R_rec_c->Fill(R_rec);
            _chi2_c -> Fill(chi2) ;
            Pt =1000.*R_rec * 2. * Bmagnet* Const;
            _Pt_rec_c->Fill(Pt);
            _Pt_diff_c->Fill(Pt-Pt_inval_si);
            nt[RREC_C]   = R_rec;
            nt[PTREC_C]  = Pt;
            if (FitSinus2(Z_cluster, R_cluster))
              {
                _Pz_rec_c->Fill(Pz);
                _Pz_diff_c->Fill(Pz-P_in_si.z());
                double Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                _P_rec_c->Fill(Ptot);
                _P_diff_c->Fill(Ptot-P_inval_si);
                nt[PREC_CZ]  = Pz;

              }
          }
      }
    _ntup->Fill(nt);
  } // end of ::analyze.

  bool ReadDPIStrawCluster::FitCircle(    vector<double> X,vector<double> Y)
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
    //int result=gmMinuit->Migrad();
    //cout << " Result: "<< result <<endl;
    bool converged = gmMinuit->fCstatu.Contains("CONVERGED");
    if (!converged)
      {
        cout <<"-----------Circle fit didn't converge---------------------------" <<endl;
        return converged;
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
    return converged;
  }
  bool  ReadDPIStrawCluster::FitSinus2( vector<double> X,vector<double> Y)
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
    TMinuit *gmMinuit2 = new TMinuit(4);
    gmMinuit2->SetPrintLevel(quiet);
    gmMinuit2->SetFCN(myfcn2);
    const int dim(4);
    const char par_name[dim][20]={"offset","radius","frequency","phase"};
    static Double_t step[dim] = {0.001,0.001,0.001,0.001};
    Double_t offset =  TMath::Sqrt(x0*x0+y0*y0);
    Double_t radius= R_rec;
    Double_t sfpar[dim]={offset,radius,0.005,0.1};
    Double_t errsfpar[dim]={0.0,0.0,0.0,0.0};
    int ierflg = 0;
    for (int ii = 0; ii<dim; ii++) {
      gmMinuit2->mnparm(ii,par_name[ii],sfpar[ii], step[ii], 0,0,ierflg);
    }
    gmMinuit2->FixParameter(0);
    gmMinuit2->FixParameter(1);
    //int result=gmMinuit2->Migrad();
    bool converged = gmMinuit2->fCstatu.Contains("CONVERGED");
    if (!converged)
      {
        cout <<"-----------Sin fit didn't converge---------------------------" <<endl;
        return converged;
      }
    for (int i = 0;i<dim;i++) {
     gmMinuit2->GetParameter(i,sfpar[i],errsfpar[i]);
    }
    Double_t p2 = sfpar[2];
    Pz = 10./(33.36*p2);
    Double_t  edm, errdef;
    Int_t nvpar, nparx,istat;
    gmMinuit2->mnstat(chi2,edm,errdef,nvpar,nparx,istat);
    return converged;
  }
}


using mu2e::ReadDPIStrawCluster;
DEFINE_ART_MODULE(ReadDPIStrawCluster);
