//=============================================================================
//
// Plugin to study the persistent data about straw hits,
// looking at tau, beta, and tan theta for each hit.
//
// Starts from ReadDPIStrawCluster_plugin.cc, adding the quantities of
// interest to these angles, and gradually eliminating the rest.
//
// $Id: BetaTauPitch_module.cc,v 1.15 2013/10/21 21:01:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:23 $
//
// Original author: Mark Fischler modifying code by Hans Wenzel
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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TVirtualFitter.h"

// Mu2e includes.
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeneralUtilities/inc/LineSegmentPCA.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {
  enum PrintLevel { quiet  =-1,
                    normal = 0,
                    verbose= 1};
  Double_t Radius;
  Double_t curv;
  Double_t zstep;
  static bool magset(false);
  static Double_t Bmagnet;
  static Double_t Const(1.49898e-4);
   TGraph *gr2;
  TGraphErrors *error;
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = error->GetN();
   f = 0;
   Double_t *x = error->GetX();
   Double_t *y = error->GetY();
   Double_t *ex = error->GetEX();
//   Double_t *ey = error->GetEY();

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
//   Double_t *ex = error->GetEX();
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
      cout<< "Panel: " << sec  <<endl;
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
  class BetaTauPitch : public art::EDAnalyzer {
  public:
    explicit BetaTauPitch(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
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
      // mf study 2
      _EnergyDep_s(0),
      _EnergyDepX_s(0),
      // --- mf
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
      // mf study 1
      _beta_c(0),
      _tanTau_c(0),
      _tanTheta_c(0),
      // --- mf
      _ntup(0)
    {
    }
    virtual ~BetaTauPitch() { }

    virtual void beginJob();

    void analyze( art::Event const& e);
    void FitCircle(vector<double> X,vector<double> Y);
    void FitSinus( vector<double> R,vector<double> Z);
    void FitSinus2( vector<double> R,vector<double> Z);
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
    // mf study 2
    TH1F* _EnergyDep_s;
    TH1F* _EnergyDepX_s;
    // --- mf
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
    // mf study 1
    TH1F* _beta_c;
    TH1F* _tanTau_c;
    TH1F* _tanTheta_c;
    // --- mf
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

    // labels for positions in the ntuple
    enum ntpos{EVT,
               PGENX    ,PGENY ,PGENZ  ,
               PINX     ,PINY  ,PINZ   ,
               POUTX    ,POUTY ,POUTZ  ,
               NINT     ,RREC  ,PTREC ,PRECZ,
               NSTRAWS  ,RREC_S,PTREC_S,PREC_SZ,
               NCLUSTERS,RREC_C,PTREC_C,PREC_CZ,
               // mf study 1
               BETA_C,TANTAU_C,TANTHETA_C,
               // --- mf
               NT_COUNT};

  }; // end of BetaTauPich class definition


  void BetaTauPitch::beginJob(){
    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint<<endl;

    art::ServiceHandle<art::TFileService> tfs;
    _hNInter       = tfs->make<TH1F>( "hNInter",   "intersection ", 100  , 0., 100. );
    _hNClusters    = tfs->make<TH1F>( "hNClusters","Number of straw clusters", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",    "Number of straw Hits", 100, 0., 100. );
    _hNEleHits     = tfs->make<TH1F>( "hNEleHits", "Number of straw Hits/conversion electron", 100, 0., 100. );
    _hNStraws      = tfs->make<TH1F>( "hNStraws",  "Number of straws/cluster", 5  , 0., 5. );
    _Polar         = tfs->make<TH1F>( "Polar",     "polar angle at production", 120, -3., 3. );
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
    // mf study 2
    _EnergyDep_s     = tfs->make<TH1F>( "EnergyDep_s", "s Energy Deposited in Straw (KeV)", 100, 0.0, 10.0);
    _EnergyDepX_s    = tfs->make<TH1F>( "EnergyDepX_s","s Energy Deposited in Straw (KeV)", 100, 0.0, 10.0);
    // --- mf
    _R_rec_c         = tfs->make<TH1F>( "R_rec_c",     "c reconstructed track radius", 100, 0., 800. );
    _Pt_rec_c        = tfs->make<TH1F>( "Pt_rec_c",    "c reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec_c         = tfs->make<TH1F>( "P_rec_c",     "c reconstructed momentum", 100, 60., 140. );
    _Pz_rec_c        = tfs->make<TH1F>( "Pz_rec_c",    "c reconstructed longitudinal momentum", 100, 0., 120. );
    _Pt_diff_c       = tfs->make<TH1F>( "Pt_diff_c",   "c delta tansverse momentum Pt", 100, -20., 20. );
    _P_diff_c        = tfs->make<TH1F>( "P_diff_c",    "c delta momentum", 100, -20., 20.);
    _Pz_diff_c       = tfs->make<TH1F>( "Pz_diff_c",   "c delta longitudinal momentum", 100, -20., 20.);
    _x0y0_c          = tfs->make<TH2F>( "x0y0_c",      "c x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2_c          = tfs->make<TH1F>( "chi2_c",      "c chi2 of 2d circle fit", 200, 0., 200. );
    _Xdiff_c         = tfs->make<TH1F>( "Xdiff_c",     "c X(stepMC)  - X(stereo hit) mm ", 100, -20., 20.);
    _Ydiff_c         = tfs->make<TH1F>( "Ydiff_c",     "c Y(stepMC)  - Y(stereo hit)  mm ", 100, -20., 20.);
    _Rdiff_c         = tfs->make<TH1F>( "Rdiff_c",     "c R(stepMC)  - R(stereo hit) mm ", 100, -50., 50.);
    _Phidiff_c       = tfs->make<TH1F>( "Phidiff_c",   "c Phi(stepMC)- Phi(stereo hit) mm ", 100, -0.5, 0.5);
    // mf study 1
    _beta_c          = tfs->make<TH1F>( "beta_c",      "s straw incidence angle beta", 45, 0., 90. );
    _tanTau_c        = tfs->make<TH1F>( "tanTau_c",    "s panel attack angle tan tau", 50, 0., 2.0 );
    _tanTheta_c      = tfs->make<TH1F>( "tanTheta_c",  "s pitch tan theta", 50, 0., 2. );
    // --- mf
    _ntup          = tfs->make<TNtuple>( "ntup", "Pattern Recognition Ntuple",
                                                 "evt:Pgenx:Pgeny:Pgenz:Pinx:Piny:Pinz:Poutx:Pouty:Poutz:Nint:"
                                                 "Rrec:Ptrec:Precz:Nstraws:Rrec_s:Ptrec_s:Prec_sz:NClusters:"
                                                 "Rrec_c:Ptrec_c:Prec_cz:"
                                                 // mf study 1
                                                 "Beta_c:TanTau_c:TanTheta_c"
                                                // --- mf
                                        );
    assert(_ntup->GetNvar()==NT_COUNT);
    //
  }

  void BetaTauPitch::analyze(art::Event const& evt)
  {
    if ( _diagLevel > 2 ) cout << "BetaTauPitch: analyze() begin"<<endl;
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
        GeomHandle<BFieldConfig> bfconf;
        B= bfconf->getDSUniformValue();
        cout << " B-field (getDSUniformValue()):  " <<B<<endl;
        Bmagnet=B.getZ();
        magset=true;
      }
    //
    // Position of the center of the tracker in mu2e coordinates.
    //
    //CLHEP::Hep2Vector tt =CLHEP::Hep2Vector( 0.0, 10200.);
    //CLHEP::Hep3Vector point =CLHEP::Hep3Vector( -3904, 0.0, 10200.);
    //CLHEP::Hep3Vector bf = bfMgr->getBField(point);
    //cout << " B-field: (getBField(center of tracker)) " <<bf<<endl;
    static int ncalls(0);
    ++ncalls;
    StrawId nsid;
    Straw str;
    StrawId sid;
    LayerId lid;
    PanelId secid;
    // mf study 1
    int panel;
    // --- mf
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

    // mf study 1
    vector<CLHEP::Hep3Vector> momentum_cluster(36);
    vector<CLHEP::Hep3Vector> planeStrawDirections(36);
    vector<double> beta_cluster;            // angle against projection of wire
    vector<double> tanTau_cluster;          // attack angle of path to panel
    vector<double> tanTheta_cluster;        // helix pitch
    // --- mf
    //cout << "[[  1 ]]\n";
    vector<CLHEP::Hep3Vector> Points3d_cluster; // x,y measurement packed in
    //double  edep[36] ;
    int nhitplane[36];
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

    // mf study 1
    tanTheta_cluster.clear();     // helix pitch
    beta_cluster.clear();         // angle against projection of wire
    tanTau_cluster.clear();       // attack angle of path to panel
    // --- mf

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
    float nt[NT_COUNT];
    //cout << "[[  2 ]]\n";
    for (int i = 0;i<NT_COUNT;i++)
      {
        nt[i]=-9999.;
      }
    //
    // Get the persistent data about pointers to StrawHits
    //

    art::Handle<StrawClusterCollection> clustersHandle;
    evt.getByLabel(_clmakerModuleLabel,clustersHandle);
    StrawClusterCollection const& clusters = *clustersHandle;

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

    //cout << "[[  3 ]]\n";
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

            //cout<< genpar.generatorId()<<endl;
            if (genpar.generatorId()== GenId::conversionGun)
              {
                nt[PGENX]=genpar.momentum().getX();
                nt[PGENY]=genpar.momentum().getY() ;
                nt[PGENZ]=genpar.momentum().getZ() ;

                foundcele=true;
                /*
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
                */
                const double PolarAngleAtProduction = simInfo.simParticle().startMomentum().vect().getTheta ();
                _Polar->Fill(PolarAngleAtProduction);
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
                for (int iplane = 0; iplane < 36 ; iplane++) {
                  nhitplane[iplane] = 0 ;
                  //edep[iplane] = 0.0 ;
                  MCPoint[iplane] = CLHEP::Hep3Vector(0.,0.,0.);
                }
                for ( size_t associatedHit=0; associatedHit<infos.size(); ++associatedHit) // Loop over associated Hits
                  {
                    StrawHitMCInfo const& info = infos.at(associatedHit);
                    StrawHit const& hit        = info.hit();
                    Straw const& str           = tracker.getStraw(hit.strawIndex());
                    sid = str.id();
                    PlaneId did = sid.getPlaneId();
                    std::vector<StepPointMC const *> const& steps = info.steps();
                    // mf study 1
                    panel = sid.getPanel();
                    const CLHEP::Hep3Vector straw_direction = str.getDirection();
                    if ((straw_direction.z() > .000001) || (straw_direction.mag() > 1.000001)) {
                      cout << "????? unexpected straw direction: \n      " << straw_direction << "\n";
                    }
                    if ((panel < 0) || (panel > 5)) {
                      cout << "????? unexpected straw panel: \n      " << panel << "\n";
                    }
                    // --- mf
                    planeStrawDirections[did] = straw_direction;
                    //cout << "Plane " << did << " Straw Direction " << straw_direction << "\n";
                    double energyAH  = 0.;
                    double energyAHX = 0.;
                    for ( size_t ks=0; ks<steps.size(); ++ks){
                      StepPointMC const& step = *(steps[ks]);
                      if (step.momentum().mag()>5)
                        {
                          MCPoint[did] =  MCPoint[did]+step.position();
                          // mf study 2
                          // cout << "Energy Step =  " << step.totalEDep() << "\n";
                          energyAH  += step.totalEDep();
                          energyAHX += step.totalEDep();
                          // mf study 1
                          momentum_cluster[did] +=step.momentum();
                          // --- mf
                          nhitplane[did]++;
                        }
                      else // step momentum is less than 5
                        {
                          energyAH  += step.totalEDep();
                          //                      cout << "delta: " << step.momentum().mag()<<endl;
                        } // end of if/else for momentum > 5
                    } // end of loops over steps in this hit
                    // mf study 2
                    // cout << "Energy for associated hit " << associatedHit << " =  " << energyAH << "\n";
                    _EnergyDep_s->Fill(1000.0*energyAH);
                    _EnergyDepX_s->Fill(1000.0*energyAHX);
                    // --- mf
                  } // end of loop over associated hits

                // Plane quantity normalization loop
                for (int iplane = 0; iplane < 36 ; iplane++) {
                  if (nhitplane[iplane] <= 0) continue;
                  double a = 1.0/double(nhitplane[iplane]);
                  MCPoint[iplane] = MCPoint[iplane]*a;
                  // mf study 1
                  momentum_cluster[iplane] *= a;
                  CLHEP::Hep3Vector projectedMomentum = momentum_cluster[iplane];
                  projectedMomentum.setZ(0);
                  double beta = std::acos(projectedMomentum.dot(planeStrawDirections[iplane])/projectedMomentum.mag());
                  // double beta = 0;
                  _beta_c->Fill(beta*180.0/3.141592653589793);
                  //              if (nt[BETA_C] == -9999.) nt[BETA_C] = beta;
                  nt[BETA_C] = beta;
                  double tanTheta = momentum_cluster[iplane].perp()/momentum_cluster[iplane].z();
                  _tanTheta_c->Fill(tanTheta);
                  //              if (nt[TANTHETA_C] == -9999.) nt[TANTHETA_C] = tanTheta;
                  nt[TANTHETA_C] = tanTheta;
                  double tanTau = tanTheta*std::sin(beta);
                  _tanTau_c->Fill(tanTau);
                  //              if (nt[TANTAU_C] == -9999.) nt[TANTAU_C] = tanTau;
                  nt[TANTAU_C] = tanTau;
                  // --- mf
                } // End of Plane quantity normalization loop
                //cout << "[[  4 ]]\n";

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

                FitCircle(X_straw, Y_straw);
                _x0y0_s->Fill(x0,y0);
                _R_rec_s->Fill(R_rec);
                _chi2_s -> Fill(chi2) ;

                // Double_t Bmagnet=B.getZ();   // magnetic field
                // Double_t Const=1.49898e-4;
                //              cout << "Pt_inval:  "<< Pt_inval_si
                //     << "  Radius:  "<< Radius
                //     << "  curv:    "<< curv
                //     << "  R_rec:  " << R_rec <<endl;
                Pt =1000.*R_rec  * 2. * Bmagnet* Const;
                _Pt_rec_s->Fill(Pt);
                _Pt_diff_s->Fill(Pt-Pt_inval_si);
                FitSinus2(R_straw, Z_straw);
                _Pz_rec_s->Fill(Pz);
                _Pz_diff_s->Fill(Pz-P_in_si.z());
                Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                _P_rec_s->Fill(Ptot);
                _P_diff_s->Fill(Ptot-P_inval_si);
                //cout<<"Sinus:  " << Pz<<endl;
                nt[NSTRAWS]= X_straw.size();
                nt[RREC_S] = R_rec;
                nt[PTREC_S]= Pt;
                nt[PREC_SZ]= Pz;
                //FitSinus2(Z_straw, R_straw);
                //cout<<"Sinus2: " << Pz<<endl;
              }   // end code done if the simparticle is conversionGun
          }       // end if on generatorINdex >= 0
      }           // end loop over simparticles
    if (!foundcele) return;       // no conversion electron found
    Int_t totalHits=0;
    _hNClusters->Fill(clusters.size());
    CLHEP::Hep3Vector dvec;
    //cout << "[[  5 ]]\n";

    for ( size_t i=0; i<clusters.size(); ++i ) { // Apparently a loop over clusters
      PlaneId did = -1;
      double hlen=9999999.;
      StrawCluster      const& cluster(clusters.at(i));
      StrawHitPtrVector const& strawHits(cluster.strawHits());
      CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
      _hNStraws->Fill(strawHits.size());
      totalHits = totalHits+strawHits.size();
      Double_t totalEnergy = 0.0;
      for( size_t j=0; j<strawHits.size(); j++ ) {  // Loop over straws in the cluster
        StrawHit const& strawhit = *strawHits[j];
        Double_t Energy = strawhit.energyDep();
        //Double_t Time   = strawhit.time();
        //Double_t deltaT = strawhit.dt();
        totalEnergy=totalEnergy+Energy;
        str = tracker.getStraw(strawhit.strawIndex());
        sid = str.id();
        lid = sid.getLayerId();
        did = sid.getPlaneId();
        secid = sid.getPanelId();
        const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
        const CLHEP::Hep3Vector dirvec = str.getDirection();
        dvec = CLHEP::Hep3Vector(dirvec.getX(),dirvec.getY(),dirvec.getZ());
        //pvec = pvec + Energy * mpvec;     // weight straw by energy deposition
        pvec = pvec + mpvec;
        if (str.getHalfLength()<hlen)
          {
            hlen=str.getHalfLength();
          }
      } // end of loop over straws in the cluster
      double a = 1./double(strawHits.size());
      pvec = pvec*a;
      //pvec = pvec/totalEnergy;             // mean weighted by energy deposition
      pstraw pstr;
      pstr.lay=lid.getLayer();
      pstr.did=did;
      pstr.sec=secid.getPanel();
      pstr.hl=hlen;
      pstr.mpx=pvec.getX();
      pstr.mpy=pvec.getY();
      pstr.mpz=pvec.getZ();
      pstr.dirx=dvec.getX();
      pstr.diry=dvec.getY();
      pstr.dirz=dvec.getZ();
      mpstraws.insert(pair<int,pstraw>(did,pstr));
    }
    //cout << "[[  6 ]]\n";
    _hNHits->Fill(totalHits);
    //cout << "[[  7 ]]\n";
    //cout << " size of pseudo straw map: " <<mpstraws.size()<<endl;

    // Loop over PAIRS of pseudostraws in the same plane, but only if they intersect
    // This will create doublets, and fill:
    //
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
      }   ///endloop over all planes

    // Fit to circle if there are at least three points
    //cout << "[[  8 ]]\n";
    _hNInter->Fill(X.size());
    //cout << "[[  9 ]]\n";
    if (X.size()>2)
      {
        FitCircle(X, Y);
        _x0y0->Fill(x0,y0);
        _R_rec->Fill(R_rec);
        _chi2 -> Fill(chi2) ;

        //      Double_t Bmagnet=10.;   // 10 KGauss magnetic field (hard wired should get from framework)
        // Double_t Const=1.49898e-4;
        //      cout << "Pt_inval:  "<< Pt_inval_si
        //     << "  Radius:  "<< Radius
        //    << "  curv:    "<< curv
        //     << "  zstep:   "<< zstep
        //     << "  R_rec:  " << R_rec <<endl;
        Pt =1000.*R_rec * 2. * Bmagnet* Const;
        _Pt_rec->Fill(Pt);
        _Pt_diff->Fill(Pt-Pt_inval_si);
        FitSinus2(R, Z);
        _Pz_rec->Fill(Pz);
        _Pz_diff->Fill(Pz-P_in_si.z());
        Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
        _P_rec->Fill(Ptot);
        _P_diff->Fill(Ptot-P_inval_si);
        nt[NINT]= X.size();
        nt[RREC] = R_rec;
        nt[PTREC]= Pt;
        nt[PRECZ]= Pz;

      }

    //cout << "[[ 10 ]]\n";
    int nclusters=0;
    double _timetodist=149.8962;
    for ( size_t i=0; i< clusters.size(); ++i )// Loop over Clusters
      {
        StrawCluster      const& cluster(clusters.at(i));
        StrawHitPtrVector const& strawHits(cluster.strawHits());
        CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
        CLHEP::Hep3Vector clusterpos =  CLHEP::Hep3Vector(0.,0.,0.);
        //totalHits = totalHits+strawHits.size();
        double totalEnergy = 0.0;
        for( size_t j=0; j<strawHits.size(); j++ ) // Loop over Straws in Cluster
          {
            StrawHit const& strawhit = *strawHits[j];
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
        double a = 1./double(strawHits.size());
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
    //cout << "[[ 11 ]]\n";
    if (X_cluster.size()>2)
      {
        FitCircle(X_cluster, Y_cluster);
        _x0y0_c->Fill(x0,y0);
        _R_rec_c->Fill(R_rec);
        _chi2_c -> Fill(chi2) ;

        //      double Bmagnet=10.;   // 10 KGauss magnetic field (hard wired should get from framework)
        // double Const=1.49898e-4;
        //cout << "Pt_inval:  "<< Pt_inval_si
        //     << "  Radius:  "<< Radius
        //     << "  curv:    "<< curv
        //     << "  zstep:   "<< zstep
        //     << "  R_rec:  " << R_rec <<endl;
        Pt =1000.*R_rec * 2. * Bmagnet* Const;
        _Pt_rec_c->Fill(Pt);
        _Pt_diff_c->Fill(Pt-Pt_inval_si);
        FitSinus2(R_cluster, Z_cluster);
        _Pz_rec_c->Fill(Pz);
        _Pz_diff_c->Fill(Pz-P_in_si.z());
        double Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
        _P_rec_c->Fill(Ptot);
        _P_diff_c->Fill(Ptot-P_inval_si);
        nt[NCLUSTERS]= X_cluster.size();
        nt[RREC_C]   = R_rec;
        nt[PTREC_C]  = Pt;
        nt[PREC_CZ]  = Pz;
      }
    //cout << "[[ 12 ]]\n";
    _ntup->Fill(nt);
    //cout << "[[ 13 ]]\n";
  } // end of ::analyze.

  void BetaTauPitch::FitSinus(    vector<double> R,     vector<double> Z)
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
  void BetaTauPitch::FitCircle(    vector<double> X,vector<double> Y)
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
        // cout <<"-----------Circle fit didn't converge---------------------------" <<endl;
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
  void BetaTauPitch::FitSinus2(    vector<double> X,vector<double> Y)
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
    //cout << " Result: "<< result <<endl;
    bool converged = gmMinuit2->fCstatu.Contains("CONVERGED");
    if (!converged)
      {
        // cout <<"-----------Sin fit didn't converge---------------------------" <<endl;
        return;
      }
    for (int i = 0;i<dim;i++) {
     gmMinuit2->GetParameter(i,sfpar[i],errsfpar[i]);
    }
    Double_t p2 = sfpar[2];
    Pz = 10./(33.36*p2);
    Double_t  edm, errdef;
    Int_t nvpar, nparx,istat;
    gmMinuit2->mnstat(chi2,edm,errdef,nvpar,nparx,istat);

  }
}


using mu2e::BetaTauPitch;
DEFINE_ART_MODULE(BetaTauPitch);
