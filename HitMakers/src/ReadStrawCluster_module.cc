//
// Plugin to test that I can read back the persistent data about straw hits.
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadStrawCluster_module.cc,v 1.23 2013/10/21 21:01:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:23 $
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
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>
#include <utility>
//#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TArc.h"
#include "TMath.h"
#include "TMinuit.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GeneralUtilities/inc/LineSegmentPCA.hh"
#include "HitMakers/inc/StrawClusterUtilities.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
using namespace std;

namespace mu2e {
  enum PrintLevel { quiet  =-1,
                    normal = 0,
                    verbose= 1};
  enum ntpos{EVT,
             PGENX    ,PGENY ,PGENZ  ,
             PINX     ,PINY  ,PINZ   ,
             POUTX    ,POUTY ,POUTZ  ,
             NINT     ,RREC  ,PTREC ,PRECZ,
             NSTRAWS  ,RREC_S,PTREC_S,PREC_SZ,
             RREC_H   ,PTREC_h,
             NCLUSTERS,RREC_C,PTREC_C,PREC_CZ,NDIM};
  float nt[NDIM];
  TGraph *gr;
  static bool magset(false);
  static Double_t Bmagnet; 
  static Double_t Const(1.49898e-4);
  TGraphErrors *error;
  Double_t R_rec,x0,y0,chi2;
  Double_t Pt,Pz;
  Double_t Pt_inval_si;
  Double_t P_inval_si;
  Double_t Pt_outval_si;
  Double_t P_outval_si;
  StrawId sid;
  LayerId lid;
  DeviceId did;
  SectorId secid;
  CLHEP::Hep3Vector  X_in;  
  CLHEP::Hep3Vector  P_in_si;
  CLHEP::Hep3Vector  P_out_si;
  double  edep[36] ;
  int nhitdev[36];
  CLHEP::Hep3Vector  MCPoint[36];
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


  struct mpoint {
    double x;
    double y;
    double z;
    double errx;
    double erry;
  };
  vector<mpoint> mpoints;
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
   //Double_t *ex = error->GetEX();
   Double_t *ey = error->GetEY();
   for (Int_t i=0;i<np;i++) {
     double_t dr = y[i] - (par[0]+par[1]*TMath::Sin(par[2]*x[i]-par[3]));
     f += (dr*dr)/(ey[i]*ey[i]);
   }
}


  //--------------------------------------------------------------------
  //
  //
  class ReadStrawCluster : public art::EDAnalyzer {
  public:
    explicit ReadStrawCluster(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel","g4run")),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _clmakerModuleLabel(pset.get<std::string>("clmakerModuleLabel")),
      _hNInter(0),
      _hNClusters(0),
      _hNStraws(0),
      _R_rec(0),
      _x0y0(0),
      _chi2(0),
      _ntup(0)
    {
    }
    virtual ~ReadStrawCluster() { }

    virtual void beginJob();
    bool FitCircle(vector<double> X,vector<double> Y);
    bool FitSinus2( vector<double> R,vector<double> Z);
    void analyze( art::Event const& e);
  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Module label of the generator module.
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
    TH1F* _hNStraws;
    TH1F* _R_rec;
    TH2F* _x0y0; 
    TH1F* _chi2;
    TH1F* _Pt_rec;
    TH1F* _P_rec;
    TH1F* _Pz_rec;
    TNtuple* _ntup;
    CLHEP::Hep3Vector B;
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
    _R_rec         = tfs->make<TH1F>( "R_rec",  "reconstructed track radius", 100, 200., 400. );
    _x0y0          = tfs->make<TH2F>( "x0y0","x0 of circle vs y0 of circle ", 500,-650.,650.,500,-650.,650.);
    _chi2          = tfs->make<TH1F>( "chi2",      "chi2 of 2d circle fit", 200, 0., 200. );
    _Pt_rec        = tfs->make<TH1F>( "Pt_rec",    "reconstructed tansverse momentum Pt", 100, 0., 160. );
    _P_rec         = tfs->make<TH1F>( "P_rec",     "reconstructed momentum", 100, 60., 140. );
    _Pz_rec        = tfs->make<TH1F>( "Pz_rec",    "reconstructed longitudinal momentum", 100, 0., 120. );
    _ntup          = tfs->make<TNtuple>( "ntup", "Pattern Recognition Ntuple", 
                      "evt:Pgenx:Pgeny:Pgenz:Pinx:Piny:Pinz:Poutx:Pouty:Poutz:Nint:Rrec:Ptrec:Precz:Nstraws:Rrec_s:Ptrec_s:Prec_sz:Rrec_h:Ptrec_h:NClusters:Rrec_c:Ptrec_c:Prec_cz");
  }

  void ReadStrawCluster::analyze(art::Event const& evt) {
    
    if (!magset)
      {
        GeomHandle<BFieldConfig> bfconf;
        B= bfconf->getDSUniformValue(); 
        cout << " B-field (getDSUniformValue()):  " <<B<<endl;
        Bmagnet=B.getZ();
        magset=true;       
      }
    cout << "ReadStrawCluster: analyze() begin"<<endl;
    for (int i = 0;i<NDIM;i++)
      {
        nt[i]=-9999.;
      }

    if ( _diagLevel > 2 ) cout << "ReadStrawCluster: analyze() begin"<<endl;
    static int ncalls(0);
    vector<double> X;       // x of cluster intersections
    vector<double> Y;       // y of cluster intersections
    vector<double> Z;       // z of cluster intersections
    vector<double> R;       // z of cluster intersections
    StrawClusterUtilities sutils=StrawClusterUtilities() ;        
    multimap<int,StrawCluster> clusterbydid;
    multimap<int,StrawCluster> clusterbystation;
    clusterbydid.clear();
    clusterbystation.clear();

    const Tracker& tracker = getTrackerOrThrow();
  
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

   // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits sims( evt,
                               _g4ModuleLabel, 
                               _makerModuleLabel,
                               "tracker",
                              0.001,
                               5 );

    if (sims.size()<1) return;  // no sim particles found
    cout << "Nr of clusters:   " << clusters.size()<<endl;
    cout << "Nr of sim particles:   " << sims.size()<<endl;




    //bool foundcele=false;
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
            cout <<"generator index:  "<<genpar.generatorId()<<endl;
            if (genpar.generatorId()== GenId::conversionGun)
              {
                cout <<"conversion:  "<<genpar.generatorId()<<endl;

                //nt[PGENX]=genpar.momentum().getX();
                //nt[PGENY]=genpar.momentum().getY() ;
                 //nt[PGENZ]=genpar.momentum().getZ() ;

                //foundcele=true;
                StepPointMC const& fstep =simInfo.firstStepPointMCinTracker();
                StepPointMC const& lstep =simInfo.lastStepPointMCinTracker();



                cout<< 0.001*fstep.momentum()<<endl;
                cout << "Pt: "<<0.001*TMath::Sqrt(fstep.momentum().x()*fstep.momentum().x()+
                                    fstep.momentum().y()*+fstep.momentum().y())<<endl;
                P_in_si =fstep.momentum();// momentum as the track enters the tracker
                cout<<" P_in_si:  "<<P_in_si<<endl;

                cout << "Pt(rh0): "<<0.001* P_in_si.rho()<<endl;
                cout<< "radius:  =" <<0.001* P_in_si.rho()/(1.49898e-3*Bmagnet*2.)<<endl;
                cout << "P: "<<0.001*TMath::Sqrt(fstep.momentum().x()*fstep.momentum().x()+
                                    fstep.momentum().y()*+fstep.momentum().y()+
                                fstep.momentum().z()*+fstep.momentum().z())<<endl;
                cout<< 0.1*fstep.position()<<endl;
                // (Helix calculates in Tesla, GeV and cm)
                // (mu2e uses Tesla, MeV and mm)
                //Helix helix = Helix(0.001*fstep.momentum(),0.1*fstep.position(),-1,Bmagnet);
                //cout << "Curvature:           " << helix.getCurvature()<<endl;
                //cout << "Pt:[GeV/c]           " << (1.49898e-3*Bmagnet)/helix.getCurvature()<<endl;
                //cout << "Radius:[cm]          " << 0.5/helix.getCurvature()<<endl;
                //cout << "Helicity:            " << helix.getHelicity() << endl;
                //cout << "cotangent of theta:  " << helix.getCotTheta() << endl;
                //cout << "phi0:                " << helix.getPhi0() << endl;
                //cout << "d0:                  " << helix.getD0() << endl;
                //cout << "Z0:                  " << helix.getZ0() << endl;
                Pt_inval_si =  P_in_si.rho();
                P_inval_si  =  P_in_si.mag();                
                //_Pt_in_si->Fill(Pt_inval_si);
                //_P_in_si ->Fill(P_inval_si);
                //_Pz_in_si->Fill(P_in_si.z());
                P_out_si= lstep.momentum();   // momentum as the track leaves the tracker
                Pt_outval_si =  P_out_si.rho();
                P_outval_si  =  P_out_si.mag();
                //_Pt_out_si->Fill(Pt_outval_si);
                //_P_out_si ->Fill(P_outval_si);
                //_Pz_out_si->Fill(P_out_si.z());
                //nt[PINX] = P_in_si.getX();
                //nt[PINY] = P_in_si.getY();
                //nt[PINZ] = P_in_si.getZ();
                //nt[POUTX]= P_out_si.getX();
                //nt[POUTY]= P_out_si.getY();
                //nt[POUTZ]= P_out_si.getZ() ;
                // Loop over all StrawsHits to which this SimParticle contributed.
                double _timetodist=149.8962;
                //_hNEleHits->Fill(infos.size());
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
                    sid = str.id();
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
                    sid = str.id();
                    did = sid.getDeviceId();
                    const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
                    const CLHEP::Hep3Vector dirvec = str.getDirection();
                    double dt =hit.dt();
                    double disttomid =dt* _timetodist;
                    CLHEP::Hep3Vector hitpos  = mpvec+disttomid*dirvec;
                    //                    CLHEP::Hep3Vector hitposm = mpvec-disttomid*dirvec;
                    X_straw.push_back(hitpos.getX());
                    Y_straw.push_back(hitpos.getY());
                    Z_straw.push_back(hitpos.getZ());
                    cout<< "Device ID: "<<did<<endl;
                    cout<< "x:  "<<hitpos.getX()<<"  y:  "<<  hitpos.getY()<<"  z:  "<<hitpos.getZ()<<endl;
                    cout<<  MCPoint[did]<<endl;
                    mpoint mp;
                    mp.x=hitpos.getX();
                    mp.y=hitpos.getY();
                    mp.z=hitpos.getZ();
                    mp.errx=50.;            // don't know what the errors are yet
                    mp.erry=50.;            // 
                    
                    mpoints.push_back(mp);
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
                    //                    _Rdiff_s-> Fill(Rmc-Rhit);
                    //_Phidiff_s->Fill(smcpos.phi()-hitpos.phi());
                    X_res_straw.push_back(smcpos.getX()-hitpos.getX());
                    Y_res_straw.push_back(smcpos.getY()-hitpos.getY());
                    //_Xdiff_s -> Fill(smcpos.getX()-hitpos.getX());
                    //_Ydiff_s -> Fill(smcpos.getY()-hitpos.getY());                
                  }                        // end loop over hits
                //nt[NSTRAWS]= X_straw.size();
                if ( X_straw.size()>4)
                  {
                    if (FitCircle(X_straw, Y_straw))
                      {
                        //_x0y0_s->Fill(x0,y0);
                        //_R_rec_s->Fill(R_rec);
                        //_chi2_s -> Fill(chi2) ;
                        cout <<x0<<"  "<<y0<<"  "<< R_rec<<endl;
                        Pt =1000.*R_rec  * 2. * Bmagnet* Const;
                        //_Pt_rec_s->Fill(Pt);
                        //_Pt_diff_s->Fill(Pt-Pt_inval_si);
                        //nt[RREC_S] = R_rec;
                        cout<<"R_rec: " << R_rec<<endl;
                        //nt[PTREC_S]= Pt;
                        if (FitSinus2(Z_straw, R_straw))
                          {
                            //  _Pz_rec_s->Fill(Pz);
                            //_Pz_diff_s->Fill(Pz-P_in_si.z());
                            // Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                            //_P_rec_s->Fill(Ptot);
                            //_P_diff_s->Fill(Ptot-P_inval_si);
                            //  nt[PREC_SZ]= Pz;
                          }
                        //Double_t p[5] ={x0,y0,R_rec,0.,0.};
                        /*if (FitHelix(p))
                          {
                            nt[RREC_H] = p[2];
                            Double_t Pth =1000.*p[2]  * 2. * Bmagnet* Const;
                            nt[PTREC_h]= Pth;
                            cout<<"Pt in:    "<<Pt_inval_si 
                                <<"  Ptrec_s:  "<< Pt<<"  diff: " << Pt_inval_si-Pt
                                <<"  Ptrec_h:  "<< Pth<<"  diff: " << Pt_inval_si-Pth
                                <<"  Rrec_h:  "<< p[2]
                                <<endl;
                          }
                        */
                      }
                  }
              }
        }
      }           // end loop over simparticles 












    /*    for ( size_t cluster=0; cluster<clusters.size(); ++cluster) // Loop over StrawClusters
      {
        StrawCluster const& scluster = clusters.at(cluster);        
        cout << "averageT:  "<<sutils.averageT(scluster,evt)
             <<"  did: "<< sutils.did(scluster,evt)
             <<"  Station: "<<sutils.Station(scluster,evt)
             <<"  sector: "<< sutils.secid(scluster,evt)
             <<" X:  "<<sutils.midX(scluster,evt)
             <<" dirX:  "<<sutils.dirX(scluster,evt)
             <<endl;
        clusterbydid.insert(pair<int,StrawCluster>(sutils.did(scluster,evt),scluster));
        clusterbystation.insert(pair<int,StrawCluster>(sutils.Station(scluster,evt),scluster));
        cout<< "cluster Device ID: "<<sutils.did(scluster,evt)<<endl;
        cout<< "dTX:"<<sutils.dTX(scluster,evt)<<endl;
        cout<< "X:  "<<sutils.midX(scluster,evt)<<endl;
        cout<< "MC:  "<< MCPoint[sutils.did(scluster,evt)]<<endl;
      }
    */
    clusterbydid=sutils.clusterbydid(clusters,evt);
    clusterbystation=sutils.clusterbystation(clusters,evt);
    cout <<"did map size: " << clusterbydid.size()<<endl;
    cout <<"station map size: " << clusterbystation.size()<<endl;
    Int_t nint = 0;
    for (int i = 0;i<36;i++)
      {
        cout << " did:  "<<i << "  Count:  "<< clusterbydid.count(i)<<endl;
        if (clusterbydid.count(i)>1) 
          {
            pair<multimap<int,StrawCluster>::iterator, multimap<int,StrawCluster>::iterator> ppp1;
            ppp1 = clusterbydid.equal_range(i);
            //            multimap<int,StrawCluster>::iterator first1 = ppp1.first;
                        multimap<int,StrawCluster>::iterator first22 = ppp1.first;
            multimap<int,StrawCluster>::iterator last1 = ppp1.second;
            last1--;
            multimap<int,StrawCluster>::iterator last2 = ppp1.second;
            for (multimap<int,StrawCluster>::iterator first1 = ppp1.first;
                 first1  != last1; ++first1)
              {
                first22=first1;
                first22++;
                for (multimap<int,StrawCluster>::iterator first2 = first22;
                     first2 != last2;
                     ++first2)
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
                      X.push_back(intersection.x());
                      Y.push_back(intersection.y());
                      //                      Z.push_back(0.5*(junk.mpz+pjunk.mpz));
                      CLHEP::Hep3Vector junkX =sutils.midX(junk,evt);
                      CLHEP::Hep3Vector pjunkX=sutils.midX(pjunk,evt);
                      Z.push_back(0.5*(junkX.z()+pjunkX.z()));        
                      Double_t Rhit=intersection.r();
                      R.push_back(Rhit);
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
   _hNInter->Fill(X.size());
   // nt[NINT] = X.size();

    if (X.size()>4)
      {
        if (FitCircle(X, Y))
          {
            _x0y0->Fill(x0,y0);
            _R_rec->Fill(R_rec);
             _chi2 -> Fill(chi2) ;
             Pt =1000.*R_rec * 2. * Bmagnet* Const;
             _Pt_rec->Fill(Pt);
            //_Pt_diff->Fill(Pt-Pt_inval_si);
            //nt[RREC] = R_rec;
            //nt[PTREC]= Pt;
             if (FitSinus2(Z, R))
                          {
                _Pz_rec->Fill(Pz);
                //_Pz_diff->Fill(Pz-P_in_si.z());
                Double_t Ptot = TMath::Sqrt(Pz*Pz+Pt*Pt);
                _P_rec->Fill(Ptot);
                //_P_diff->Fill(Ptot-P_inval_si);
                //nt[PRECZ]= Pz;
                 }
          }
      }



   Int_t sint = 0;
    for (int i = 0;i<18;i++)
      {
        cout << " station:  "<<i << "  Count:  "<< clusterbystation.count(i)<<endl;
        if (clusterbystation.count(i)>1) 
          {
            pair<multimap<int,StrawCluster>::iterator, multimap<int,StrawCluster>::iterator> ppp1;
            ppp1 = clusterbystation.equal_range(i);
            multimap<int,StrawCluster>::iterator first22 = ppp1.first;
            multimap<int,StrawCluster>::iterator last1 = ppp1.second;
            last1--;
            multimap<int,StrawCluster>::iterator last2 = ppp1.second;
            for (multimap<int,StrawCluster>::iterator first1 = ppp1.first;
                 first1 != last1;
                 ++first1)
              {
                first22=first1;
                first22++;
                for (multimap<int,StrawCluster>::iterator first2 = first22;
                     first2 != last2;
                     ++first2)
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
                      sint ++;
                      break;
                    }  // end switch 
                } // end for first2
            }// end for first1
        }// end count >1
          // cout << "Number of elements with key: "<<i<<"  " << m.count(i) << endl;
          //pair<multimap< int,straw>::iterator, multimap<int,straw>::iterator> ppp;

      }   ///endloop over all devices
    cout<<"sint:  "<<nint<<endl;

    ++ncalls;  
  } // end of ::analyze.


  bool ReadStrawCluster::FitCircle(    vector<double> X,vector<double> Y)
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
      cout <<   "x[i]:  "<<X[i]
           <<   "  y[i]" <<Y[i]<<endl;
      ex[i] = 5.0 ; 
      ey[i] = 5.0 ;
    }
    error = new TGraphErrors(n,x,y,ex,ey);
    TMinuit *gmMinuit = new TMinuit(3); 
    gmMinuit->SetPrintLevel(normal);
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
 
  bool  ReadStrawCluster::FitSinus2( vector<double> X,vector<double> Y)
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
    int result=gmMinuit2->Migrad();
    cout << " Result: "<< result <<endl;
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


using mu2e::ReadStrawCluster;
DEFINE_ART_MODULE(ReadStrawCluster);
