//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no tracks (makes sense only for single CE events)
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TValidationModule.hh"

ClassImp(TValidationModule)
//-----------------------------------------------------------------------------
TValidationModule::TValidationModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;
  fFillDioHist     = 1;

  fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
}

//-----------------------------------------------------------------------------
TValidationModule::~TValidationModule() {
}


//-----------------------------------------------------------------------------
void TValidationModule::BookCaloHistograms(CaloHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------
  HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);

  for (int i=0; i<4; i++) {
    HBook1F(Hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
    HBook1F(Hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

    HBook1F(Hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

    HBook1F(Hist->fR700  [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRWE700[i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
  }
}

//-----------------------------------------------------------------------------
void TValidationModule::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),150, 0, 300,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TValidationModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TValidationModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 400,  90  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 104.5,105.5,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
  HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 400,  90  ,110. ,Folder);
  Hist->fPDio->Sumw2(kTRUE);

  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fChi2DofC   ,"chi2dc"   ,Form("%s: track chi2/N calc" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 200,-20  ,20 ,Folder);
  HBook1F(Hist->fChi2Match  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook2F(Hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
  HBook2F(Hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDucVsPath  ,"duc_vs_path",Form("%s: T-C Duc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvcVsPath  ,"dvc_vs_path",Form("%s: T-C Dvc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
  HBook2F(Hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

  HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
  HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
  HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
  HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
  HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
  HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

  HBook2F(Hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
  HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
  HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
  HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

  HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
}

//-----------------------------------------------------------------------------
void TValidationModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(Hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  //  char  name[200];
  for (int i=0; i<2; i++) {
    sprintf(name,"ncch_%i",i);
    HBook1F(Hist->fNCaloCrystalHits[i],name,Form("%s: N(calo crystal hits) [%i]",Folder,i),500,0,1000,Folder);
    sprintf(name,"ncch_vs_vane_%i",i);
    HBook2F(Hist->fNCaloHitsVsVane[i],name,Form("%s: N(calo crystal hits) vs vane[%i]",Folder,i),4,0,4,200,0,200,Folder);
    sprintf(name,"ncch_vs_row_%i",i);
    HBook2F(Hist->fNCaloHitsVsRow[i],name,Form("%s: N(calo crystal hits) vs row [%i]",Folder,i),20,0,20,200,0,200,Folder);
    sprintf(name,"ncch_vs_col_%i",i);
    HBook2F(Hist->fNCaloHitsVsCol[i],name,Form("%s: N(calo crystal hits) vs col [%i]",Folder,i),50,0,50,200,0,200,Folder);
  }

  for (int i=0; i<4; i++) {
    HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]",Folder,i), 300, 0,150,Folder);
    HBook2F(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,500,0,100,Folder);
    HBook2F(Hist->fNCrVsR      [i],Form("ncr_vs_r_%i",i),Form("%s: N Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,100,0,100,Folder);

    HBook2F(Hist->fNCrystalHitsVsR[i],Form("ncrh_vs_r_%i",i),Form("%s: N Crystal Hits[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
    HBook2F(Hist->fNHitCrystalsVsR[i],Form("nhcr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
  }

  HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  HBook1F(Hist->fECal,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
  HBook1F(Hist->fECalOverEKin,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
}

//-----------------------------------------------------------------------------
void TValidationModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
}
//_____________________________________________________________________________
void TValidationModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

//-----------------------------------------------------------------------------
// book crystal histograms
//-----------------------------------------------------------------------------
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
  book_event_histset[ 2] = 1;	        // events without reconstructed tracks
  book_event_histset[ 3] = 1;	        // events with a reconstructed cluster
  book_event_histset[ 4] = 1;	        // events without reconstructed clusters
  book_event_histset[ 5] = 1;	        // events w/o reconstructed tracks, |costh|<0.4
  book_event_histset[ 6] = 1;	        // events with tracks passing "Set C" cuts
  book_event_histset[ 7] = 1;	        // events with E(cluster) > 60 MeV
  book_event_histset[ 8] = 1;	        // events with the highest energy cluster on the 1st disk
  book_event_histset[ 9] = 1;	        // events with the highest energy cluster on the 2nd disk
  book_event_histset[10] = 0;	        // 
  book_event_histset[11] = 1;	        // selection cuts
  book_event_histset[12] = 1;	        // 
  book_event_histset[13] = 1;	        // 
  book_event_histset[14] = 1;	        // 
  book_event_histset[15] = 1;	        // 
  book_event_histset[16] = 1;	        // 
  book_event_histset[17] = 1;	        // 
  book_event_histset[18] = 1;	        // 
  book_event_histset[19] = 0;	        // 
  book_event_histset[20] = 0;	        // 
  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 
					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
  book_track_histset[  2] = 1;		// all tracks e- passing Set C cuts, events with clusters 
  book_track_histset[  3] = 1;		// all tracks e- passing Set C cuts, events w/o  clusters
  book_track_histset[  4] = 1;		// all tracks e- passing Set C cuts, events with clusters, no closest
  book_track_histset[  5] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
  book_track_histset[  6] = 1;		// all tracks e-  E/P > 0.4
  book_track_histset[  7] = 1;		// all tracks e- passing Set C cuts, E/P > 0.4
  book_track_histset[  8] = 1;		// Set C tracks e- , |xslope| < 3.
  book_track_histset[  9] = 1;		// all  tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 10] = 1;		// Set C tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 11] = 1;		// all tracks with P > 103.5
  book_track_histset[ 12] = 1;		// tracks with fcons < 1.e-4
  book_track_histset[ 13] = 1;		// "Set C" tracks with 100 <= P < 110 
  book_track_histset[ 14] = 1;		// tracks with fcons < 1.e-2
  book_track_histset[ 15] = 1;		// tracks intersecting the 1st disk
  book_track_histset[ 16] = 1;		// tracks intersecting the 2nd disk
  book_track_histset[ 17] = 1;		// tracks with no calorimeter intersections
  book_track_histset[ 18] = 1;		// Set C tracks, T0 > 700
  book_track_histset[ 19] = 1;		// Set C tracks with E/P > 0

  book_track_histset[ 20] = 1;		// tracks with Nhits >= 20
  book_track_histset[ 21] = 1;		// tracks with Nhits >= 20 and chi/Ndof < 3
  book_track_histset[ 22] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100
  book_track_histset[ 23] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) > 0 (interesting for muons)
  book_track_histset[ 24] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) < 0 (interesting for electrons)

  book_track_histset[ 25] = 1;		// Set C tracks, 0 < E/P < 1.15,  -2 < DT < 4, chi2(match) < 100, P>100
  book_track_histset[ 26] = 1;		// [25] + LLHR_CAL > 0 - interesting for muons
  book_track_histset[ 27] = 1;		// [25] + LLHR_CAL < 0 - interesting for electrons
  book_track_histset[ 28] = 1;		// Set C tracks, E/P > 1.1
  book_track_histset[ 29] = 1;		// Set C tracks, E/P > 0, P > 100 *precursor for TRK_25*

  book_track_histset[ 30] = 1;		// tracks with Nhits >= 25
  book_track_histset[ 31] = 1;		// tracks with Nhits >= 25 and chi/Ndof < 3

  book_track_histset[ 40] = 1;		// all tracks, alg_mask = 1
  book_track_histset[ 41] = 1;		// Set "C" tracks, alg_mask = 1
  book_track_histset[ 42] = 1;		// Set "C" tracks, alg_mask = 1, T > 700

  book_track_histset[ 50] = 1;		// all tracks, alg_mask = 2
  book_track_histset[ 51] = 1;		// Set "C" tracks, alg_mask = 2
  book_track_histset[ 52] = 1;		// Set "C" tracks, alg_mask = 2, T > 700

  book_track_histset[ 60] = 1;		// all tracks, alg_mask = 3
  book_track_histset[ 61] = 1;		// Set "C" tracks, alg_mask = 3
  book_track_histset[ 62] = 1;		// Set "C" tracks, alg_mask = 3, T > 700

  

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-
  book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
  book_cluster_histset[3] = 1;		// clusters in events w/track passing SetC cuts and |dt|<2.5ns 
  book_cluster_histset[4] = 1;		// clusters > 10 MeV
  book_cluster_histset[5] = 1;		// clusters > 60 MeV
  book_cluster_histset[6] = 1;		// clusters disk#0
  book_cluster_histset[7] = 1;		// clusters disk#1

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book calorimeter histograms
//-----------------------------------------------------------------------------
  int book_calo_histset[kNCaloHistSets];
  for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

  book_calo_histset[0] = 1;		// all crystals
  book_calo_histset[1] = 1;		// all crystals, e > 0
  book_calo_histset[2] = 1;		// all crystals, e > 0.1
  book_calo_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNCaloHistSets; i++) {
    if (book_calo_histset[i] != 0) {
      sprintf(folder_name,"cal_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalo[i] = new CaloHist_t;
      BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
//   book_genp_histset[1] = 1;		// all crystals, e > 0
//   book_genp_histset[2] = 1;		// all crystals, e > 0.1
//   book_genp_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TValidationModule::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th, dio_wt, xv, yv, rv, zv, p;
  double            e, m, r;
  TLorentzVector    mom;

  fParticle->Momentum(mom);

  p      = mom.P();

  cos_th = mom.Pz()/p;

  dio_wt = TStntuple::DioWeightAl(p);

  xv = fParticle->Vx()+3904.;
  yv = fParticle->Vy();
  rv = sqrt(xv*xv+yv*yv);
  zv = fParticle->Vz();

  Hist->fEleMom->Fill(p);
  Hist->fDioMom->Fill(p,dio_wt);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks->Fill  (fNTracks[0]);
  Hist->fNStrawHits[0]->Fill(fNStrawHits);
  Hist->fNStrawHits[1]->Fill(fNStrawHits);

  double emax   = -1;
  double t0_cls = -1;
  double dt     = 9999.;

  TStnCluster* cluster(0);
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  TStnTrack* track(0);
  if (fNTracks > 0) track = fTrackBlock->Track(0);

  if (cluster) {
    emax   = cluster->Energy();
    t0_cls = cluster->Time();
  }

  double t0_trk = -1;
  if (track) {
    t0_trk = track->fT0;
  }

  if (track && cluster) {
    dt = t0_cls-t0_trk;
  }

  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  TStrawHitData*  sh;
  int n_good_hits = 0;
  for (int i=0; i<fNStrawHits; i++ ) {
    sh = fStrawDataBlock->Hit(i);
    dt = t0_cls-sh->Time() + 15;
    Hist->fDtClS->Fill(dt);
    Hist->fSHTime->Fill(sh->Time());

    if (fabs(dt+15.)< 50) n_good_hits += 1;
  }

  Hist->fNGoodSH->Fill(n_good_hits);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(fNGenp);
//-----------------------------------------------------------------------------
// crystals - count crystals with E > 1MeV
//-----------------------------------------------------------------------------
  TCalHitData* cch;

  int n_cch_1mev = 0;

  if (fCalorimeterType == 1) {
//-----------------------------------------------------------------------------
// vane calorimeter
//-----------------------------------------------------------------------------
    int  nhits_vane[2][4], nhits_row [2][20], nhits_col[2][50];
    int  crystal_id, vane_id, local_id, vane_row, vane_col;

    for (int i=0; i<4; i++) {
      nhits_vane[0][i] = 0;
      nhits_vane[1][i] = 0;
    }
      
    for (int i=0; i<20; i++) {
      nhits_row[0][i] = 0;
      nhits_row[1][i] = 0;
    }

    for (int i=0; i<50; i++) {
      nhits_col[0][i] = 0;
      nhits_col[1][i] = 0;
    }
      
    for (int ic=0; ic<fNCalHits; ic++) {
      cch        = fCalDataBlock->CalHitData(ic);
      crystal_id = cch->ID();

      if (cch->Energy() > 1.) {
	n_cch_1mev += 1;
      }
      // for each crystal determine its row and column
      // the following is for vanes
      vane_id  = crystal_id/484.;
      local_id = crystal_id-vane_id*484;
      vane_row = local_id/44;
      vane_col = local_id-vane_row*44;
      
      nhits_vane[0][vane_id ] += 1;
      nhits_row [0][vane_row] += 1;
      nhits_col [0][vane_col] += 1;
      
      if (cch->Energy() > 1.) {
	nhits_row [1][vane_row] += 1;
	nhits_col [1][vane_col] += 1;
	nhits_vane[1][vane_id ] += 1;
      }
    }

    Hist->fNCaloCrystalHits[0]->Fill(fNCalHits);
    Hist->fNCaloCrystalHits[1]->Fill(n_cch_1mev);

    for (int iv=0; iv<4; iv++) {
      Hist->fNCaloHitsVsVane[0]->Fill(iv,nhits_vane[0][iv]);
      Hist->fNCaloHitsVsVane[1]->Fill(iv,nhits_vane[1][iv]);
    }

    for (int ir=0; ir<20; ir++) {
      Hist->fNCaloHitsVsRow[0]->Fill(ir,nhits_row[0][ir]);
      Hist->fNCaloHitsVsRow[1]->Fill(ir,nhits_row[1][ir]);
    }

    for (int ic=0; ic<50; ic++) {
      Hist->fNCaloHitsVsCol[0]->Fill(ic,nhits_col[0][ic]);
      Hist->fNCaloHitsVsCol[1]->Fill(ic,nhits_col[1][ic]);
    }
  }
  else if (fCalorimeterType == 2) {
//-----------------------------------------------------------------------------
// disk calorimeter
//-----------------------------------------------------------------------------
    int      ndisks, n_hit_crystals[4], n_hit_crystals_tot;
    double   etot[4];

    TCalHitData* hit;

    //    TDisk*       disk;
    //    TEvdCrystal* cr;

    ndisks = fDiskCalorimeter->NDisks();

    int   bin, hit_id, idisk, nhits, nhits_r[4][100], n_hit_crystals_r[4][100];

    for (int id=0; id<ndisks; id++) {
      n_hit_crystals[id] = 0;
      etot[id]           = 0;

      for (int ib=0; ib<100; ib++) {
	nhits_r         [id][ib] = 0;
	n_hit_crystals_r[id][ib] = 0;
      }
    }

    nhits = fCalDataBlock->NHits();

    for (int i=0; i< nhits; i++) {
      hit    = fCalDataBlock->CalHitData(i);

      hit_id = hit->ID();
      idisk  = fDiskCalorimeter->DiskNumber(hit_id);
      r      = fDiskCalorimeter->CrystalRadius(hit_id);
      e      = hit->Energy(); 

      etot          [idisk] += e;
      n_hit_crystals[idisk] += 1;

      Hist->fECrVsR[idisk]->Fill(r,e);
      Hist->fNCrVsR[idisk]->Fill(r,1);

      bin  = (int) (r/10.);

      nhits_r         [idisk][bin] += 1;
//-----------------------------------------------------------------------------
// this is not correct, one needs to check whether this crystal has been hit,
// for the moment, to get going, ignore that
//-----------------------------------------------------------------------------
      n_hit_crystals_r[idisk][bin] += 1;
    }

    n_hit_crystals_tot = 0;

    double ecal = 0;
    for (int id=0; id<ndisks; id++) {
      n_hit_crystals_tot += n_hit_crystals[id];
      ecal += etot[id];
//-----------------------------------------------------------------------------
// fill 'per-disk' histograms
//-----------------------------------------------------------------------------
      Hist->fETot[id]->Fill(etot[id]);

//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
      for (int ib=0; ib<100; ib++) {
	r = (ib+0.5)*10.;
	Hist->fNCrystalHitsVsR[id]->Fill(r,nhits_r         [id][ib]);
	Hist->fNHitCrystalsVsR[id]->Fill(r,n_hit_crystals_r[id][ib]);
      }
    }

    Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
    Hist->fECal->Fill(ecal);

    double ekin(-1.);
//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
    if (fSimp) {
      p    = fSimp->fStartMom.P();
      m    = 105.658; // in MeV
      ekin = sqrt(p*p+m*m)-m;
    }
    Hist->fECalOverEKin->Fill(ecal/ekin);
  }
}

//-----------------------------------------------------------------------------
void TValidationModule::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  // determine crystal coordinates

  TDisk* disk = Cr->Disk();

  int idisk = disk->SectionID();
  // time needs to be defiend
  //    t  = Cr->Time();
  e     = Cr->Energy();
  r     = Cr->Radius();
  nhits = Cr->NHits();

  Hist->fVaneID->Fill(idisk);

  Hist->fEnergy  [idisk]->Fill(e);
  Hist->fNHits   [idisk]->Fill(nhits);
  //    Hist->fTime    [idisk]->Fill(t);
  Hist->fRadius  [idisk]->Fill(r);
  Hist->fRadiusWE[idisk]->Fill(r,e);
    
  e700 = 0;
  n700 = 0;
  for (int i=0; i<nhits; i++) {
    hit  = Cr->CalHitData(i);
    t   = hit->Time();
    Hist->fTime[idisk]->Fill(t);
    if (t > 700.) {
      n700 += 1;
      e700 += hit->Energy();
      Hist->fT700[idisk]->Fill(t);
    }
  }

  Hist->fE700   [idisk]->Fill(e700);
  Hist->fN700   [idisk]->Fill(n700);

  if (n700 > 0) {
    Hist->fR700  [idisk]->Fill(r);
    Hist->fRWE700[idisk]->Fill(r,e700);
  }
}


//-----------------------------------------------------------------------------
void TValidationModule::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fVaneID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY->Fill(Cluster->fSigY);
  Hist->fSigZ->Fill(Cluster->fSigZ);
  Hist->fSigR->Fill(Cluster->fSigR);
  Hist->fNCr0->Fill(Cluster->fNCrystals);
  Hist->fNCr1->Fill(Cluster->fNCr1);
  Hist->fFrE1->Fill(Cluster->fFrE1);
  Hist->fFrE2->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TValidationModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
}

//-----------------------------------------------------------------------------
void TValidationModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TValidationModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

  TLorentzVector  mom;
  double          chi2c, r;
  int             itrk;
  TrackPar_t*     tp;
					// pointer to local track parameters
  itrk = Track->Number();
  tp   = fTrackPar+itrk;

  Hist->fP[0]->Fill (Track->fP);
  Hist->fP[1]->Fill (Track->fP);
  Hist->fP[2]->Fill (Track->fP);
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fPDio->Fill(Track->fP,tp->fDioWt);

  Hist->fFitMomErr->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  Hist->fDpFront ->Fill(tp->fDpF);
  Hist->fDpFront0->Fill(tp->fDp0);
  Hist->fDpFront2->Fill(tp->fDp2);
  Hist->fDpFSt   ->Fill(tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,tp->fDpF);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fNDof->Fill(Track->NActive()-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));
  Hist->fNActive->Fill(Track->NActive());
  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  //  printf("TValidationModule::FillTrackHistograms: track charge is not defined yet\n");
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fAlgMask->Fill(Track->AlgMask());

  chi2c = Track->fChi2C/(Track->NActive()-5.);
  Hist->fChi2DofC->Fill(chi2c);

  //  int nh, nst_with_nh[10];
					// 2014-04-29: currently not saved

  //  for (int i=0; i<10; i++) nst_with_nh[i] = 0;

//   for (int i=0; i<40; i++) {
//     Hist->fNHVsStation->Fill(i,Track->fNHPerStation[i]);
//     nh = Track->fNHPerStation[i];
//     if ((nh >= 0) && (nh < 10)) {
//       nst_with_nh[nh] += 1;
//     }
//     else {
//       printf(">>> ERROR : nh = %20i, IGNORE \n",nh);
//     }
//  }

//   for (int i=0; i<10; i++) {
//     Hist->fNHVsNSt->Fill(i,nst_with_nh[i]);
//   }
  //-----------------------------------------------------------------------------
  // track-cluster matching part: 
  // - for residuals, determine intersection with the most energetic cluster
  // - for track -only parameters use intersection with lowest trajectory length
  //-----------------------------------------------------------------------------
  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-only
  //  TStnTrack::InterData_t*    vr = Track->fVMaxEp; // residuals

  if (vt) {
    Hist->fVaneID->Fill(vt->fID  );
    Hist->fXTrk->Fill  (vt->fXTrk);
    Hist->fYTrk->Fill  (vt->fYTrk);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    Hist->fRTrk->Fill  (r);

    Hist->fZTrk->Fill  (vt->fZTrk);
  }
  else {
//-----------------------------------------------------------------------------
// fill histograms with numbers easy to recognize as dummy
//-----------------------------------------------------------------------------
    Hist->fVaneID->Fill(-1.);
    Hist->fXTrk->Fill  (999.);
    Hist->fYTrk->Fill  (999.);
    Hist->fRTrk->Fill  (999.);
    Hist->fZTrk->Fill  (-1. );
  }

//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
  double ekin(-1.);
  if (fSimp) {
    double p, m;
    //    p    = fSimp->fStartMom.P();
    p = Track->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  Hist->fECl->Fill(tp->fEcl);
  Hist->fEClEKin->Fill(tp->fEcl/ekin);
  Hist->fEp->Fill(tp->fEp);
  Hist->fEpVsPath->Fill(tp->fPath,tp->fEp);

  Hist->fDx->Fill(tp->fDx);
  Hist->fDy->Fill(tp->fDy);
  Hist->fDz->Fill(tp->fDz);

  Hist->fDt->Fill(tp->fDt);
  Hist->fChi2Match->Fill(tp->fChi2Match);

  Hist->fDu->Fill    (tp->fDu);
  Hist->fDv->Fill    (tp->fDv);
  Hist->fDvVsDu->Fill(tp->fDu,tp->fDv);

  Hist->fPath->Fill(tp->fPath);
  Hist->fDuVsPath->Fill(tp->fPath,tp->fDu);

  //  double cu[4] = { -60.3049, -0.749111,    0.00522242,  -7.52018e-06};
  double cu[4] = { -59.5174, -0.541226, 0.00414309, -5.84989e-06 };
  double cv[4] = {  6.44161, 0.0722353, -0.000653084, 1.14054e-06};

  double x = tp->fPath;
  double corr_u = cu[0]+cu[1]*x+cu[2]*x*x+cu[3]*x*x*x;
  double corr_v = cv[0]+cv[1]*x+cv[2]*x*x+cv[3]*x*x*x;

  //  double duc = tp->fDu-0.34*(tp->fPath-350.);
  double duc = tp->fDu-corr_u;
  double dvc = tp->fDv-corr_v;
  
  Hist->fDucVsPath->Fill(tp->fPath,duc);
  Hist->fDvcVsPath->Fill(tp->fPath,dvc);

  Hist->fDvVsPath->Fill(tp->fPath,tp->fDv);
  Hist->fDtVsPath->Fill(tp->fPath,tp->fDt);

  Hist->fDuVsTDip->Fill(Track->fTanDip,tp->fDu);
  Hist->fDvVsTDip->Fill(Track->fTanDip,tp->fDv);

  Hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fRSlope->Fill(Track->RSlope());
  Hist->fXSlope->Fill(Track->XSlope());

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
  Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

  llhr_cal = Track->LogLHRCal();
  Hist->fLogLHRCal->Fill(llhr_cal);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  Hist->fEpVsDt->Fill(tp->fDt,tp->fEp);
  Hist->fLogLHRDeDx->Fill(llhr_dedx);
  Hist->fLogLHRXs->Fill(llhr_xs);
  Hist->fLogLHRTrk->Fill(llhr_trk);
  Hist->fLogLHR->Fill(llhr);

  Hist->fPdgCode->Fill(Track->fPdgCode);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  Hist->fNEPlVsNHPl->Fill(tp->fNEPl,tp->fNHPl);
  Hist->fNDPlVsNHPl->Fill(tp->fNDPl,tp->fNHPl);
  Hist->fChi2dVsNDPl->Fill(tp->fNDPl,Track->Chi2Dof());
  Hist->fDpFVsNDPl  ->Fill(tp->fNDPl,tp->fDpF);

  float        fre1(-1), fre2(-1);
  int          icl;
  TStnCluster* cl;

  if (Track->fVMinS) {
    icl = Track->fVMinS->fClusterIndex;
    if (icl >= 0) {
      cl = fClusterBlock->Cluster(icl);
      fre1 = cl->fFrE1;
      fre2 = cl->fFrE2;
    }
  }

  Hist->fFrE1->Fill(fre1);
  Hist->fFrE2->Fill(fre2);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TValidationModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrackBlock"    ,"TStnTrackBlock"   ,&fTrackBlock  );
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"  ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("StrawDataBlock","TStrawDataBlock"  ,&fStrawDataBlock);
  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"     ,"TSimpBlock"       ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);
//-----------------------------------------------------------------------------
// initialize likelihood histograms
// TRK 19: "Set C" plus reconstructed and matched cluster
//-----------------------------------------------------------------------------
//   const char   *mu2e_hist_dir, *muo_template_fn, *ele_template_fn;
//   //  char         fn[200];

//   mu2e_hist_dir   = gEnv->GetValue("mu2e.HistDir",gSystem->Getenv("MU2E_HIST_DIR"));
//   ele_template_fn = gEnv->GetValue("mu2e.EleTemplates",gSystem->Getenv("_none_"));
//   muo_template_fn = gEnv->GetValue("mu2e.MuoTemplates",gSystem->Getenv("_none_"));

//   //   sprintf(fn,"%s/v4_2_1/e00s1212.track_ana.hist",mu2e_hist_dir);

//   TH1* he_dt = gh1(ele_template_fn,"TrackAna","trk_19/dt");
//   TH1* he_ep = gh1(ele_template_fn,"TrackAna","trk_19/ep");

//   fLogLH->SetEleDtHist(he_dt);
//   fLogLH->SetEleEpHist(he_ep);

//   //  sprintf(fn,"%s/v4_2_1/m00s1212.track_ana.hist",mu2e_hist_dir);

//   TH1* hm_dt = gh1(muo_template_fn,"TrackAna","trk_19/dt");
//   TH1* hm_ep = gh1(muo_template_fn,"TrackAna","trk_19/ep");

//   fLogLH->SetMuoDtHist(hm_dt);
//   fLogLH->SetMuoEpHist(hm_ep);

//   //  sprintf(fn,"%s/v4_2_1/e00s1212.track_ana.hist",mu2e_hist_dir);

//   TH1* he_xs = gh1(ele_template_fn,"TrackAna","trk_1/xslope");
//   fLogLH->SetEleXsHist(he_xs);

//   //  sprintf(fn,"%s/v4_2_1/m00s1212.track_ana.hist",mu2e_hist_dir);

//   TH1* hm_xs = gh1(muo_template_fn,"TrackAna","trk_1/xslope");
//   fLogLH->SetMuoXsHist(hm_xs);

  fLogLH->Init("v4_2_4");

  return 0;
}


//_____________________________________________________________________________
int TValidationModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TValidationModule::FillHistograms() {

  double       cos_th (-2.),  cl_e(-1.);
  int          disk_id(-1), alg_mask, nsh, nactive;
  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  TStnCluster  *cl0;

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  if (fNClusters > 0) {
    cl0     = fClusterBlock->Cluster(0);
    cl_e    = cl0->Energy();
    disk_id = cl0->DiskID();
  }
  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
  else                FillEventHistograms(fHist.fEvent[2]);

  if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3]);
  else                FillEventHistograms(fHist.fEvent[4]);

  if ((fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
    FillEventHistograms(fHist.fEvent[5]); 
  }

  if (fNGoodTracks > 0) {
    FillEventHistograms(fHist.fEvent[6]); 

    TLorentzVector    mom;
    
    fParticle->Momentum(mom);

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fNTracks[0]));
    }
  }

  if (cl_e > 60.) {
    FillEventHistograms(fHist.fEvent[7]); 
    if (GetDebugBit(34)) {
      if (fNTracks[0] <= 0) {
	GetHeaderBlock()->Print(Form(" bit:034 cl_e = %10.3f",cl_e));
      }
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8]);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9]);
//-----------------------------------------------------------------------------
// Dave's ladder for all tracks
// 1. N(straw hits) > 20
//-----------------------------------------------------------------------------
  if (fSimp) {
    nsh    = fSimp->NStrawHits();
    pfront = fSimp->fMomTrackerFront;
  }
  else {
    nsh    = -1;
    pfront = -1.e6;
  }
  
  if (nsh >= 20) {
    FillEventHistograms(fHist.fEvent[11]);
    if (pfront > 100.) {
      FillEventHistograms(fHist.fEvent[12]);
      
      ce_pitch = 0.7; // kludge
      if ((ce_pitch > 0.577) && (ce_pitch < 1.)) {
	FillEventHistograms(fHist.fEvent[13]);

	if (fNTracks[0] > 0) {
	  FillEventHistograms(fHist.fEvent[14]);

					// here we have a track reconstructed

	  TStnTrack* trk = fTrackBlock->Track(0);

	  fcons = trk->fFitCons;
	  t0    = trk->T0();
	  reco_pitch = trk->fTanDip;
	  sigp       = trk->fFitMomErr;
	  sigt       = trk->fT0Err;
	  nactive    = trk->fNActive;
	  p          = trk->fP;
					// fit quality
	  if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	    FillEventHistograms(fHist.fEvent[15]);
	    if (t0 > 700) {
	      FillEventHistograms(fHist.fEvent[16]);
	      if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		FillEventHistograms(fHist.fEvent[17]);
		if (p > 103.5) {
		  FillEventHistograms(fHist.fEvent[18]);
		}
	      }
	    }
	  }

	  alg_mask = trk->AlgMask();

	  if ((alg_mask == 1) || (alg_mask == 3)) {
//-----------------------------------------------------------------------------
// track reconstructed with TrkPatRec 
//-----------------------------------------------------------------------------
	    FillEventHistograms(fHist.fEvent[24]);
	    if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	      FillEventHistograms(fHist.fEvent[25]);
	      if (t0 > 700) {
		FillEventHistograms(fHist.fEvent[26]);
		if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		  FillEventHistograms(fHist.fEvent[27]);
		  if (p > 103.5) {
		    FillEventHistograms(fHist.fEvent[28]);
		  }
		}
	      }
	    }
	  }
	  else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// track reconstructed with CalPatRec, but not with TrkPatRec
//-----------------------------------------------------------------------------
//	    int x=0;
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// the same ladder for TrkPatRec tracks 
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk);

      if (GetDebugBit(32) && (trk->fVMinS != NULL)) {
	if (trk->fVMinS->fChi2Match > 100) {
	  GetHeaderBlock()->Print(Form("bit032: chi2(match) = %10.3lf",trk->fVMinS->fChi2Match));
	}
      }

      if (GetDebugBit(35) && (trk->fP > 106.)) {
	GetHeaderBlock()->Print(Form("bit035: P = %10.3lf",trk->fP));
      }

					// events with at least one  cluster
      if (fNClusters > 0) {
	FillTrackHistograms(fHist.fTrack[2],trk);
      }
      else {
					// events without a cluster
	FillTrackHistograms(fHist.fTrack[3],trk);
      }
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	FillTrackHistograms(fHist.fTrack[4],trk);
      }
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	FillTrackHistograms(fHist.fTrack[5],trk);
      }
//-----------------------------------------------------------------------------
// TRK 8: good track, |xslope| < 3
//-----------------------------------------------------------------------------
      if (fabs(trk->XSlope()) < 3.) {
	FillTrackHistograms(fHist.fTrack[8],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
// TRK 7 : events with a "Set C" track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
    if ((trk->Ep() > 0.4) && ( trk->Ep() < 1.2)) {
      FillTrackHistograms(fHist.fTrack[6],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[7],trk);
      }
    }
//----------------------------------------------------------------------------
//TRK  9: events with track and with no EM cluster      with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//----------------------------------------------------------------------------
    if (cl_e > 60) {
      FillTrackHistograms(fHist.fTrack[9],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[10],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 103.5 MeV
//-----------------------------------------------------------------------------
    if (trk->P() > 103.5) FillTrackHistograms(fHist.fTrack[11],trk);
//-----------------------------------------------------------------------------
// TRK 12: tracks with fcon < 1e-4
// TRK 13: "Set C" tracks with 100 <= P < 110 
// TRK 14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------
    if (trk->fFitCons < 1.e-4) FillTrackHistograms(fHist.fTrack[12],trk);

    if ((trk->fIDWord == 0) && (trk->P() >= 100.) && (trk->P() < 110.)) {
      FillTrackHistograms(fHist.fTrack[13],trk);
    }

    if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk);

    TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only
//-----------------------------------------------------------------------------
// TRK 15: tracks which have intersection with the 1st disk
// TRK 16: tracks which have intersection with the 2nd disk
// TRK 17: tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
    if (vt) {
      if      (vt->fID == 0) {
	FillTrackHistograms(fHist.fTrack[15],trk);
      }
      else if (vt->fID == 1) {
	FillTrackHistograms(fHist.fTrack[16],trk);
      }
    }
    else {
      FillTrackHistograms(fHist.fTrack[17],trk);
    }
//-----------------------------------------------------------------------------
// TRK 18: Set "C" tracks with T0 > 700
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->T0() > 700)) {
      FillTrackHistograms(fHist.fTrack[18],trk);
    }
//-----------------------------------------------------------------------------
// TRK 19: Set "C" tracks with an associated cluster
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->Ep() > 0)) {
      FillTrackHistograms(fHist.fTrack[19],trk);

      if (trk->LogLHRCal() < 0) {
	if (GetDebugBit(29)) {
	  GetHeaderBlock()->Print(Form(" bit:029 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f chi2_tcm = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt(), trk->fVMinS->fChi2Match));
	}
      }

      if (tp->fDu < -80.) {
	if (GetDebugBit(33)) {
	  GetHeaderBlock()->Print(Form(" bit:033 DU = %10.3f dv = %10.3f ep = %10.3f dt = %10.3f",
				       tp->fDu, tp->fDv, tp->fEp, trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK 20: tracks with >= 20 hits
// TRK 21: tracks with >= 20 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->fNActive >= 20) {
      FillTrackHistograms(fHist.fTrack[20],trk);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[21],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 22: Set "C" tracks with an associated cluster and chi2(match) < 100
// TRK 23: Set "C" tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) > 0
//         this is interesting to see which muons are getting misidentified
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (tp->fEp > 0) && (tp->fChi2Match < 100.)) {
      FillTrackHistograms(fHist.fTrack[22],trk);
      if    (trk->LogLHRCal() > 0) {
	FillTrackHistograms(fHist.fTrack[23],trk);
	
	if (trk->fP < 80.) {
	  if (GetDebugBit(36)) {
	    GetHeaderBlock()->Print(Form(" bit:036 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	  }
	}
      }
      else {
//-----------------------------------------------------------------------------
// TRK 24: Set "C" tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) < 0
//         this set allows to see which electrons are getting misidentified
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[24],trk);
      }
//-----------------------------------------------------------------------------
// TRK 25: Set "C" tracks, 0 < E/p < 1.15,  |dt_corr| < 3, chi2(match) < 100, P>100
//-----------------------------------------------------------------------------
      double dt_corr = trk->Dt(); // -1.;
      if ( (fabs(dt_corr) < 3.) && (tp->fEp < 1.15) && (trk->fP > 100.)) {
	FillTrackHistograms(fHist.fTrack[25],trk);
//-----------------------------------------------------------------------------
// more details on the calorimeter-based likelihood 
// TRK 26 : TRK 25 events with LLHR_CAL > 0 ( interesting for muons)
// TRK 27 : TRK 25 events with LLHR_CAL < 0 ( interesting for electrons)
//-----------------------------------------------------------------------------
	double llhr_cal = trk->LogLHRCal();
	if (llhr_cal > 0) {
	  FillTrackHistograms(fHist.fTrack[26],trk);

	  if (llhr_cal > 5.) {
	    if (GetDebugBit(37)) {
	      GetHeaderBlock()->Print(Form(" bit:037 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	    }
	  }
	}
	else {
	  FillTrackHistograms(fHist.fTrack[27],trk);
	}
      }

    }
//-----------------------------------------------------------------------------
// TRK 28 : events with a "Set C" track and a cluster E/P > 1.1
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if (trk->Ep() > 1.1) {
	FillTrackHistograms(fHist.fTrack[28],trk);

	if (GetDebugBit(28)) {
	  GetHeaderBlock()->Print(Form(" bit:028 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK 29 : events with a "Set C" track, E/P>0 and P>100 : precursor for TRK 25
// in effect startign from 2014-06-17 10:31am
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if ((tp->fEp > 0) && (trk->fP > 100.)) {
	FillTrackHistograms(fHist.fTrack[29],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 30: tracks with >= 25 hits
// TRK 31: tracks with >= 25 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->fNActive >= 25) {
      FillTrackHistograms(fHist.fTrack[30],trk);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[31],trk);
      }
    }
//-----------------------------------------------------------------------------
// split tracks by the algorithm mask: 1 , 2 , or 3
//-----------------------------------------------------------------------------
    alg_mask = trk->AlgMask();
    if      (alg_mask == 1) {
//-----------------------------------------------------------------------------
// TrkPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[40],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[41],trk);
	// print run/event numbers :
	if (GetDebugBit(6)) {
	  double ep = trk->Ep();
	  if ((ep > 0.8) && (ep < 1.1)) {
	    GetHeaderBlock()->Print(Form(" bit:006 trk_41: track E/P = %8.3f",trk->Ep()));
	  }
	}
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[42],trk);
      }
    }
    else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// CalPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[50],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[51],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[52],trk);
      }
    }
    else if (alg_mask == 3) {
//-----------------------------------------------------------------------------
// TrkPatRec+CalPatRec tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[60],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[61],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[62],trk);
      }
    }
  }
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
  TStnCluster* cl;
  int id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()    > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl);
  }
//-----------------------------------------------------------------------------
// calorimeter histograms
//-----------------------------------------------------------------------------
  TDisk*         disk;
  TStnCrystal*   cr;

  if (fCalorimeterType == 2) {
    int nd = fDiskCalorimeter->NDisks();

    for (int i=0; i<nd; i++) {
      disk = fDiskCalorimeter->Disk(i);
      for (int ic=0; ic<disk->NCrystals(); ic++) {
	cr = disk->Crystal(ic);
	FillCaloHistograms(fHist.fCalo[0],cr);

	if (cr->Energy() > 0) {
	  FillCaloHistograms(fHist.fCalo[1],cr);
	}
	if (cr->Energy() > 0.1) {
	  FillCaloHistograms(fHist.fCalo[2],cr);
	}
	if (cr->Energy() > 1.0) {
	  FillCaloHistograms(fHist.fCalo[3],cr);
	}
      }
    }
  }
  //-----------------------------------------------------------------------------
  // radial distributions for crystals
  //-----------------------------------------------------------------------------
  static int first_entry(1);

  if (first_entry == 1) {
    if (fCalorimeterType == 2) {
      int nd = fDiskCalorimeter->NDisks();
	
      for (int i=0; i<nd; i++) {
	disk = fDiskCalorimeter->Disk(i);
	for (int ic=0; ic<disk->NCrystals(); ic++) {
	  cr = disk->Crystal(ic);

	  fHist.fCrystalR[i]->Fill(cr->Radius());
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
  first_entry = 0;
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TValidationModule::Event(int ientry) {

  double                xs, p;
  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TLorentzVector        mom;

  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  //  fStrawDataBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();

  TGenParticle* genp;
  int           pdg_code, generator_code;

  fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fParticle = genp;
      break;
    }
  }
					// may want to revisit the definition of fSimp
  fSimp     = fSimpBlock->Particle(0);

  fParticle->Momentum(mom);
					// this is a kludge, to be removed at the next 
					// ntupling 
  //  fEleE     = fParticle->Energy();
  p         = mom.P();
  fEleE     = sqrt(p*p+0.511*0.511);


  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  fNTracks[0] = fTrackBlock->NTracks();
  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
  fNStrawHits = fStrawDataBlock->NHits();

  fDiskCalorimeter->InitEvent(fCalDataBlock);

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  int ntrk = fNTracks[0];

  TrackPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
					// assume less 20 tracks
    tp             = fTrackPar+itrk;

    track          = fTrackBlock->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
      }
    }
//-----------------------------------------------------------------------------
// process hit masks
//-----------------------------------------------------------------------------
    int i1, i2, n1(0) ,n2(0), ndiff(0);
    int nbits = track->fHitMask.GetNBits();
    for (int i=0; i<nbits; i++) {
      i1 = track->HitMask()->GetBit(i);
      i2 = track->ExpectedHitMask()->GetBit(i);
      n1 += i1;
      n2 += i2;
      if (i1 != i2) ndiff += 1;
    }
//-----------------------------------------------------------------------------
// define additional parameters
//-----------------------------------------------------------------------------
    tp->fNHPl = n1;
    tp->fNEPl = n2;
    tp->fNDPl = ndiff;

    tp->fDpF   = track->fP     -track->fPFront;
    tp->fDp0   = track->fP0    -track->fPFront;
    tp->fDp2   = track->fP2    -track->fPFront;
    tp->fDpFSt = track->fPFront-track->fPStOut;

    if (fFillDioHist == 0) tp->fDioWt = 1.;
    else                   tp->fDioWt = TStntuple::DioWeightAl(fEleE);
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Match = -1.e6;
    tp->fPath      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt - 0.22; // - 1.;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
      tp->fPath      = vr->fPath;
    }

    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
      GetHeaderBlock()->Print(Form(" TValidationModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f\n ",tp->fEp,track->fEp));
    }
//-----------------------------------------------------------------------------
// PID likelihoods
//-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    if (GetDebugBit(7)) {
      if ((id_word == 0) && (llhr_cal > 20)) {
	GetHeaderBlock()->Print(Form("bit:007: dt = %10.3f ep = %10.3f",track->Dt(),tp->fEp));
      }
    }

    if (GetDebugBit(8)) {
      if ((id_word == 0) && (llhr_cal < -20)) {
	GetHeaderBlock()->Print(Form("bit:008: p = %10.3f dt = %10.3f ep = %10.3f",
				     track->P(),track->Dt(),tp->fEp));
      }
    }

    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
  }

  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  //  fDiskCalorimeter->InitEvent(fCalDataBlock);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TValidationModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fDpF > 1.) {
	GetHeaderBlock()->Print(Form("pF pRec, fDpf = %10.3f  %10.3f  %10.3f",
				     trk->fPFront, trk->Momentum()->P(),tp->fDpF));
      }
    }
//-----------------------------------------------------------------------------
// bit 9: Set C tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(9) == 1) {
      double ep = trk->Ep();
      if (trk->fIDWord == 0) { 
	if (((ep > 0.42) && (ep < 0.46)) || ((ep > 0.35) && (ep < 0.39))) {
	  GetHeaderBlock()->Print(Form("bit:009 ep = %10.3f e = %10.3f p = %10.3f",
				       trk->fEp,trk->fEp*trk->fP,trk->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 10: Set C tracks with Ecl > 80
//-----------------------------------------------------------------------------
    if (GetDebugBit(10) == 1) {
      double ecl = trk->ClusterE();
      if (trk->fIDWord == 0) { 
	if (ecl > 60) {
	  GetHeaderBlock()->Print(Form("bit:010 e = %10.3f p = %10.3f",
				       ecl,trk->fP));
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if (ntrk > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %i5",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TValidationModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TValidationModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

