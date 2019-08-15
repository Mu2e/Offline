///////////////////////////////////////////////////////////////////////////////
// cloned from Vadim's code
//
// 2015-07-10 P.Murat: default condiguration is stored in ParticleID/fcl/prolog.fcl
//
// assume that both electron and muon track reconstruction have been attempted,
// so expect two track collections on input
// need calculation od dEdX prob to be symmetric over the electron and muon hypotheses
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>
#include <sstream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

//ROOTs
#include "TH1F.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/PIDProduct.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "BTrkData/inc/Doublet.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/TrkFitDirection.hh"

#include "ParticleID/inc/PIDUtilities.hh"
#include "RecoDataProducts/inc/AvikPIDProductCollection.hh"

#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"

#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"

using CLHEP::Hep3Vector;

namespace mu2e {

  class DoubletAmbigResolver;

  class AvikPID : public art::EDProducer {

  private:

    enum { nbounds = 11 };
    static float _pathbounds[nbounds];

    AvikPIDProduct _pid;

    int    _debugLevel;
    int    _verbosity;
    int    _diagLevel;

    int    _processed_events;
    int    _evtid;

    int    _ele_trkid;
    int    _muo_trkid;

    string _trkPatRecDemModuleLabel;
    string _trkPatRecDmmModuleLabel;

    string _eleTemplates;
    string _muoTemplates;

    string _eleDedxTemplateFile;
    string _muoDedxTemplateFile;

    ProditionsHandle<Mu2eDetector> _mu2eDetector_h;

    TrkParticle     _tpart;
    TrkFitDirection _fdir;
    std::string     _iname; // data instance name

    TH1D* _heletemp[nbounds];
    TH1D* _hmuotemp[nbounds];

    int   _templatesnbins ;
    float _templateslastbin ;
    float _templatesbinsize ;

    const KalRepPtrCollection* _listOfEleTracks;
    const KalRepPtrCollection* _listOfMuoTracks;

    int    _trkid;
    double _trkmom;

    double _logDedxProbEle;
    double _logDedxProbMuo;
//-----------------------------------------------------------------------------
// Vadim's global slopes
//-----------------------------------------------------------------------------
    double _drdsVadimEle;
    double _drdsVadimEleErr;

    double _drdsVadimMuo;
    double _drdsVadimMuoErr;
//-----------------------------------------------------------------------------
// Avik's sums
//-----------------------------------------------------------------------------
    double _sumAvikEle;
    double _sumAvikMuo;

    int    _nMatched;
    int    _nMatchedAll;
    double _sq2AvikEle;
    double _sq2AvikMuo;
//-----------------------------------------------------------------------------
// same-sign slopes
//-----------------------------------------------------------------------------
    int    _ele_nusedSsH;  // Nhits used to calculate the SS slopes
    int    _muo_nusedSsH;
    double _drdsSsEle;
    double _drdsSsEleErr;
    double _drdsSsMuo;
    double _drdsSsMuoErr;
    double _logRatioSs;
//-----------------------------------------------------------------------------
// Avik doesn't calculate OS slopes! ...
//-----------------------------------------------------------------------------
    double   _drdsOsEle;
    double   _drdsOsEleErr;
    double   _drdsOsMuo;
    double   _drdsOsMuoErr;
    int      _ele_nusedOsH;             // Nhits used to calculate the OS slopes
    int      _muo_nusedOsH;
    double   _logRatioOs;

    double   _ele_resSumOs;             // d(dxdz)^alpha sums
    double   _muo_resSumOs;
    int      _ele_nusedOsD;             // Ndoublets used to calculate the sums
    int      _muo_nusedOsD;
//-----------------------------------------------------------------------------
// power coefficients
//-----------------------------------------------------------------------------
    double   _pow1;
    double   _pow2;
    double   _bound;

    double   _maxDeltaDxDzOs;

    TMinuit*               _minuit;

    fhicl::ParameterSet    _darPset;         // parameter set for doublet ambig resolver
    DoubletAmbigResolver*  _dar;

    TTree*                 _pidtree;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit AvikPID(fhicl::ParameterSet const& pset);
    virtual ~AvikPID();

    void beginJob();
    void beginRun(art::Run &run);
    void beginSubRun(art::SubRun & lblock );
    virtual void produce(art::Event& event);
    void endJob();


    static  void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);
    static  int  findlowhist(float d);

    bool calculateVadimSlope(const KalRep* KRep, double *Slope, double *Eslope);

    double calculateDedxProb(std::vector<double>* GasPaths ,
                             std::vector<double>* EDeps    ,
                             TH1D**               Templates);

    void   doubletMaker(const KalRep* ele_Trk, const KalRep* muo_Trk);

    //    double calculateAvikSums();

    int    CalculateSlope(vector<double>& Fltlen, vector<double>& Resid,
                          double&         Slope , double&         SlopeErr);

    int    AddHits(const Doublet* Multiplet, vector<double>& Fltlen, vector<double>& Resid);

    int    AddSsMultiplets(const vector<Doublet>* ListOfDoublets,
                           vector<double>&        Fltlen        ,
                           vector<double>&        Resid         );

    int    AddOsMultiplets(const vector<Doublet>* ListOfDoublets,
                           vector<double>&        Fltlen        ,
                           vector<double>&        Resid         );

    void   calculateSsSums(const vector<Doublet>* ListOfDoublets, double& Drds, double& DrdsErr, int& NUsed);

    void   calculateOsSums(const vector<Doublet>* ListOfDoublets,
                           double& Drds, double& DrdsErr, int& NUsedHits,
                           double& Sum , int& NUsedDoublets);

    double weightedResidual     (double Dr);

    double weightedSlopeResidual(double Dr);

    // Save directory from beginJob so that we can go there in endJob.
    //    TDirectory* _directory;


  };


  TGraphErrors *error;

  float AvikPID::_pathbounds[AvikPID::nbounds] = {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

//-----------------------------------------------------------------------------
// Vadim's fitting: compute sum of squares of residuals
//-----------------------------------------------------------------------------
  void AvikPID::myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //minimisation function computing the sum of squares of residuals
    Int_t np = error->GetN();
    f = 0;
    Double_t *x = error->GetX();
    Double_t *y = error->GetY();
    Double_t *ey = error->GetEY();
    //   Double_t *ey = error->GetEY();

    for (Int_t i=0;i<np;i++) {
      Double_t dr = y[i] - (par[0]+par[1]*x[i]);
      f += (dr*dr)/(ey[i]*ey[i]);
    }
  }

//-----------------------------------------------------------------------------
  int AvikPID::findlowhist(float d){

    if (d<=AvikPID::_pathbounds[0]) return 0;

    else if (d>_pathbounds[0] && d<=_pathbounds[1]) return 0;
    else if (d>_pathbounds[1] && d<=_pathbounds[2]) return 1;
    else if (d>_pathbounds[2] && d<=_pathbounds[3]) return 2;
    else if (d>_pathbounds[3] && d<=_pathbounds[4]) return 3;
    else if (d>_pathbounds[4] && d<=_pathbounds[5]) return 4;
    else if (d>_pathbounds[5] && d<=_pathbounds[6]) return 5;
    else if (d>_pathbounds[6] && d<=_pathbounds[7]) return 6;
    else if (d>_pathbounds[7] && d<=_pathbounds[8]) return 7;
    else if (d>_pathbounds[8] && d<=_pathbounds[9]) return 8;
    else if (d>_pathbounds[9] && d<=_pathbounds[10]) return 9;
    else if (d>_pathbounds[10]) return 10;

    else {
      cout<<"Out of  bounds. Should never end up here\n";
    }

    return -9999.;
  }


//-----------------------------------------------------------------------------
  AvikPID::AvikPID(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _debugLevel(pset.get<int>("debugLevel")),
    _diagLevel (pset.get<int>("diagLevel" )),

    _trkPatRecDemModuleLabel(pset.get<string>("trkPatRecDemModuleLabel")),
    _trkPatRecDmmModuleLabel(pset.get<string>("trkPatRecDmmModuleLabel")),

    _eleDedxTemplateFile(pset.get<std::string>("EleDedxTemplateFile")),
    _muoDedxTemplateFile(pset.get<std::string>("MuoDedxTemplateFile")),

    _darPset            (pset.get<fhicl::ParameterSet>("DoubletAmbigResolver")),

    _pidtree(0)
  {
    _processed_events = -1;

    _iname = _fdir.name() + _tpart.name();
    produces<AvikPIDProductCollection>();

    // location-independent files
    ConfigFileLookupPolicy configFile;

    _eleTemplates = configFile(_eleDedxTemplateFile);
    _muoTemplates = configFile(_muoDedxTemplateFile);

    char name[50];

    TFile* eleDedxTemplateFile = TFile::Open(_eleTemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempe%d",i);
      eleDedxTemplateFile->GetObject(name,_heletemp[i]);
    }

    TFile* muoDedxTemplateFile = TFile::Open(_muoTemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempm%d",i);
      muoDedxTemplateFile->GetObject(name,_hmuotemp[i]);
    }
//-----------------------------------------------------------------------------
// all electron and muon De/Dx template histograms are supposed to have the
// same limits and number of bins
//-----------------------------------------------------------------------------
    _templatesnbins   = _heletemp[0]->GetNbinsX();
    _templateslastbin = _heletemp[0]->GetBinLowEdge(_templatesnbins)+_heletemp[0]->GetBinWidth(1);
    _templatesbinsize = _heletemp[0]->GetBinWidth(1);
//-----------------------------------------------------------------------------
// Avik function parameters for dRdz slope residuals
//-----------------------------------------------------------------------------
    _pow1    = 1.5;
    _bound   = 0.0625;
    _pow2    = 0.25;

    _maxDeltaDxDzOs = 0.5;

    _minuit  = NULL;

    _dar     = new DoubletAmbigResolver(_darPset,0,0,0);

  }


//-----------------------------------------------------------------------------
  AvikPID::~AvikPID() {
    if (_minuit) delete _minuit;
    delete _dar;
  }


//-----------------------------------------------------------------------------
  void AvikPID::beginJob() {

    // histograms

    art::ServiceHandle<art::TFileService> tfs;

    if (_diagLevel) {
      _pidtree = tfs->make<TTree>("PID", "PID info");

      _pidtree->Branch("trkid"          , &_trkid            , "trkid/I");
      _pidtree->Branch("p"              , &_trkmom           , "trkmom/D");
      _pidtree->Branch("drdsVadimEle"   , &_drdsVadimEle     , "drdsVadimEle/D");
      _pidtree->Branch("drdsVadimEleErr", &_drdsVadimEleErr  , "drdsVadimEleErr/D");
      _pidtree->Branch("drdsVadimMuo"   , &_drdsVadimMuo     , "drdsVadimMuo/D");
      _pidtree->Branch("drdsVadimMuoErr", &_drdsVadimMuoErr  , "drdsVadimMuoErr/D");
      _pidtree->Branch("logDedxProbEle" , &_logDedxProbEle   , "logDedxProbEle/D");
      _pidtree->Branch("logDedxProbMuo" , &_logDedxProbMuo   , "logDedxProbMuo/D");
    }

    _minuit = new TMinuit(2);

  }

//-----------------------------------------------------------------------------
  void AvikPID::beginRun(art::Run & run){
    if (_debugLevel >= 2) cout << "AvikPID: From beginRun: " << run.id().run() << endl;


  }

//-----------------------------------------------------------------------------
  void AvikPID::beginSubRun(art::SubRun & lblock ) {
    if (_debugLevel>=2) cout << "AvikPID: From beginSubRun. " << endl;
  }

//-----------------------------------------------------------------------------
  void AvikPID::endJob(){
    if (_debugLevel>=2) cout << "AvikPID: From endJob. " << endl;
  }

//-----------------------------------------------------------------------------
  void AvikPID::doubletMaker(const KalRep* ele_Trk, const KalRep* muo_Trk) {

    art::Handle<mu2e::KalRepPtrCollection> eleKrepsHandle;
    art::Handle<mu2e::KalRepPtrCollection> muoKrepsHandle;

    TrkHitVector const& ele_hot_list = ele_Trk->hitVector();
    TrkHitVector const& muo_hot_list = muo_Trk->hitVector();
//-----------------------------------------------
// ELECTRONS:
//-----------------------------------------------
    int ele_nhits  = 0;  // total number of hits
    int ele_ncount = 0;  // number of hits in doublets/triplets etc

    for (auto it=ele_hot_list.begin(); it<ele_hot_list.end(); it++) ele_nhits+=1;

    int   ele_iamb0 = 0;

    float ele_resall   [ele_nhits];
    int   ele_resgood  [ele_nhits];
    float ele_res      [ele_nhits];
    int   ele_panelall   [ele_nhits];
    int   ele_planeall   [ele_nhits];
    int   ele_layall   [ele_nhits];
    int   ele_Nall     [ele_nhits];
    int   ele_strawall [ele_nhits];
    int   ele_straw    [ele_nhits];
    int   ele_iamball  [ele_nhits];
    int   ele_nmatchall[ele_nhits] = {0};   // number of matches for each hit

    int   ele_nnlet[6] = {ele_nhits};       // array giving number of singlets, doublets, triplets, ...

//     float h1_fltlen = ele_Trk->firstHit()->kalHit()->hit()->fltLen() - 10;
//     float hn_fltlen = ele_Trk->lastHit ()->kalHit()->hit()->fltLen() - 10;
//     float entlen = std::min(h1_fltlen,hn_fltlen);

//    BbrVectorErr momerr = ele_Trk->momentumErr(entlen);
//     CLHEP::Hep3Vector fitmom = ele_Trk->momentum(entlen);
//     CLHEP::Hep3Vector momdir = fitmom.unit();

//     HepVector momvec(3);
//     for (int i=0; i<3; i++) momvec[i] = momdir[i];

    int k = 0;
    Hep3Vector pos;
    double     hitres, hiterr;

    for (auto it=ele_hot_list.begin(); it<ele_hot_list.end(); it++) {

      const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) (*it);

      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

      ele_panelall[k]=straw->id().getPanel();
      ele_planeall[k]=straw->id().getPlane();
      ele_layall[k]=straw->id().getLayer();
      ele_Nall[k]=straw->id().getStraw();
      ele_strawall[k]= straw->id().asUint16();
      ele_iamball[k]= hit->ambig();
      hit->hitPosition(pos);
      ele_resgood[k] = hit->resid(hitres,hiterr,1);
      ele_resall[k]=hitres;

      k += 1;
    }
//-----------------------------------------------------------------------------
// loop over the electron track hits
//-----------------------------------------------------------------------------
    for(int i=0; i<ele_nhits; i++) {

      if(ele_iamball[i]==0) ele_iamb0+=1;

      for(int j=0; j<i; j++) {     //this loop checks current hit against all previous hits in list

        if((ele_panelall[i]==ele_panelall[i-(j+1)]) &&
           (ele_planeall[i]==ele_planeall[i-(j+1)]) &&
           (ele_layall[i]!=ele_layall[i-(j+1)]) &&
           (abs(ele_Nall[i]-ele_Nall[i-(j+1)])<=2) &&
           (ele_iamball[i]!=ele_iamball[i-(j+1)]) &&
           (ele_resgood[i]==1) && (ele_resgood[i-(j+1)]==1)) {   //doublet conditions

          ele_nmatchall[i]+=1;    //number of matches to other hits
          if((ele_nmatchall[i-(j+1)]==0) && (ele_nmatchall[i]==1)) {  //i.e. has found a doublet
            ele_res[ele_ncount+1]=ele_resall[i];   // i is current hit
            ele_straw[ele_ncount+1]=ele_strawall[i];
            ele_res[ele_ncount]=ele_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
            ele_straw[ele_ncount]=ele_strawall[i-(j+1)];
            ele_ncount+=2;
            ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            ele_nnlet[1]+=1;
            ele_nnlet[0]-=2;
          }

          else if (ele_nmatchall[i-(j+1)]==0) {
            ele_res[ele_ncount]=ele_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
            ele_straw[ele_ncount]=ele_strawall[i-(j+1)];
            ele_ncount+=1;
            ele_nnlet[ele_nmatchall[i]]+=1;
            ele_nnlet[ele_nmatchall[i]-1]-=1;
            ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            ele_nnlet[0]-=1;
          }

          else {
            ele_res[ele_ncount]=ele_resall[i];   // i is current hit
            ele_straw[ele_ncount]=ele_strawall[i];
            ele_ncount+=1;
            ele_nnlet[ele_nmatchall[i-(j+1)]+ele_nmatchall[i]]+=1;
            ele_nnlet[ele_nmatchall[i-(j+1)]]-=1;
            ele_nnlet[ele_nmatchall[i]-1]-=1;
            ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// MUONS:
//-----------------------------------------------------------------------------
    int muo_nhits  = 0; // total number of hits
    int muo_ncount = 0; // number of special hits (i.e. events in doublets or triplets etc)

    for(auto it=muo_hot_list.begin(); it<muo_hot_list.end(); it++) muo_nhits+=1;

    int   muo_iamb0=0;
    float muo_resall   [muo_nhits];
    float muo_res      [muo_nhits];
    int   muo_resgood  [muo_nhits];
    int   muo_panelall   [muo_nhits];
    int   muo_planeall   [muo_nhits];
    int   muo_layall   [muo_nhits];
    int   muo_Nall     [muo_nhits];
    int   muo_strawall [muo_nhits];
    int   muo_straw    [muo_nhits];
    int   muo_iamball  [muo_nhits];
    int   muo_nmatchall[muo_nhits] = {0};   //number of matches for each hit
    int   muo_nnlet    [6] = {muo_nhits};    //array giving number of singlets, doublets, triplets, ...
    //    float muo_chi2 = muo_Trk->chisq();

    k = 0;
    for(auto it=muo_hot_list.begin(); it<muo_hot_list.end(); it++) {

      const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) (*it);

      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

      muo_panelall[k]  = straw->id().getPanel();
      muo_planeall[k]  = straw->id().getPlane();
      muo_layall[k]  = straw->id().getLayer();
      muo_Nall[k]    = straw->id().getStraw();
      muo_strawall[k]= straw->id().asUint16();
      muo_iamball[k] = hit->ambig();
      hit->hitPosition(pos);
      muo_resgood[k] = hit->resid(hitres,hiterr,1);
      muo_resall [k] = hitres;
      k             += 1;
    }

    for(int i=0; i<muo_nhits; i++) {  //loop goes through each hit in the track

      if(muo_iamball[i]==0) muo_iamb0+=1;

      for(int j=0; j<i; j++){     //this loop checks current hit against all previous hits in list

        if((muo_panelall[i]==muo_panelall[i-(j+1)]) &&
           (muo_planeall[i]==muo_planeall[i-(j+1)]) &&
           (muo_layall[i]!=muo_layall[i-(j+1)]) &&
           (abs(muo_Nall[i]-muo_Nall[i-(j+1)])<=2) &&
           (muo_iamball[i]!=muo_iamball[i-(j+1)]) &&
           (muo_resgood[i]==1) && (muo_resgood[i-(j+1)]==1)) {   //doublet conditions

          muo_nmatchall[i]+=1;    //number of matches to other hits

          if((muo_nmatchall[i-(j+1)]==0) && (muo_nmatchall[i]==1)) {  //i.e. has found a doublet
            muo_res[muo_ncount+1]=muo_resall[i];   // i is current hit
            muo_straw[muo_ncount+1]=muo_strawall[i];
            muo_res[muo_ncount]=muo_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
            muo_straw[muo_ncount]=muo_strawall[i-(j+1)];
            muo_ncount+=2;
            muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            muo_nnlet[1]+=1;
            muo_nnlet[0]-=2;
          }
          else if (muo_nmatchall[i-(j+1)]==0) {
            muo_res[muo_ncount]=muo_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
            muo_straw[muo_ncount]=muo_strawall[i-(j+1)];
            muo_ncount+=1;
            muo_nnlet[muo_nmatchall[i]]+=1;
            muo_nnlet[muo_nmatchall[i]-1]-=1;
            muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            muo_nnlet[0]-=1;
          }

          else {
            muo_res[muo_ncount]=muo_resall[i];   // i is current hit
            muo_straw[muo_ncount]=muo_strawall[i];
            muo_ncount+=1;
            muo_nnlet[muo_nmatchall[i-(j+1)]+muo_nmatchall[i]]+=1;
            muo_nnlet[muo_nmatchall[i-(j+1)]]-=1;
            muo_nnlet[muo_nmatchall[i]-1]-=1;
            muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// ANALYSIS:
//-----------------------------------------------------------------------------
    float res_ele    [ele_ncount];
    float res_muo    [muo_ncount];
    float res_ele_sq [ele_ncount];
    float res_muo_sq [muo_ncount];

    float resall_ele [ele_nhits] = {0};
    float resall_muo [muo_nhits] = {0};
    float res_ele_sq2[ele_nhits] = {0};
    float res_muo_sq2[muo_nhits] = {0};

    float res_ele_sum  = 0;
    float res_muo_sum  = 0;
    float res_ele_sum2 = 0;
    float res_muo_sum2 = 0;

    for(int i=0; i<ele_ncount; i++) {
      res_ele   [i]=0.0;
      res_ele_sq[i]=0.0;
    }

    for(int i=0; i<muo_ncount; i++) {
      res_muo   [i]=0.0;
      res_muo_sq[i]=0.0;
    }

    for (int i=0; i<ele_nhits; i++) {
      resall_ele [i] = 0.;
      res_ele_sq2[i] = 0.;
    }

    for (int i=0; i<muo_nhits; i++) {
      resall_muo [i] = 0.;
      res_muo_sq2[i] = 0.;
    }

    for(int i=0; i<ele_ncount; i++) {
      for(int j=0; j<muo_ncount; j++) {
        if(ele_straw[i]==muo_straw[j]) {   // if a hit is part of doublets in both hypotheses...
          res_ele[i]=ele_res[i];           // enter residuals for both respectively
          res_muo[j]=muo_res[j];
        }
      }
    }

    for(int i=0; i<ele_nhits; i++) {
      for(int j=0; j<muo_nhits; j++) {
        if(ele_strawall[i]==muo_strawall[j]) {   // if hit exists in both hypotheses...
          resall_ele[i]=ele_resall[i];           // enter residuals for both respectively
          resall_muo[j]=muo_resall[j];
        }
      }
    }

    for (int i=0; i<ele_ncount; i++) {
      if (abs(res_ele[i])>1.5) {   //if for either hypothesis residuals are too high
        res_ele[i]=0;
        for (int j=0; j<muo_ncount; j++) {
          if (ele_straw[i]==muo_straw[j]) res_muo[j]=0;   //set them both to zero
        }
      }
    }

    for (int i=0; i<muo_ncount; i++) {
      if (abs(res_muo[i])>1.5) {   //if for either hypothesis residuals are too high
        res_muo[i]=0;
        for (int j=0; j<ele_ncount; j++) {
          if (muo_straw[i]==ele_straw[j]) res_ele[j]=0;   //set them both to zero
        }
      }
    }

    for (int i=0; i<ele_nhits; i++) {
      if (abs(resall_ele[i])>1.5) {   //if for either hypothesis residuals are too high
        resall_ele[i]=0;
        for (int j=0; j<muo_nhits; j++) {
          if (ele_strawall[i]==muo_strawall[j]) resall_muo[j]=0;   //set them both to zero
        }
      }
    }

    for (int i=0; i<muo_nhits; i++) {
      if (abs(resall_muo[i])>1.5) {   // if for either hypothesis residuals are too high
        resall_muo[i]=0;
        for (int j=0; j<ele_nhits; j++) {
          if (muo_strawall[i]==ele_strawall[j]) resall_ele[j]=0;   //set them both to zero
        }
      }
    }

    float matchhits = 0.0;
    float matchhits_all = 0.0;

    for(int i=0; i<ele_ncount; i++) {
      if (res_ele[i]!=0) matchhits+=1.0;

      res_ele_sq[i] = weightedResidual(res_ele[i]);
      res_ele_sum  += res_ele_sq[i];
    }

    for(int i=0; i<muo_ncount; i++) {
      res_muo_sq[i] = weightedResidual(res_muo[i]);
      res_muo_sum  += res_muo_sq[i];
    }

    for(int i=0; i<ele_nhits; i++) {
      if (resall_ele[i]!=0) matchhits_all+=1.0;

      res_ele_sq2[i] = weightedResidual(resall_ele[i]);
      res_ele_sum2  += res_ele_sq2[i];
    }

    for(int i=0; i<muo_nhits; i++) {
      res_muo_sq2[i] = weightedResidual(resall_muo[i]);
      res_muo_sum2  += res_muo_sq2[i];
    }

    _sumAvikEle  = res_ele_sum;
    _sumAvikMuo  = res_muo_sum;

    _sq2AvikEle  = res_ele_sum2;
    _sq2AvikMuo  = res_muo_sum2;

    _nMatched    = matchhits;
    _nMatchedAll = matchhits_all;

    float ratio, ratio2;
    float logratio, logratio2, logratio3;

    if ((_sumAvikEle != 0) && (_sumAvikMuo != 0)) ratio = _sumAvikMuo/_sumAvikEle;
    else                                          ratio = 1.;

    if ((_sq2AvikEle != 0) && (_sq2AvikMuo != 0)) ratio2 = _sq2AvikMuo/_sq2AvikEle;
    else                                          ratio2 = 1.;

    logratio  = log(ratio);
    logratio2 = log(ratio2);
    logratio3 = (matchhits/matchhits_all)*log(ratio) + log(ratio2);

    /*
       printf("ANALYSIS \n");

       for(int i=0; i<ele_ncount; i++) {
       printf( "res_ele is:  %4.4f %s %4.4f \n", res_ele[i], "  res_ele_sq is: ", res_ele_sq[i]);
       }

       for(int i=0; i<muo_ncount; i++) {
       printf( "res_muo is:  %4.4f %s %4.4f \n", res_muo[i], "  res_muo_sq is: ", res_muo_sq[i]);
       }
    */

    if (_debugLevel > 0) {
      printf( "res_ele_sum  is:  %8.4f res_muo_sum  is: %8.4f res_ele_sum2 is:  %8.4f res_muo_sum2 is: %8.4f\n",
	      res_ele_sum ,res_muo_sum, res_ele_sum2,res_muo_sum2);
      printf("logratio  is: %8.4f logratio2 is: %8.4f logratio3 is: %8.4f\n",
	     logratio, logratio2, logratio3);
    }
  }

//-----------------------------------------------------------------------------
  double AvikPID::calculateDedxProb(std::vector<double>* GasPaths ,
                                    std::vector<double>* EDeps    ,
                                    TH1D**               Templates) {

    static const double _minpath = 0.5;
    static const double _maxpath = 10.;

    double thisprob = 1;

    for (unsigned int ipath = 0; ipath < GasPaths->size(); ipath++){
      double thispath = GasPaths->at(ipath);
      double thisedep = EDeps->at(ipath);

      double tmpprob = 0;
      if (thispath > _minpath && thispath<=_maxpath){
          int lowhist = findlowhist(thispath);

          PIDUtilities util;
          TH1D* hinterp = util.th1dmorph(Templates[lowhist],
                                         Templates[lowhist+1],
                                         _pathbounds[lowhist],
                                         _pathbounds[lowhist+1],
                                         thispath,1,0);
//-----------------------------------------------------------------------------
// probability for this edep
//-----------------------------------------------------------------------------
          int thisedepbin = -999;
          if (thisedep > _templateslastbin) thisedepbin = _templatesnbins;
          else                              thisedepbin = int(thisedep/_templatesbinsize)+1;

          tmpprob=hinterp->GetBinContent(thisedepbin);

          hinterp->Delete();
        }

      if (tmpprob> 0) thisprob = thisprob * tmpprob;
    }

    return thisprob;
  }

//-----------------------------------------------------------------------------
// the straight line fit should be done explicitly
//-----------------------------------------------------------------------------
  bool AvikPID::calculateVadimSlope(const KalRep  *KRep  ,
                                    double        *Slope ,
                                    double        *Eslope) {

    std::vector<double> res, flt, eres, eflt;
    mu2e::TrkStrawHit*  hit;
    double              resid, residerr, aresd, normflt, normresd;
    TrkHitVector const&   hotList = KRep->hitVector();

    for (auto ihot=hotList.begin(); ihot != hotList.end(); ++ihot) {
      hit = (mu2e::TrkStrawHit*) (*ihot);
      if (hit->isActive()) {
//-----------------------------------------------------------------------------
// 'unbiased' residual, signed with the radius - if the drift radius is greater than
// the track-to-wire distance, the residual is positive (need to double check the sign!)
// use active hits only
//-----------------------------------------------------------------------------
        hit->resid(resid,residerr,true);

        aresd    = (hit->poca().doca()>0?resid:-resid);
        normflt  = hit->fltLen() -  KRep->flt0();
        normresd = aresd/residerr;

        res.push_back(normresd);
        flt.push_back(normflt);
        eres.push_back(1.);
        eflt.push_back(0.1);
      }
    }

    error = new TGraphErrors(res.size(),flt.data(),res.data(),eflt.data(),eres.data());

    _minuit->SetPrintLevel(-1);
    _minuit->SetFCN(myfcn);

    const int dim(2);
    const char      par_name[dim][20]= {"offset","slope"};
    static Double_t step    [dim]    = {0.001,0.001};
    Double_t        sfpar   [dim]    = {0.0,0.005};
    Double_t        errsfpar[dim]    = {0.0,0.0};

    int ierflg = 0;
    for (int ii = 0; ii<dim; ii++) {
      _minuit->mnparm(ii,par_name[ii],sfpar[ii], step[ii], 0,0,ierflg);
    }

    _minuit->FixParameter(0);
    _minuit->Migrad();

    bool converged = _minuit->fCstatu.Contains("CONVERGED");
    if (!converged)
      {
        cout <<"-----------TOF Linear fit did not converge---------------------------" <<endl;
        return converged;
      }
    for (int i = 0;i<dim;i++) {
      _minuit->GetParameter(i,sfpar[i],errsfpar[i]);
    }

    *Slope  = sfpar[1];
    *Eslope = errsfpar[1];

    delete error;

    return converged;
  }

//-----------------------------------------------------------------------------
  int AvikPID::AddHits(const Doublet* Multiplet, vector<double>& Fltlen, vector<double>& Resid) {
    //    int nhits = Multiplet->fNStrawHits;
    int ihit;
    double res, flt, reserr;

    for (int i=0; i<2; i++) {
      ihit = Multiplet->fHitIndex[i];
      flt  = Multiplet->fHit[ihit]->fltLen();
      Multiplet->fHit[ihit]->resid(res, reserr, true);
      res  = (Multiplet->fHit[ihit]->poca().doca()>0?res:-res);
      Fltlen.push_back(flt);
      Resid.push_back(res);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// doublet ambiguity resolver best combinations: 0:(++) 1:(+-) 2:(--) 3:(-+)
// so 0 and 2 correspond to the SS doublet, 1 and 3 - to the OS doublet
// see KalmanTests/src/DoubletAmbigResolver.cc for details
// use SS doublets, require the best slope to be close to that of the track
//-----------------------------------------------------------------------------
  int AvikPID::AddSsMultiplets(const vector<Doublet>* ListOfDoublets,
                               vector<double>&        Fltlen        ,
                               vector<double>&        Resid         ) {

    const mu2e::Doublet  *multiplet;
    double               dxdzresid;

    int ndblts = ListOfDoublets->size();

    for (int i=0; i<ndblts; i++) {
      multiplet = &ListOfDoublets->at(i);
      int nhits = multiplet->fNStrawHits;
      if ((nhits > 1) && multiplet->isSameSign()) {
        dxdzresid = multiplet->bestDxDzRes();
//-----------------------------------------------------------------------------
// always require the local doublet slope to be close to that of the track
//-----------------------------------------------------------------------------
        if (fabs(dxdzresid) < .1) {
            AddHits(multiplet,Fltlen,Resid);
        }
      }
    }

    return 0;
  }


//-----------------------------------------------------------------------------
// doublet ambiguity resolver best combinations: 0:(++) 1:(+-) 2:(--) 3:(-+)
// so 0 and 2 correspond to the SS doublet, 1 and 3 - to the OS doublet
// see KalmanTests/src/DoubletAmbigResolver.cc for details
//-----------------------------------------------------------------------------
  int AvikPID::AddOsMultiplets(const vector<Doublet>* ListOfDoublets,
                               vector<double>&        Fltlen        ,
                               vector<double>&        Resid         ) {

    const mu2e::Doublet  *multiplet, *mj;
    int                  best, bestj, sid, sidj;
    double               trkdxdz, bestdxdz, dxdzresid, trkdxdzj, bestdxdzj, dxdzresidj;

    int      ndblts = ListOfDoublets->size();

    for (int i=0; i<ndblts; i++) {
      multiplet = &ListOfDoublets->at(i);
      sid       = multiplet->fStationId/2;
      int nhits = multiplet->fNStrawHits;
      if (nhits > 1) {
        best      = multiplet->fIBest;
        trkdxdz   = multiplet->fTrkDxDz;
        bestdxdz  = multiplet->fDxDz[best];
        dxdzresid = fabs(trkdxdz - bestdxdz);
//-----------------------------------------------------------------------------
// always require the local doublet slope to be close to that of the track
//-----------------------------------------------------------------------------
        if (dxdzresid < .1) {
          if ((best == 1) || (best == 3)) {
//-----------------------------------------------------------------------------
// OS doublet
//-----------------------------------------------------------------------------
            AddHits(multiplet,Fltlen,Resid);
          }
          else {
//-----------------------------------------------------------------------------
// SS multiplet
// check if there is a OS doublet in the same station and add this SS doublet
// only if the OS one exists - in hope that the OS one would keep coordinates
// in place
//-----------------------------------------------------------------------------
            for (int j=0; j<ndblts; j++) {
              mj   = &ListOfDoublets->at(j);
              sidj = mj->fStationId/2;
              if (sid == sidj) {
                int nhj = mj->fNStrawHits;
                if (nhj > 1) {
//-----------------------------------------------------------------------------
// don't use single hit multiplets. However, in principle they could be used
//-----------------------------------------------------------------------------
                  bestj      = mj->fIBest;
                  trkdxdzj   = mj->fTrkDxDz;
                  bestdxdzj  = mj->fDxDz[bestj];
                  dxdzresidj = fabs(trkdxdzj - bestdxdzj);
//-----------------------------------------------------------------------------
// always require the local doublet slope to be close to that of the track
//-----------------------------------------------------------------------------
                  if ((dxdzresidj < .1) && ((bestj == 1) || (bestj == 3))) {
                    // best multiplet is SS
                    AddHits(multiplet,Fltlen,Resid);
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }

    return 0;
  }


//-----------------------------------------------------------------------------
// calculate parameters of the straight line fit
//-----------------------------------------------------------------------------
  int AvikPID::CalculateSlope(vector<double>& Fltlen, vector<double>& Res, double& Slope, double& Err) {

    double dl, dr, fltSum(0), resSum(0), fltMean, resMean, fltVar(0), resVar(0), fltRes(0), fltDev, resDev, rCoeff;

    int n = Fltlen.size();

    if (n > 1) {
      for (int i=0; i<n; i++) {
        fltSum += Fltlen[i];
        resSum += Res[i];
      }

      fltMean = fltSum/n;
      resMean = resSum/n;

      for (int i=0; i<n; i++) {
        dl      = Fltlen[i] - fltMean;
        dr      = Res   [i] - resMean;
        fltVar += dl*dl;
        resVar += dr*dr;
        fltRes += dl*dr;
      }

      fltDev = sqrt(fltVar/n);              // sigxx
      resDev = sqrt(resVar/n);              // sigyy
      rCoeff = fltRes/sqrt(fltVar*resVar);  //

      Slope  = rCoeff*resDev/fltDev;
      Err    = 0.1/sqrt(fltVar);
    }
    else {
//-----------------------------------------------------------------------------
// degenerate case
//-----------------------------------------------------------------------------
      Slope = 1.e6;
      Err   = 1.e6;
    }
    return 0;
  }

//--------------------------------------------------------------------------------
  void AvikPID::calculateSsSums(const vector<Doublet>* ListOfDoublets, double&  Drds, double& DrdsErr, int& NUsed) {

    vector<double> fltLen, resid;

    AddSsMultiplets(ListOfDoublets,fltLen,resid);

    NUsed = fltLen.size();

    CalculateSlope(fltLen,resid,Drds,DrdsErr);
  }

//------------------------------------------------------------------------------
// calculate sums over the local doublet residuals
// residuals are weighted, as Avik is trying to de-weight the tails
//-----------------------------------------------------------------------------
  void AvikPID::calculateOsSums(const vector<Doublet>* ListOfDoublets,
                                double& Drds, double& DrdsErr, int& NUsedHits,
                                double& Sum , int& NUsedDoublets) {

    vector<double> fltLen;
    vector<double> resid;
    const Doublet* multiplet;
//-----------------------------------------------------------------------------
// calculate slopes using OS doublets
//-----------------------------------------------------------------------------
    AddOsMultiplets(ListOfDoublets,fltLen,resid);
    CalculateSlope(fltLen,resid,Drds,DrdsErr);
    NUsedHits = fltLen.size();


    NUsedDoublets = 0;
    Sum           = 0.;

    int nnlets    = ListOfDoublets->size();
    int ndblts    = 0.;

    for (int i=0; i<nnlets; i++) {
      multiplet    = &ListOfDoublets->at(i);
      int nhits    = multiplet->fNStrawHits;
      if (nhits >= 2) {
        ndblts += 1;
        if (! multiplet->isSameSign()) {
          double ddxdz    = multiplet->bestDxDzRes();
          if (fabs(ddxdz) < _maxDeltaDxDzOs) {
//          double residsq = 0;
//          if (abs(dxdzresid) < _bound) residsq = (1/pow(_bound, _pow1))*pow(abs(dxdzresid), _pow1);
//          else residsq   = (1/pow(_bound, _pow1))*pow(_bound, (_pow1 - _pow2))*pow(abs(dxdzresid), _pow2);
            double dr2     = weightedSlopeResidual(ddxdz);
            Sum           += dr2;
            NUsedDoublets += 1.;
          }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
  void AvikPID::produce(art::Event& event) {

    auto detmodel = _mu2eDetector_h.getPtr(event.id());

    art::Handle<mu2e::KalRepPtrCollection> eleHandle, muoHandle;

    double         firsthitfltlen, lasthitfltlen, entlen;
    double         path, eprob, muprob;

    int const      max_ntrk(100);
    int            n_ele_trk, n_muo_trk, ele_unique[max_ntrk], muo_unique[max_ntrk];
    int            found, ele_nhits, muo_nhits, ncommon;

    const mu2e::ComboHit               *esh, *msh;
    //    mu2e::DoubletAmbigResolver::Data_t r;

    vector<double>         gaspaths;
    vector<double>         edeps;

    vector<Doublet>        ele_listOfDoublets, muo_listOfDoublets;

    const KalRep            *ele_Trk, *muo_Trk;

    TrkHitVector        ele_hots, muo_hots;
    const mu2e::TrkStrawHit *hit, *ehit, *mhit;

    _evtid = event.id().event();

    ++_processed_events;

    if (_processed_events%100 == 0) {
      if (_debugLevel >= 1) cout << "AvikPID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;
    }

    if (_debugLevel >= 2) cout << "AvikPID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;

    unique_ptr<AvikPIDProductCollection> pids(new AvikPIDProductCollection );

    art::Selector  ele_selector(art::ProcessNameSelector("")         &&
                                art::ModuleLabelSelector(_trkPatRecDemModuleLabel));

    art::Selector  muo_selector(art::ProcessNameSelector("")         &&
                                art::ModuleLabelSelector(_trkPatRecDmmModuleLabel));

    event.get(ele_selector,eleHandle);
    event.get(muo_selector,muoHandle);

    if (! eleHandle.isValid()) {
      printf("TAnaDump::printKalRepPtrCollection: no ELE KalRepPtrCollection for module, BAIL OUT\n");
      goto END;
    }

    if (! muoHandle.isValid()) {
      printf("TAnaDump::printKalRepPtrCollection: no MUO KalRepPtrCollection for module, BAIL OUT\n");
      goto END;
    }

    _listOfEleTracks = eleHandle.product();
    _listOfMuoTracks = muoHandle.product();

    n_ele_trk = _listOfEleTracks->size();
    n_muo_trk = _listOfMuoTracks->size();

    if (_debugLevel > 0) {
      printf("Event: %8i : n_ele_trk: %2i n_muo_trk: %2i\n",_evtid,n_ele_trk,n_muo_trk);
    }
//-----------------------------------------------------------------------------
// proceed further
//-----------------------------------------------------------------------------
    for (int i=0; i<max_ntrk; i++) {
      ele_unique[i] = 1;
      muo_unique[i] = 1;
    }

    if ((n_ele_trk > max_ntrk/2) || (n_muo_trk > max_ntrk/2)) {
      printf("Event: %8i : n_ele_trk: %2i n_muo_trk: %2i. BAIL OUT\n",
             _evtid,n_ele_trk,n_muo_trk);
      goto END;
    }

    for (int i=0; i<n_ele_trk; i++) {
      _ele_trkid     = i;
      ele_Trk        = _listOfEleTracks->at(i).get();
      ele_hots       = ele_Trk->hitVector();
      ele_nhits      = ele_Trk->nActive();
//-----------------------------------------------------------------------------
// electron track hit doublets
//-----------------------------------------------------------------------------
      _dar->findDoublets  (ele_Trk,&ele_listOfDoublets);
//-----------------------------------------------------------------------------
// dE/dX for the electron track hits
// calculate dE/dX, clear vectors, start forming a list of hits from the electron track
//-----------------------------------------------------------------------------
      firsthitfltlen = ele_Trk->firstHit()->kalHit()->hit()->fltLen() - 10;
      lasthitfltlen  = ele_Trk->lastHit()->kalHit()->hit()->fltLen() - 10;
      entlen         = std::min(firsthitfltlen,lasthitfltlen);
      _trkmom        = ele_Trk->momentum(entlen).mag();

      gaspaths.clear();
      edeps.clear();

      for (auto ihot=ele_hots.begin(); ihot != ele_hots.end(); ++ihot) {
        const TrkStrawHit* hit = dynamic_cast<const mu2e::TrkStrawHit*> (*ihot);
        if (hit && hit->isActive()) {
//-----------------------------------------------------------------------------
// hit charges: '2.*' here because KalmanFit reports half-path through gas.
//-----------------------------------------------------------------------------
          const DetStrawElem* strawelem = detmodel->strawElem(hit->straw());
          path = 2.*strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction(hit->fltLen()));
          gaspaths.push_back(path);
          edeps.push_back(hit->comboHit().energyDep());
        }
      }
//-----------------------------------------------------------------------------
// calculate ddR/ds slope for the electron tracks
//-----------------------------------------------------------------------------
      calculateVadimSlope(ele_Trk,&_drdsVadimEle,&_drdsVadimEleErr);

      calculateOsSums(&ele_listOfDoublets,_drdsOsEle,_drdsOsEleErr,_ele_nusedOsH,_ele_resSumOs,_ele_nusedOsD);
      calculateSsSums(&ele_listOfDoublets,_drdsSsEle,_drdsSsEleErr,_ele_nusedSsH);
//-----------------------------------------------------------------------------
// loop over muon tracks
//-----------------------------------------------------------------------------
      for (int j=0; j<n_muo_trk; j++) {
        _muo_trkid     = j;
        muo_Trk        = _listOfMuoTracks->at(j).get();
        muo_hots       = muo_Trk->hitVector();
        muo_nhits      = muo_Trk->nActive();
//-----------------------------------------------------------------------------
// check if muon and electron tracks are the same track - use 50% of active hits
// of the same hits criterion
//-----------------------------------------------------------------------------
        ncommon = 0;
        for(auto ite=ele_hots.begin(); ite<ele_hots.end(); ite++) {
          ehit = dynamic_cast<const mu2e::TrkStrawHit*> (*ite);
          if (ehit && ehit->isActive()) {
            for(auto itm=muo_hots.begin(); itm<muo_hots.end(); itm++) {
              mhit = dynamic_cast<const mu2e::TrkStrawHit*> (*itm);
              if (mhit && mhit->isActive()) {
                if (&ehit->comboHit() == &mhit->comboHit()) {
                  ncommon += 1;
                  break;
                }
              }
            }
          }
        }

        if (ncommon > (ele_nhits+muo_nhits)/4.) {
//-----------------------------------------------------------------------------
// consider the two tracks to be the same
// 2. look at the muon track hits and add active muon track hits not present
//    in the electron track list
//-----------------------------------------------------------------------------
          ele_unique[i] = 0;
          muo_unique[j] = 0;
//-----------------------------------------------------------------------------
// list of muon doublets
//-----------------------------------------------------------------------------
          _dar->findDoublets  (muo_Trk,&muo_listOfDoublets);

          for (auto ihot=muo_hots.begin(); ihot != muo_hots.end(); ++ihot) {
            const TrkStrawHit* hit = dynamic_cast<const mu2e::TrkStrawHit*> (*ihot);
            if (hit && hit->isActive()) {
              msh = &hit->comboHit();
//-----------------------------------------------------------------------------
// check if 'hit' is unique for the muon track
//-----------------------------------------------------------------------------
              found = 0;
              for (auto ehot=ele_hots.begin(); ehot != ele_hots.end(); ++ehot) {
                ehit = dynamic_cast<const mu2e::TrkStrawHit*> (*ehot);
                if (ehit && ehit->isActive()) {
                  esh  = &ehit->comboHit();
                  if (esh == msh) {
                    found = 1;
                    break;
                  }
                }
              }

              if (found == 0) {
//-----------------------------------------------------------------------------
// straw hit present in the list of active muon track hits, but not in the list
// of active electron track hits, add it to the list of hits used in de/dx calculation
//-----------------------------------------------------------------------------
                const DetStrawElem* strawelem = detmodel->strawElem(hit->straw());
                path = 2.*strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction(hit->fltLen()));
                gaspaths.push_back(path);
                edeps.push_back(hit->comboHit().energyDep());
              }
            }
          }

          eprob  = calculateDedxProb(&gaspaths, &edeps, _heletemp);
          muprob = calculateDedxProb(&gaspaths, &edeps, _hmuotemp);

          _logDedxProbEle = log(eprob );
          _logDedxProbMuo = log(muprob);
//-----------------------------------------------------------------------------
// calculate Vadim's ddR/ds slopes and SS and OS ddR/ds slopes for the muon tracks
// also: OS sums of slope residuals
//-----------------------------------------------------------------------------
          calculateVadimSlope(muo_Trk,&_drdsVadimMuo,&_drdsVadimMuoErr);
          calculateSsSums(&muo_listOfDoublets,_drdsSsMuo,_drdsSsMuoErr,_muo_nusedSsH);
          calculateOsSums(&muo_listOfDoublets,_drdsOsMuo,_drdsOsMuoErr,_muo_nusedOsH,_muo_resSumOs,_muo_nusedOsD);
//-----------------------------------------------------------------------------
// calculate Avik's sums
//-----------------------------------------------------------------------------
          doubletMaker(ele_Trk,muo_Trk);
//-----------------------------------------------------------------------------
// calculate OS slopes - these make sense only when both - electron and muon -
// versions of he track are present
//-----------------------------------------------------------------------------
          _pid.init(_ele_trkid     , _muo_trkid       ,
                    _logDedxProbEle, _logDedxProbMuo  ,
                    _drdsVadimEle  , _drdsVadimEleErr ,
                    _drdsVadimMuo  , _drdsVadimMuoErr ,
                    _nMatched      , _nMatchedAll     ,
                    _sumAvikEle    , _sumAvikMuo      ,
                    _sq2AvikEle    , _sq2AvikMuo      ,
                    _drdsOsEle     , _drdsOsEleErr    ,
                    _drdsOsMuo     , _drdsOsMuoErr    ,
                    _ele_nusedOsH  , _muo_nusedOsH    ,
                    _drdsSsEle     , _drdsSsEleErr    ,
                    _drdsSsMuo     , _drdsSsMuoErr    ,
                    _ele_nusedSsH  , _muo_nusedSsH    ,
                    _ele_resSumOs  , _muo_resSumOs    ,
                    _ele_nusedOsD  , _muo_nusedOsD    );

          pids->push_back(_pid);
          break;
        }
        else {
//-----------------------------------------------------------------------------
// electron and muon tracks are different, proceed with the next muon track
//-----------------------------------------------------------------------------
        }
      }
//-----------------------------------------------------------------------------
// if electron track is unique, still can use the energy depositions to calculate
// dE/dx probability under muon hypothesis
//-----------------------------------------------------------------------------
      if (ele_unique[i] == 1) {
        _muo_trkid       = -1;
        _drdsVadimMuo    =  1.e6;
        _drdsVadimMuoErr =  1.e6;
        _sumAvikMuo      = -1;
        _sq2AvikMuo      = -1;
        _drdsOsMuo       =  1.e6;
        _drdsOsMuoErr    =  1.e6;
        _muo_nusedSsH    = -1;
        _drdsSsMuo       =  1.e6;
        _drdsSsMuoErr    =  1.e6;
        _muo_nusedOsH    = -1;
        _muo_resSumOs    =  1.e6;
        _muo_nusedOsD    = -1;
        _nMatched        = -1;
        _nMatchedAll     = -1;

        eprob  = calculateDedxProb(&gaspaths, &edeps, _heletemp);
        muprob = calculateDedxProb(&gaspaths, &edeps, _hmuotemp);

        _logDedxProbEle = log(eprob );
        _logDedxProbMuo = log(muprob);

        _pid.init(_ele_trkid     , _muo_trkid       ,
                  _logDedxProbEle, _logDedxProbMuo  ,
                  _drdsVadimEle  , _drdsVadimEleErr ,
                  _drdsVadimMuo  , _drdsVadimMuoErr ,
                  _nMatched      , _nMatchedAll     ,
                  _sumAvikEle    , _sumAvikMuo      ,
                  _sq2AvikEle    , _sq2AvikMuo      ,
                  _drdsOsEle     , _drdsOsEleErr    ,
                  _drdsOsMuo     , _drdsOsMuoErr    ,
                  _ele_nusedSsH  , _muo_nusedSsH    ,
                  _drdsSsEle     , _drdsSsEleErr    ,
                  _drdsSsMuo     , _drdsSsMuoErr    ,
                  _ele_nusedOsH  , _muo_nusedOsH    ,
                  _ele_resSumOs  , _muo_resSumOs    ,
                  _ele_nusedOsD  , _muo_nusedOsD    );

        pids->push_back(_pid);
      }
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
      if (_diagLevel) {
        _pidtree->Fill();
      }
    }

//-----------------------------------------------------------------------------
// finally: loop over muon tracks one more time and store PID records for unique ones
//-----------------------------------------------------------------------------
    for (int j=0; j<n_muo_trk; j++) {
      if (muo_unique[j] == 1) {
        _muo_trkid     = j;
        muo_Trk        = _listOfMuoTracks->at(j).get();
        muo_hots       = muo_Trk->hitVector();
        muo_nhits      = muo_Trk->nActive();
//-----------------------------------------------------------------------------
// lists of doublets need to be recreated even when they existed - their
// presence depends on ambiguity resolvers used by the reconstruction modules
//-----------------------------------------------------------------------------
        _dar->findDoublets(muo_Trk,&muo_listOfDoublets);
//-----------------------------------------------------------------------------
// dE/dX for the muoon track hits
// calculate dE/dX, clear vectors, start forming a list of hits from the electron track
//-----------------------------------------------------------------------------
        firsthitfltlen = muo_Trk->firstHit()->kalHit()->hit()->fltLen() - 10;
        lasthitfltlen  = muo_Trk->lastHit()->kalHit()->hit()->fltLen() - 10;
        entlen         = std::min(firsthitfltlen,lasthitfltlen);
        _trkmom        = muo_Trk->momentum(entlen).mag();

        gaspaths.clear();
        edeps.clear();

        for (auto ihot=muo_hots.begin(); ihot != muo_hots.end(); ++ihot) {
          hit = dynamic_cast<const mu2e::TrkStrawHit*> (*ihot);
          if (hit && hit->isActive()) {
//-----------------------------------------------------------------------------
// hit charges: '2.*' here because KalmanFit reports half-path through gas.
//-----------------------------------------------------------------------------

            const DetStrawElem* strawelem = detmodel->strawElem(hit->straw());
            path = 2.*strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction(hit->fltLen()));
            gaspaths.push_back(path);
            edeps.push_back(hit->comboHit().energyDep());
          }
        }

        eprob  = calculateDedxProb(&gaspaths, &edeps, _heletemp);
        muprob = calculateDedxProb(&gaspaths, &edeps, _hmuotemp);

        _logDedxProbEle = log(eprob );
        _logDedxProbMuo = log(muprob);

//-----------------------------------------------------------------------------
// calculate Vadim's ddR/ds slopes and SS slopes for the muon track
//-----------------------------------------------------------------------------
        calculateVadimSlope(muo_Trk,&_drdsVadimMuo,&_drdsVadimMuoErr);
        calculateSsSums(&muo_listOfDoublets,_drdsSsMuo,_drdsSsMuoErr,_muo_nusedSsH);
        calculateOsSums(&muo_listOfDoublets,_drdsOsMuo,_drdsOsMuoErr,_muo_nusedOsH,_muo_resSumOs,_muo_nusedOsD);

        _ele_trkid       = -1;
        _drdsVadimEle    =  1.e6;
        _drdsVadimEleErr =  1.e6;
        _sumAvikEle      = -1;
        _sq2AvikEle      = -1;
        _drdsOsEle       =  1.e6;
        _drdsOsEleErr    =  1.e6;
        _ele_nusedSsH    = -1;
        _drdsSsEle       =  1.e6;
        _drdsSsEleErr    =  1.e6;
        _ele_nusedOsH    = -1;
        _ele_resSumOs    = -1;
        _ele_nusedOsD    = -1;
        _nMatched        = -1;
        _nMatchedAll     = -1;

        _pid.init(_ele_trkid     , _muo_trkid       ,
                  _logDedxProbEle, _logDedxProbMuo  ,
                  _drdsVadimEle  , _drdsVadimEleErr ,
                  _drdsVadimMuo  , _drdsVadimMuoErr ,
                  _nMatched      , _nMatchedAll     ,
                  _sumAvikEle    , _sumAvikMuo      ,
                  _sq2AvikEle    , _sq2AvikMuo      ,
                  _drdsOsEle     , _drdsOsEleErr    ,
                  _drdsOsMuo     , _drdsOsMuoErr    ,
                  _ele_nusedSsH  , _muo_nusedSsH    ,
                  _drdsSsEle     , _drdsSsEleErr    ,
                  _drdsSsMuo     , _drdsSsMuoErr    ,
                  _ele_nusedOsH  , _muo_nusedOsH    ,
                  _ele_resSumOs  , _muo_resSumOs    ,
                  _ele_nusedOsD  , _muo_nusedOsD    );

        pids->push_back(_pid);
      }
    }
//-----------------------------------------------------------------------------
// end of the routine
//-----------------------------------------------------------------------------
  END:;
    event.put(std::move(pids));
  }

//-----------------------------------------------------------------------------
// Avik's weighted residual
//-----------------------------------------------------------------------------
  double AvikPID::weightedResidual(double R) {
    double wr(1.e10);

    double ar = fabs(R);
    if (ar < 0.2) wr = pow(ar, 1.7);
    else          wr = pow(0.2,1.3)*pow(ar,0.4);

    return wr;
  }

//-----------------------------------------------------------------------------
  double AvikPID::weightedSlopeResidual(double DDrDz) {

    double res(0), ar;

    ar = fabs(DDrDz);

    if (ar < _bound) res = (1/pow(_bound, _pow1))*pow(ar, _pow1);
    else             res = (1/pow(_bound, _pow1))*pow(_bound, (_pow1 -_pow2))*pow(ar, _pow2);

    return res;
  }


} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::AvikPID);
