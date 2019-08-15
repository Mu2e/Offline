// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////////
// 2016-06-02 P.Murat 
// given a list of tracks reconstructed under a certain hypothesis, calculates 
// PID variables for this track
// As the number of hypotheses can be large, to determine the best one, 
// the calculated probabilities need to be combined into likelihoods 
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
#include "RecoDataProducts/inc/AvikPIDNewProductCollection.hh"

#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"

#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"

using CLHEP::Hep3Vector;

namespace mu2e {

  class AvikPIDNew : public art::EDProducer {

  private:

    enum { kNbounds = 11 };
    static float _pathbounds[kNbounds];

    AvikPIDNewProduct _pid;

    int    _debugLevel;
    int    _verbosity;
    int    _diagLevel;

    int    _processed_events;
    int    _evtid;

    string _trkRecModuleLabel;

    string _eleDedxTemplates;
    string _eleDedxTemplateFile;
    string _muoDedxTemplates;
    string _muoDedxTemplateFile;

    ProditionsHandle<Mu2eDetector> _mu2eDetector_h;

    TrkParticle     _tpart;
    TrkFitDirection _fdir;
    std::string     _iname;		// data instance name

					// electron and muon dE/dX template histograms
    TH1D* _heletemp[kNbounds];
    TH1D* _hmuotemp[kNbounds];

    int   _templatesnbins ;
    float _templateslastbin ;
    float _templatesbinsize ;

    const KalRepPtrCollection* _listOfTracks;

    int    _trkid;
    int    _nMatched;
    int    _nMatchedAll;
    int    _nusedSsH;	       // Nhits used to calculate the SS slopes
    int    _nusedOsH;	       // Nhits used to calculate the OS slopes
    int    _nusedOsD;	       // Ndoublets used to calculate the sums

    double _trkmom;

    double _logDedxProbEle;
    double _logDedxProbMuo;
//-----------------------------------------------------------------------------
// Vadim's slopes : all hits, SS doublets, OS doublets
//-----------------------------------------------------------------------------
    double   _drdsVadim;
    double   _drdsVadimErr;
    double   _drdsSs;
    double   _drdsSsErr;
    double   _drdsOs;
    double   _drdsOsErr;
//-----------------------------------------------------------------------------
// Avik's sums
//-----------------------------------------------------------------------------
    double   _sumAvik;
    double   _sq2Avik;
    double   _resSumOs;             // d(dxdz)^alpha sums
//-----------------------------------------------------------------------------
// power coefficients
//-----------------------------------------------------------------------------
    double                 _pow1;
    double                 _pow2;
    double                 _bound;
    double                 _maxDeltaDxDzOs;

    TMinuit*               _minuit;

    fhicl::ParameterSet    _darPset;         // parameter set for doublet ambig resolver
    DoubletAmbigResolver*  _dar;

    TTree*                 _pidtree;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit AvikPIDNew(fhicl::ParameterSet const& pset);
    virtual ~AvikPIDNew();

    virtual void beginJob   ();
    virtual void beginRun   (art::Run &run);
    virtual void beginSubRun(art::SubRun & lblock );
    virtual void produce    (art::Event& event);
    virtual void endJob     ();


    static  void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);
    static  int  findlowhist(float d);

    bool   calculateVadimSlope(const KalRep* KRep, double *Slope, double *Eslope);

    double calculateDedxProb  (std::vector<double>* GasPaths ,
			       std::vector<double>* EDeps    ,
			       TH1D**               Templates);

    void   doubletMaker(const KalRep* Trk);

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

  float AvikPIDNew::_pathbounds[AvikPIDNew::kNbounds] = {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

//-----------------------------------------------------------------------------
// Vadim's fitting: compute sum of squares of residuals
//-----------------------------------------------------------------------------
  void AvikPIDNew::myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
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
  int AvikPIDNew::findlowhist(float d) {

    if (d<=AvikPIDNew::_pathbounds[0]) return 0;

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
      printf ("AvikPIDNew::findlowhist ERROR: d=%12.5e is out of  bounds. BAIL OUT\n",d);
    }

    return -9999;
  }


//-----------------------------------------------------------------------------
  AvikPIDNew::AvikPIDNew(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _debugLevel             (pset.get<int>                ("debugLevel"          )),
    _diagLevel              (pset.get<int>                ("diagLevel"           )),
    _trkRecModuleLabel      (pset.get<string>             ("trkRecModuleLabel"   )),
    _eleDedxTemplateFile    (pset.get<std::string>        ("eleDedxTemplateFile" )),
    _muoDedxTemplateFile    (pset.get<std::string>        ("muoDedxTemplateFile" )),
    _darPset                (pset.get<fhicl::ParameterSet>("DoubletAmbigResolver")),
    _pidtree                (0)
  {
    _processed_events = -1;

    _iname = _fdir.name() + _tpart.name();
    produces<AvikPIDNewProductCollection>();

    // location-independent files
    ConfigFileLookupPolicy configFile;

    _eleDedxTemplates = configFile(_eleDedxTemplateFile);
    _muoDedxTemplates = configFile(_muoDedxTemplateFile);

    TFile* eleDedxTemplateFile = TFile::Open(_eleDedxTemplates.c_str());
    TFile* muoDedxTemplateFile = TFile::Open(_muoDedxTemplates.c_str());

    char name[50];

    for (int i = 0; i < kNbounds; i++){
      sprintf(name,"htempe%d",i);
      eleDedxTemplateFile->GetObject(name,_heletemp[i]);
      sprintf(name,"htempm%d",i);
      muoDedxTemplateFile->GetObject(name,_hmuotemp[i]);
    }
//-----------------------------------------------------------------------------
// all dE/dX template histograms are supposed to have the same limits and Nbins
//-----------------------------------------------------------------------------
    _templatesnbins   = _heletemp[0]->GetNbinsX();
    _templateslastbin = _heletemp[0]->GetBinLowEdge(_templatesnbins)+_heletemp[0]->GetBinWidth(1);
    _templatesbinsize = _heletemp[0]->GetBinWidth(1);
//-----------------------------------------------------------------------------
// Avik function parameters for dRdz slope residuals
//-----------------------------------------------------------------------------
    _pow1           = 1.5;
    _bound          = 0.0625;
    _pow2           = 0.25;

    _maxDeltaDxDzOs = 0.5;

    _minuit  = NULL;
    _dar     = new DoubletAmbigResolver(_darPset,0,0,0);

  }


//-----------------------------------------------------------------------------
  AvikPIDNew::~AvikPIDNew() {
    if (_minuit) delete _minuit;
    delete _dar;
  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::beginJob() {

    // histograms

    art::ServiceHandle<art::TFileService> tfs;

    if (_diagLevel) {
      _pidtree = tfs->make<TTree>("PID", "PID info");

      _pidtree->Branch("trkid"         , &_trkid         , "trkid/I");
      _pidtree->Branch("p"             , &_trkmom        , "trkmom/D");
      _pidtree->Branch("drdsVadim"     , &_drdsVadim     , "drdsVadim/D");
      _pidtree->Branch("drdsVadimErr"  , &_drdsVadimErr  , "drdsVadimErr/D");
      _pidtree->Branch("logDedxProbEle", &_logDedxProbEle, "logDedxProbEle/D");
      _pidtree->Branch("logDedxProbMuo", &_logDedxProbMuo, "logDedxProbMuo/D");
    }

    _minuit = new TMinuit(2);

  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::beginRun(art::Run & run){
    if (_debugLevel >= 2) cout << "AvikPIDNew: From beginRun: " << run.id().run() << endl;


  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::beginSubRun(art::SubRun & lblock ) {
    if (_debugLevel>=2) cout << "AvikPIDNew: From beginSubRun. " << endl;
  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::endJob(){
    if (_debugLevel>=2) cout << "AvikPIDNew: From endJob. " << endl;
  }


//-----------------------------------------------------------------------------
// Avik's weighted residual
//-----------------------------------------------------------------------------
  double AvikPIDNew::weightedResidual(double R) {
    double wr(1.e10);

    double ar = fabs(R);
    if (ar < 0.2) wr = pow(ar, 1.7);
    else          wr = pow(0.2,1.3)*pow(ar,0.4);

    return wr;
  }

//-----------------------------------------------------------------------------
  double AvikPIDNew::weightedSlopeResidual(double DDrDz) {

    double res(0), ar;

    ar = fabs(DDrDz);

    if (ar < _bound) res = (1/pow(_bound, _pow1))*pow(ar, _pow1);
    else             res = (1/pow(_bound, _pow1))*pow(_bound, (_pow1 -_pow2))*pow(ar, _pow2);

    return res;
  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::doubletMaker(const KalRep* Trk) {

    art::Handle<mu2e::KalRepPtrCollection> krepsHandle;

    TrkHitVector const& hot_list = Trk->hitVector();
//-----------------------------------------------
// ELECTRONS:
//-----------------------------------------------
    int nhits  = 0;  // total number of hits
    int ncount = 0;  // number of hits in doublets

    //    for (auto it=ele_hot_list.begin(); it<ele_hot_list.end(); it++) ele_nhits+=1;

    nhits = hot_list.size();

    int   iamb0(0);			// number of hits with undefined drift direction

    //    float resall   [nhits];
    int   resgood  [nhits];
    //    float res      [nhits];
    int   panelall [nhits];
    int   planeall [nhits];
    int   layall   [nhits];
    int   Nall     [nhits];
    //    int   strawall [nhits];
    //    int   straw    [nhits];
    int   iamball  [nhits];
    int   nmatchall[nhits] = {0};   // number of matches for each hit

    int   nnlet[6] = {nhits};       // array giving number of singlets, doublets, triplets, ...

    Hep3Vector pos;
    double     hitres, hiterr;

    for (int i=0; i<nhits; i++) {

      const mu2e::TrkStrawHit* hit = dynamic_cast<const mu2e::TrkStrawHit*>(hot_list.at(i));

      if (hit) {
	mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

	panelall[i] = straw->id().getPanel();
	planeall[i] = straw->id().getPlane();
	layall  [i] = straw->id().getLayer();
	Nall    [i] = straw->id().getStraw();
	iamball [i] = hit->ambig();
	resgood [i] = hit->resid(hitres,hiterr,1);
      }
    }
//-----------------------------------------------------------------------------
// loop over the track hits
//-----------------------------------------------------------------------------
    for(int i=0; i<nhits; i++) {

      if (iamball[i] == 0) iamb0 += 1;

      for(int j=0; j<i; j++) {     //this loop checks current hit against all previous hits in list

        if((panelall[i]==panelall[i-(j+1)]) &&
           (planeall[i]==planeall[i-(j+1)]) &&
           (layall[i]!=layall[i-(j+1)]) &&
           (abs(Nall[i]-Nall[i-(j+1)])<=2) &&
           (iamball[i]!=iamball[i-(j+1)]) &&
           (resgood[i]==1) && (resgood[i-(j+1)]==1)) {   // doublet conditions

          nmatchall[i]+=1;    // number of matches to other hits
          if((nmatchall[i-(j+1)]==0) && (nmatchall[i]==1)) {  // i.e. has found a doublet
	    //            res[ncount+1]=resall[i];   // i is current hit
	    //            straw[ncount+1]=strawall[i];
	    //            res[ncount]=resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    //            straw[ncount]=strawall[i-(j+1)];
            ncount+=2;
            nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            nnlet[1]+=1;
            nnlet[0]-=2;
          }

          else if (nmatchall[i-(j+1)]==0) {
	    //            res[ncount]=resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    //            straw[ncount]=strawall[i-(j+1)];
            ncount+=1;
            nnlet[nmatchall[i]]+=1;
            nnlet[nmatchall[i]-1]-=1;
            nmatchall[i-(j+1)]+=1;    //number of matches to other hits
            nnlet[0]-=1;
          }

          else {
	    //            res[ncount]=resall[i];   // i is current hit
	    //            straw[ncount]=strawall[i];
            ncount+=1;
            nnlet[nmatchall[i-(j+1)]+nmatchall[i]]+=1;
            nnlet[nmatchall[i-(j+1)]]-=1;
            nnlet[nmatchall[i]-1]-=1;
            nmatchall[i-(j+1)]+=1;    //number of matches to other hits
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// ANALYSIS:
//-----------------------------------------------------------------------------
    float res    [ncount];
    float res_sq [ncount];

    float res_all [nhits] = {0};
    float res_sq2[nhits] = {0};

    float res_sum  = 0;
    float res_sum2 = 0;

    for(int i=0; i<ncount; i++) {
      res   [i] = 0.0;
      res_sq[i] = 0.0;
    }

    for (int i=0; i<nhits; i++) {
      res_all[i] = 0.;
      res_sq2[i] = 0.;
    }

    for (int i=0; i<ncount; i++) {
      if (abs(res[i])>1.5) {   // residual is too large
        res[i]=0;
      }
    }

    for (int i=0; i<nhits; i++) {
      if (abs(res_all[i])>1.5) {   // if for either hypothesis residuals are too high
        res_all[i]=0;
      }
    }

    float matchhits = 0.0;
    float matchhits_all = 0.0;

    for(int i=0; i<ncount; i++) {
      if (res[i]!=0) matchhits+=1.0;

      res_sq[i] = weightedResidual(res[i]);
      res_sum  += res_sq[i];
    }

    for(int i=0; i<nhits; i++) {
      if (res_all[i] != 0) matchhits_all+=1.0;

      res_sq2[i] = weightedResidual(res_all[i]);
      res_sum2  += res_sq2[i];
    }

    _sumAvik     = res_sum;
    _sq2Avik     = res_sum2;
    _nMatched    = matchhits;
    _nMatchedAll = matchhits_all;

    if (_debugLevel > 0) {
      printf( "res_sum  is:  %8.4f res_sum2 is:  %8.4f \n",res_sum,res_sum2);
    }
  }

//-----------------------------------------------------------------------------
  double AvikPIDNew::calculateDedxProb(std::vector<double>* GasPaths ,
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
  bool AvikPIDNew::calculateVadimSlope(const KalRep  *KRep  ,
				       double        *Slope ,
				       double        *Eslope) {
    
    std::vector<double>   res, flt, eres, eflt;
    mu2e::TrkStrawHit*    hit;
    double                resid, residerr, aresd, normflt, normresd;
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
  int AvikPIDNew::AddHits(const Doublet* Multiplet, vector<double>& Fltlen, vector<double>& Resid) {
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
  int AvikPIDNew::AddSsMultiplets(const vector<Doublet>* ListOfDoublets,
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
  int AvikPIDNew::AddOsMultiplets(const vector<Doublet>* ListOfDoublets,
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
  int AvikPIDNew::CalculateSlope(vector<double>& Fltlen, vector<double>& Res, double& Slope, double& Err) {

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
  void AvikPIDNew::calculateSsSums(const vector<Doublet>* ListOfDoublets, double&  Drds, double& DrdsErr, int& NUsed) {

    vector<double> fltLen, resid;

    AddSsMultiplets(ListOfDoublets,fltLen,resid);
    NUsed = fltLen.size();
    CalculateSlope(fltLen,resid,Drds,DrdsErr);
  }
//------------------------------------------------------------------------------
// calculate sums over the local doublet residuals
// residuals are weighted, as Avik is trying to de-weight the tails
//-----------------------------------------------------------------------------
  void AvikPIDNew::calculateOsSums(const vector<Doublet>* ListOfDoublets,
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
            double dr2     = weightedSlopeResidual(ddxdz);
            Sum           += dr2;
            NUsedDoublets += 1.;
          }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
  void AvikPIDNew::produce(art::Event& event) {

    auto detmodel = _mu2eDetector_h.getPtr(event.id());

    art::Handle<mu2e::KalRepPtrCollection> handle;

    double         path, dedx_prob_ele, dedx_prob_muo;

    int const      max_ntrk(100);
    int            n_trk; 

    vector<Doublet>                    listOfDoublets;

    vector<double>         gaspaths;
    vector<double>         edeps;
    const KalRep            *trk;

    const TrkHitVector*      hots;
    const TrkStrawHit*       hit ;


    _evtid = event.id().event();

    ++_processed_events;

    if (_processed_events%100 == 0) {
      if (_debugLevel >= 1) cout << "AvikPIDNew: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;
    }

    if (_debugLevel >= 2) cout << "AvikPIDNew: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;

    unique_ptr<AvikPIDNewProductCollection> pids(new AvikPIDNewProductCollection );

    art::Selector  selector(art::ProcessNameSelector("") && art::ModuleLabelSelector(_trkRecModuleLabel));

    event.get(selector,handle);

    if (! handle.isValid()) {
      printf("TAnaDump::printKalRepPtrCollection: no KalRepPtrCollection for module, BAIL OUT\n");
      goto END;
    }

    _listOfTracks = handle.product();

    n_trk = _listOfTracks->size();

    if (_debugLevel > 0) {
      printf("Event: %8i : n_trk: %2i\n",_evtid,n_trk);
    }
//-----------------------------------------------------------------------------
// proceed further
//-----------------------------------------------------------------------------
//    for (int i=0; i<max_ntrk; i++) unique[i] = 1;

    if (n_trk > max_ntrk/2) {
      printf("Event: %8i : n_trk: %2i BAIL OUT\n", _evtid,n_trk);
      goto END;
    }

    for (int i=0; i<n_trk; i++) {
      _trkid     = i;
      trk        = _listOfTracks->at(i).get();
      hots       = &trk->hitVector();
      int nh     = hots->size();
//-----------------------------------------------------------------------------
// track hit doublets
//-----------------------------------------------------------------------------
      _dar->findDoublets(trk,&listOfDoublets);
//-----------------------------------------------------------------------------
// dE/dX: use only 'active' hits
// calculate dE/dX, clear vectors, start forming a list of hits from the track
//-----------------------------------------------------------------------------
      double        firsthitfltlen(1.e6), lasthitfltlen(1.e6), entlen;

      const TrkHit  *first(nullptr), *last(nullptr); 

      for (int ih=0; ih<nh; ++ih) {
        const TrkHit* hit =  dynamic_cast<const TrkHit*> (hots->at(ih));
      	if (hit   != nullptr) {
      	  if (first == nullptr) first = hit;
      	  last = hit;
      	}
      }

      // first = dynamic_cast<const TrkHit*> (trk->firstHit()->kalHit()->hit());
      // last  = dynamic_cast<const TrkHit*> (trk->lastHit ()->kalHit()->hit());

      // if (dynamic_cast<const TrkStrawHit*> (first) == nullptr) { 
      // 	printf("ERROR in AvikPIDNew::produce for Event: %8i : first hit is not a TrkStrawHit, test fltLen*\n",_evtid);
      // 	double len = first->fltLen();
      // 	printf("first->fltLen() = %10.3f\n",len);
      // }

      // if (dynamic_cast<const TrkStrawHit*> (last ) == nullptr) { 
      // 	printf("ERROR in AvikPIDNew::produce for Event: %8i : last  hit is not a TrkStrawHit, test fltLen*\n",_evtid);
      // 	double len = last->fltLen();
      // 	printf("last->fltLen() = %10.3f\n",len);
      // }

      if (first) firsthitfltlen = first->fltLen() - 10;
      if (last ) lasthitfltlen  = last->fltLen()  - 10;

      entlen  = std::min(firsthitfltlen,lasthitfltlen);
      _trkmom = trk->momentum(entlen).mag();

      gaspaths.clear();
      edeps.clear();

      for (int ih=0; ih<nh; ++ih) {
        hit =  dynamic_cast<TrkStrawHit*> (hots->at(ih));
        if (hit && hit->isActive()) {
//-----------------------------------------------------------------------------
// hit charges: '2.*' here because KalmanFit reports half-path through gas.
//-----------------------------------------------------------------------------
	  const Straw* straw = &hit->straw();
          const DetStrawElem* strawelem = detmodel->strawElem(*straw);
          path = 2.*strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction(hit->fltLen()));
          gaspaths.push_back(path);
          edeps.push_back(hit->comboHit().energyDep());
        }
      }

      dedx_prob_ele   = calculateDedxProb(&gaspaths, &edeps, _heletemp);
      _logDedxProbEle = log(dedx_prob_ele);
      dedx_prob_muo   = calculateDedxProb(&gaspaths, &edeps, _hmuotemp);
      _logDedxProbMuo = log(dedx_prob_muo);
//-----------------------------------------------------------------------------
// calculate ddR/ds slope for the electron tracks
//-----------------------------------------------------------------------------
      calculateVadimSlope(trk,&_drdsVadim,&_drdsVadimErr);

      calculateOsSums(&listOfDoublets,_drdsOs,_drdsOsErr,_nusedOsH,_resSumOs,_nusedOsD);
      calculateSsSums(&listOfDoublets,_drdsSs,_drdsSsErr,_nusedSsH);
//-----------------------------------------------------------------------------
// calculate Avik's sums
//-----------------------------------------------------------------------------
      doubletMaker(trk);
//-----------------------------------------------------------------------------
// save...
//-----------------------------------------------------------------------------
      _pid.init(_trkid         , _nMatched       , _nMatchedAll   ,
		_nusedOsH      , _nusedSsH   , _nusedOsD  , 
		_logDedxProbEle, _logDedxProbMuo ,
		_drdsVadim  , _drdsVadimErr,
		_drdsOs     , _drdsOsErr   ,
		_drdsSs     , _drdsSsErr   ,
		_sumAvik    , _sq2Avik     , _resSumOs);
      
      pids->push_back(_pid);
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
      if (_diagLevel > 0) _pidtree->Fill();
    }
//-----------------------------------------------------------------------------
// end of the routine
//-----------------------------------------------------------------------------
  END:;
    event.put(std::move(pids));
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::AvikPIDNew);
