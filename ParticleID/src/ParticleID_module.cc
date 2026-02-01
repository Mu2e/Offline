///////////////////////////////////////////////////////////////////////////////
//  $Id:
//  $Author:
//  $Date:
//
//  Original author Vadim Rusu
//
// 2015-01-23 P.Murat: default condiguration is stored in ParticleID/fcl/prolog.fcl
///////////////////////////////////////////////////////////////////////////////
// C++ includes.
#include <iostream>
#include <string>
#include <sstream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
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
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TStyle.h"

#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/PIDProduct.hh"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetector.hh"

#include "Offline/ParticleID/inc/PIDUtilities.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

using namespace std;

namespace mu2e {

TGraphErrors *error;

  //compute sum of squares of residuals
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
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

  static const int nbounds=11;
float pathbounds[nbounds]= {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

int findlowhist(float d){

  if (d<=pathbounds[0]) return 0;

  else if (d>pathbounds[0] && d<=pathbounds[1]) return 0;
  else if (d>pathbounds[1] && d<=pathbounds[2]) return 1;
  else if (d>pathbounds[2] && d<=pathbounds[3]) return 2;
  else if (d>pathbounds[3] && d<=pathbounds[4]) return 3;
  else if (d>pathbounds[4] && d<=pathbounds[5]) return 4;
  else if (d>pathbounds[5] && d<=pathbounds[6]) return 5;
  else if (d>pathbounds[6] && d<=pathbounds[7]) return 6;
  else if (d>pathbounds[7] && d<=pathbounds[8]) return 7;
  else if (d>pathbounds[8] && d<=pathbounds[9]) return 8;
  else if (d>pathbounds[9] && d<=pathbounds[10]) return 9;
  else if (d>pathbounds[10]) return 10;

  else
    {    cout<<"Out of  bounds. Should never end up here\n"; return -9999.;}

}

  class ParticleID : public art::EDProducer {

  public:
    explicit ParticleID(fhicl::ParameterSet const& pset);
    virtual ~ParticleID() { }
    void beginJob();
    void beginRun(art::Run &run);
    void beginSubRun(art::SubRun & lblock );
    virtual void produce(art::Event& event);
    void endJob();

  private:

    PIDProduct _pid;

    int _processed_events;
    int _evtid;

    std::string _fitterModuleLabel;
    std::string _electronTemplateFile;
    std::string _muonTemplateFile;

    ProditionsHandle<Mu2eDetector> _mu2eDetector_h;

    TrkParticle _tpart;
    TrkFitDirection _fdir;

    int _trkid;
    double _trkmom;
    double _residualsSlope;
    double _residualsSlopeError;
    double _logeprob;
    double _logmprob;

    int  _debugLevel;
    int  _verbosity;
    int  _diagLevel;
    bool _doDisplay;

    std::string _electrontemplates;
    std::string _muontemplates;

    TH1D* _heletemp[nbounds];
    TH1D* _hmuotemp[nbounds];

    int   _templatesnbins ;
    float _templateslastbin ;
    float _templatesbinsize ;

    TTree *       _pidtree;
    TCanvas*      _plotCanvas;

    bool calculateSlope(std::vector<double>vresd,std::vector<double>vflt,
                        std::vector<double>evresd,std::vector<double>evflt,
                        double * slope,
                        double * eslope);

    double calculateProb(std::vector<double>gaspaths, std::vector<double>edeps, TH1D** templates);

    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob.
    TDirectory* _directory;
  };

  ParticleID::ParticleID(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _fitterModuleLabel   (pset.get<string>("fitterModuleLabel")),
    _electronTemplateFile(pset.get<string>("ElectronTemplates")),
    _muonTemplateFile    (pset.get<string>("MuonTemplates"    )),
    _tpart    ((TrkParticle::type)            (pset.get<int>("fitparticle" ))),
    _fdir     ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _debugLevel(pset.get<int>("debugLevel", 0   )),
    _verbosity(pset.get<int> ("verbosity" , 0   )),
    _diagLevel(pset.get<int> ("diagLevel" , 0   )),
    _doDisplay(pset.get<bool>("doDisplay" ,false)),
    _pidtree(0),
    _plotCanvas(0)
  {
    _processed_events = -1;

    produces<PIDProductCollection>();

    // location-independent files
    ConfigFileLookupPolicy configFile;

    _electrontemplates = configFile(_electronTemplateFile);
    _muontemplates     = configFile(_muonTemplateFile    );

    char name[50];
    TFile* electrontemplatefile = TFile::Open(_electrontemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempe%d",i);
      electrontemplatefile->GetObject(name,_heletemp[i]);
    }

    TFile* muontemplatefile = TFile::Open(_muontemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempm%d",i);
      muontemplatefile->GetObject(name,_hmuotemp[i]);
    }

    _templatesnbins = _heletemp[0]->GetNbinsX();
    _templateslastbin = _heletemp[0]->GetBinLowEdge(_templatesnbins)+_heletemp[0]->GetBinWidth(1);
    _templatesbinsize = _heletemp[0]->GetBinWidth(1);

  }

  void ParticleID::beginJob(){

    // histograms

    art::ServiceHandle<art::TFileService> tfs;

    if (_diagLevel) {
      _pidtree = tfs->make<TTree>("PID", "PID info");
      _pidtree->Branch("trkid"   , &_trkid              , "trkid/I");
      _pidtree->Branch("slope"   , &_residualsSlope     , "slope/D");
      _pidtree->Branch("errslope", &_residualsSlopeError, "errslope/D");
      _pidtree->Branch("p"       , &_trkmom             , "trkmom/D");
      _pidtree->Branch("logeprob", &_logeprob           , "logeprob/D");
      _pidtree->Branch("logmprob", &_logmprob           , "logmprob/D");
    }

    if(_doDisplay) {
      // If needed, create the ROOT interactive environment. See note 1.
      if ( !gApplication ){
        int    tmp_argc(0);
        char** tmp_argv(0);
        _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
      }

      gStyle->SetPalette(1);
      gROOT->SetStyle("Plain");
      gStyle->SetOptFit(1);
      gStyle->SetMarkerStyle(22);

      _plotCanvas = new TCanvas("plots","PID Plots",600,400);
    }
  }

  void ParticleID::beginRun(art::Run & run){
    if (_verbosity>=2) cout << "ParticleID: From beginRun: " << run.id().run() << endl;

  }

  void ParticleID::beginSubRun(art::SubRun & lblock ) {
    if (_verbosity>=2) cout << "ParticleID: From beginSubRun. " << endl;
  }

  void ParticleID::endJob(){
    if (_verbosity>=2) cout << "ParticleID: From endJob. " << endl;
    if (_doDisplay) delete _plotCanvas;

  }

////////// Produce ///////////

  void ParticleID::produce(art::Event& event) {

    auto detmodel = _mu2eDetector_h.getPtr(event.id());

    _evtid = event.id().event();
    ++_processed_events;
    if (_processed_events%100 == 0) {
      if (_verbosity>=1) cout << "ParticleID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;
    }
    if (_verbosity>=2) cout << "ParticleID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;

    unique_ptr<PIDProductCollection> pids(new PIDProductCollection );

   art::Handle<KalRepPtrCollection> trksHandle;
   event.getByLabel(_fitterModuleLabel,trksHandle);
   const KalRepPtrCollection* const trks = trksHandle.product();

   if (!trksHandle.isValid()) {
     if (_verbosity>=1) cout << "ParticleID : " << "no" << " obj for " << _fitterModuleLabel.c_str() << " of event " << _evtid << endl;
   }

   if (trks->size() >0) {
     if (_verbosity>=1) cout << "ParticleID : " << trks->size() << " obj for " << _fitterModuleLabel.c_str() << " of event " << _evtid << endl;
   }

   for ( size_t i=0; i< trks->size(); ++i ){

     _trkid = i;
     const KalRep* krep = trks->at(i).get();

     double firsthitfltlen = krep->firstHit()->kalHit()->hit()->fltLen() - 10;
     double lasthitfltlen = krep->lastHit()->kalHit()->hit()->fltLen() - 10;
     double entlen = std::min(firsthitfltlen,lasthitfltlen);
     _trkmom = krep->momentum(entlen).mag();

     std::vector<KalSite*> kalsites= krep->siteList();

     std::vector<double> vresd;
     std::vector<double> vflt;
     std::vector<double> evflt;
     std::vector<double> evresd;
     std::vector<double> gaspaths;
     std::vector<double> edeps;

     for(unsigned isite=0;isite<kalsites.size();isite++){
       KalSite* ksite = kalsites[isite];

       if (ksite->type() == KalSite::hitSite){

         TrkStrawHit* hit = dynamic_cast<TrkStrawHit*>(ksite->kalHit()->hit());
         double resid, residerr;
         hit->resid(resid,residerr,true);

         bool activehit = hit->isActive();
         if (activehit){
           double aresd = (hit->poca().doca()>0?resid:-resid);
           double normflt = hit->fltLen() -  krep->flt0();
           double normresd = aresd/residerr;

           vresd.push_back(normresd);
           vflt.push_back(normflt);
           evresd.push_back(1.);
           evflt.push_back(0.1);

           //      2. * here because KalmanFit reports half the path through gas.
           const DetStrawElem* strawelem = detmodel->strawElem(hit->straw());
           gaspaths.push_back(2. * strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction( hit->fltLen() )));

           edeps.push_back(hit->comboHit().energyDep());

         }
       }
     }

     calculateSlope(vresd,vflt,evresd,evflt,&_residualsSlope,&_residualsSlopeError);

     double eprob = calculateProb(gaspaths, edeps, _heletemp);
     double muprob = calculateProb(gaspaths, edeps, _hmuotemp);

     _logeprob = log(eprob);
     _logmprob = log(muprob);

     if(_diagLevel)
       _pidtree->Fill();

     _pid.clear();
     _pid.SetTrkID(_trkid);
     _pid.SetResidualsSlope(_residualsSlope);
     _pid.SetResidualsSlopeError(_residualsSlopeError);
     _pid.SetLogEProb(_logeprob);
     _pid.SetLogMProb(_logmprob);

     pids->push_back(_pid);

   }  // end of trks loop
   event.put(std::move(pids));

  }

  double ParticleID::calculateProb(std::vector<double>gaspaths, std::vector<double>edeps, TH1D** templates){

    static const double _minpath = 0.5;
    static const double _maxpath = 10.;

    double thisprob = 1;

    for (unsigned int ipath = 0; ipath < gaspaths.size(); ipath++){
      double thispath = gaspaths.at(ipath);
      double thisedep = edeps.at(ipath);

      double tmpprob = 0;
      if (thispath > _minpath && thispath<=_maxpath){
          int lowhist = findlowhist(thispath);

          PIDUtilities util;
          TH1D* hinterp = util.th1dmorph(templates[lowhist],templates[lowhist+1],pathbounds[lowhist],pathbounds[lowhist+1],thispath,1,0);
          //what is the probability for this edep
          int thisedepbin = -999;
          if (thisedep > _templateslastbin) thisedepbin = _templatesnbins;
          else
            thisedepbin = int(thisedep/_templatesbinsize)+1;

          tmpprob=hinterp->GetBinContent(thisedepbin);

          if (_doDisplay){
            templates[lowhist]->Draw();
            templates[lowhist+1]->Draw("same");
            hinterp->SetLineColor(2);
            hinterp->Draw("same");
            _plotCanvas->WaitPrimitive();
          }
          hinterp->Delete();
        }

      if (tmpprob> 0)
        thisprob = thisprob * tmpprob;

    }

    return thisprob;

  }

  bool ParticleID::calculateSlope(std::vector<double>vresd,std::vector<double>vflt,std::vector<double>evresd,std::vector<double>evflt,  double * slope, double * eslope) {

    error = new TGraphErrors(vresd.size(),vflt.data(),vresd.data(),evflt.data(),evresd.data());

    TMinuit *gmMinuit = new TMinuit(2);
    gmMinuit->SetPrintLevel(-1);
    gmMinuit->SetFCN(myfcn);
    const int dim(2);
    const char par_name[dim][20]={"offset","slope"};
    static Double_t step[dim] = {0.001,0.001};
    Double_t sfpar[dim]={0.0,0.005};
    Double_t errsfpar[dim]={0.0,0.0};
    int ierflg = 0;
    for (int ii = 0; ii<dim; ii++) {
      gmMinuit->mnparm(ii,par_name[ii],sfpar[ii], step[ii], 0,0,ierflg);
    }
    gmMinuit->FixParameter(0);
    gmMinuit->Migrad();
    bool converged = gmMinuit->fCstatu.Contains("CONVERGED");
    if (!converged)
      {
        cout <<"-----------TOF Linear fit did not converge---------------------------" <<endl;
        return converged;
      }
    for (int i = 0;i<dim;i++) {
      gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
    }

    *slope = sfpar[1];
    *eslope = errsfpar[1];

    if (_doDisplay){
      error->Draw("AP");
      _plotCanvas->WaitPrimitive();
    }

    delete error;

    return converged;
  }

} // end namespace mu2e

using mu2e::ParticleID;
DEFINE_ART_MODULE(ParticleID)
