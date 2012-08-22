//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: kalmanFitv2_module.cc,v 1.1 2012/08/22 17:30:37 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/08/22 17:30:37 $
//

#include <iostream>
#include <string>
#include <memory>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/ITracker.hh"
//#include "GeometryService/inc/WorldG4.hh"
//#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"

#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

#include <iostream>
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
//#include "TrkBase/TrkRecoTrk.hh"
//#include "BaBar/PdtPid.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "DchGeom/DchDetector.hh"
#include "BaBar/ErrLog.hh"
#include "BField/BFieldFixed.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetSet.hh"

#include "KalmanTestsI/inc/DchGDchP.hh"
#include "KalmanTestsI/inc/KalFitI.hh"

#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <fstream>

#include "CLHEP/Random/RandGaussQ.h"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

using namespace std; 

namespace mu2e 
{
  class kalmanFitv2 : public art::EDProducer
  {
    public:
    explicit kalmanFitv2(fhicl::ParameterSet const&);
    virtual ~kalmanFitv2() { }
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    void endJob();
    void produce(art::Event& e);

    private:
    CLHEP::RandGaussQ _gaussian;


    // Label of the module that created the data products.
    std::string _g4ModuleLabel;
    // Label of the generator.
    std::string _generatorModuleLabel;
    // Instance names of data products
    std::string _targetStepPoints;


    string _hitMakerModuleLabel;

    
    double _cutdp;
    int _oneturn;
    double _aveSpcRes;
    bool _usemchitdist;
    
    int _debugLvl;

    DchDetector* dchdet;
    std::string _materialdb;
    KalFitI _kfit;

    // End: run time parameters

    bool FillMCInfo(art::Event& event);

    std::vector<hitIndex> strawhits;
    HelixTraj seed;
    double time0;
    double s0;
    std::vector<double> hitflt;

    CLHEP::Hep3Vector gen_pos;
    CLHEP::Hep3Vector gen_mom;

    // Pointers to histograms, ntuples, TGraphs.
    void FillHistos(TrkDef& trkdef,TrkKalFit& myfit); 

    int nbadfit;
    TH2D *hpull;
    TH2D *hdev;
    TH1D *hchi2;
    TH1D *hprob;
    TH1D *hchi2_0;
    TH1D *hmom;
    TH1D *hpars[5];
    TH1D *hnhits;
    TH1D *hnhitstotal;
    TH1D *hnturns;
    TH1D *heloss;
    TH1D *hdistres;
    TH1D *htargeteloss;
    TTree *treefit;
    TrkExchangePar startpar;
    TrkExchangePar recopar;



    struct {
      int fit;
      float chi2,chi2out;
      int nhits,nhitstot,nturn;
      float momin,momout;
      float fitmom, fitmomerr;
      float t0,t0fit,errt0;
      void clear(){
	fit=0;chi2=-1;chi2out=-1;nhits=0;nhitstot=0;nturn=0;momin=0;momout=0;
	fitmom=0; fitmomerr=0;t0=0;t0fit=0;errt0=-1;};
    } recoinfo;
    struct {
      float elosstgt;
      int ntgt;
      void clear(){elosstgt=0;ntgt=0;};
    } eventinfo;


  };

  kalmanFitv2::kalmanFitv2(fhicl::ParameterSet const& pset) :
    _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget")),
    _hitMakerModuleLabel(pset.get<std::string>("hitMakerModuleLabel", "makeDcH")),
    _cutdp(pset.get<double>("cutdp",1000.)),
    _oneturn(pset.get<int>("oneturn",0)),
    _aveSpcRes(pset.get<double>("aveSpcRes",0.2)),
    _usemchitdist(pset.get<bool>("mchitdist",false)),
    _debugLvl(pset.get<int>("debuglevel",1)),
    _materialdb(pset.get<std::string>("materialdb", "")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    seed(TrkParams(HelixTraj::NHLXPRM)),
    nbadfit(0),
    startpar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1)),
    recopar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1))
  {
    produces<CLHEP::HepVector>("test");
  }

  void kalmanFitv2::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    
    hpull  =tfs->make<TH2D>("hpull","pull;param;#delta parameter/#sigma",5,0,5,10000,-500,500);
    hdev   =tfs->make<TH2D>("hdev","dev",5,0,5,200000,-1,1);
    hchi2  =tfs->make<TH1D>("hchi2","chi2/ndf;#chi^{2}/ndf",100000,0,1000);
    hprob  =tfs->make<TH1D>("hprob","prob",10000,0,1);
    hchi2_0=tfs->make<TH1D>("hchi2_0","chi2",100000,0,1000);
    hmom   =tfs->make<TH1D>("hmom",";#delta P, MeV/c",10000,-100,100);

    hpars[0]=tfs->make<TH1D>("d0",";d0, mm",400,-30,30);
    hpars[1]=tfs->make<TH1D>("phi0",";#phi0, rad",400,-0.05,0.05);
    hpars[2]=tfs->make<TH1D>("omega",";#omega, 1/mm",400,-2e-4,2e-4);
    hpars[3]=tfs->make<TH1D>("z0",";z0, mm",400,-30,30);
    hpars[4]=tfs->make<TH1D>("tan",";tan",400,-0.05,0.05);

    hnhits     =tfs->make<TH1D>("hnhits",";number of hits per turn",500,0,500);
    hnhitstotal=tfs->make<TH1D>("hnhitstotal",";total number of hits",500,0,500);
    hnturns    =tfs->make<TH1D>("hnturns",";number of turns",10,0,10);
    heloss     =tfs->make<TH1D>("heloss",";E loss from first to last cell, MeV",10000,0,10);
    hdistres   =tfs->make<TH1D>("hdistres",";hit dist resol, mm",10000,-10,10);
    htargeteloss     =tfs->make<TH1D>("htargeteloss",";E loss at target, MeV",10000,0,10);

    treefit=tfs->make<TTree>("trfit","track fit params");
    treefit->Branch("startpar","TrkExchangePar",&startpar);
    treefit->Branch("recopar","TrkExchangePar",&recopar);
    treefit->Branch("fitinfo",&recoinfo,"fit/I:chi2/F:chi2out/F:nhits/I:nhitstot/I:nturn/I:momin/F:momout/F:fitmom/F:fitmomerr/F:t0/F:t0fit/F:errt0/F");
    treefit->Branch("eventinfo",&eventinfo,"elosstgt/F:ntgt/I");

    
  }

  void kalmanFitv2::beginRun(art::Run&){
    std::cout<<"init dchdetector"<<std::endl;
    DchGDchP gdch(_materialdb);
    dchdet=new DchDetector(gdch,true);
  }

  void kalmanFitv2::produce(art::Event& event) 
  {
    if(_debugLvl>9) 
      ErrMsg::ErrLoggingLevel = debugging;


    int iev=event.id().event();
    if((iev%1000)==0)cout<<"event="<<iev<<endl;

    art::Handle<mu2e::StrawHitMCTruthCollection> cellmchits;
    art::Handle<mu2e::StrawHitCollection> cellhits;
    if(!(event.getByLabel(_hitMakerModuleLabel,cellhits)&&event.getByLabel(_hitMakerModuleLabel,cellmchits))){
      cout<<"no dc hits "
	  <<event.getByLabel(_hitMakerModuleLabel,cellhits)<<" "
	  <<event.getByLabel(_hitMakerModuleLabel,cellmchits)<<endl;
      return;
    }

    if(!FillMCInfo(event)) return;

    TrkDef trkdef(cellhits.product(), strawhits,seed,TrkParticle(),TrkFitDirection(),time0,0.0);

    _kfit._flt0=s0;
    _kfit._hitflt=hitflt;
    TrkKalFit myfit;
    _kfit.makeTrack(trkdef,myfit);

    if(myfit._fit.success()){
     TrkHotList* hots = myfit._krep->hotList();
   
     std::cout<<"active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
      std::vector<hitIndex> strawhitsinactive;
      
      _kfit.reActivateHitsbyTurn(trkdef,myfit);
      
      std::cout<<"Re:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
    }

    if(myfit._fit.success()){
      auto_ptr<HepVector> recomom(new HepVector(6));


      KalRep *kalrep = myfit._krep;
      double fltlen=kalrep->firstHit()->globalLength();
      HepVector hvpos(3),hvmom(3);
      HepSymMatrix ehpos(3),ehmom(3);
      HepMatrix ehposmom(3,3);
      kalrep->getAllWeights(fltlen,hvpos,hvmom,ehpos,ehmom,ehposmom);
      for(int i=1;i<=3;i++) {
	(*recomom)(i)=hvpos(i);
	(*recomom)(i+3)=hvmom(i);
      }
    
      event.put(recomom,"test");
    }

    FillHistos(trkdef,myfit);

    
  }

  bool kalmanFitv2::FillMCInfo(art::Event& event){
    int iev=event.id().event();

    GeomHandle<ITracker> itr;
    CellGeometryHandle *itwp = itr->getCellGeometryHandle();
    //const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    startpar=TrkExchangePar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
    recopar=TrkExchangePar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
    eventinfo.clear();

    art::Handle<GenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel, genParticles);

    // Steps in the stopping target
    art::Handle<StepPointMCCollection> stHitsHandle;
    event.getByLabel(_g4ModuleLabel,_targetStepPoints,stHitsHandle);
    StepPointMCCollection const& sthits(*stHitsHandle);


    art::Handle<mu2e::PtrStepPointMCVectorCollection> cellptr;
    event.getByLabel(_hitMakerModuleLabel,"StrawHitMCPtr",cellptr);
    PtrStepPointMCVectorCollection const* hits_mcptr = cellptr.product();
    if(hits_mcptr->size()==0) {
      cout<<"Empty ev = "<<event.id().event()<<endl;
      return false;
    }

    art::Handle<mu2e::StrawHitMCTruthCollection> cellmchits;
    art::Handle<mu2e::StrawHitCollection> cellhits;
    if(!(event.getByLabel(_hitMakerModuleLabel,cellhits)&&event.getByLabel(_hitMakerModuleLabel,cellmchits))){
      cout<<"no dc hits "
	  <<event.getByLabel(_hitMakerModuleLabel,cellhits)<<" "
	  <<event.getByLabel(_hitMakerModuleLabel,cellmchits)<<endl;
      return false;
    }

    //select first mc
    int firstmcstep=-1;
    double time0_min=1e10;
    for(unsigned int i=0;i<hits_mcptr->size();i++){
      const PtrStepPointMCVector&  mcptr=hits_mcptr->at(i);
      double tcell=mcptr.at(0)->time();
      if(mcptr.size()&&
	 mcptr.at(0)->trackId().asInt()==1&&
	 tcell<time0_min){time0_min=tcell;firstmcstep=i;}
    }
    if(firstmcstep<0){
      cout<<"ev = "<<event.id().event()
          <<" Not electron in all cells, ncell="<<hits_mcptr->size();
      if(hits_mcptr->size()&&hits_mcptr->at(0).size())
	cout<<" id= "<< hits_mcptr->at(0).at(0)->trackId().asInt();
      cout<<endl;
      return false;
    }
    
    const PtrStepPointMCVector& mcptr=hits_mcptr->at(firstmcstep);
    
    double sumEDep=0.;
    if (true/*_fromOrigin*/) {
      gen_pos = mcptr.at(0)->simParticle()->genParticle()->position();
      //      gen_pos.setX(gen_pos.getX());
      gen_pos.setZ(gen_pos.getZ()+(12000-itr->z0()));
      gen_mom = mcptr.at(0)->simParticle()->genParticle()->momentum().vect();
      //cout<<"gen pos "<<gen_pos<<"gen mom "<<gen_mom<<endl;


      set<int> hitFoils;
      for ( StepPointMCCollection:: const_iterator i=sthits.begin();
	    i !=sthits.end(); ++i ){
	StepPointMC const& hit = *i;
	//	cout<<hit.trackId().asInt()<<" "<<hit.simParticle()->generatorIndex()<<" "<<hit.eDep()<<endl;
	if ( hit.trackId().asInt()!=1){
	  const SimParticle* sp=hit.simParticle().get();
	  while(sp&&sp->id().asInt()!=1) sp = sp->parent().get();
	  if(!sp||sp->id().asInt()!=1) continue;
	}
	sumEDep += hit.eDep();
	hitFoils.insert(hit.volumeId());
      }
      htargeteloss->Fill(sumEDep);
      eventinfo.elosstgt=sumEDep;
      eventinfo.ntgt=hitFoils.size();
    }


    time0=mcptr.at(0)->time();
    CLHEP::Hep3Vector pos=mcptr.at(0)->position();
    CLHEP::Hep3Vector mom=mcptr.at(0)->momentum();

    double mom0=mom.mag();

    HepVector startPar(5);
    s0=0;
    double charge(-1.);

    if(_debugLvl>0) {
      cout<<"pos "<<pos<<" "<<mom<<" |p| "<<mom.mag()<<" targetloss "<<sumEDep<<endl;
    }

    TrkHelixUtils::helixFromMom( startPar, s0, 
				 HepPoint(pos.x(),pos.y(),pos.z()),
				 mom,
				 charge,
				 *_kfit.bField());
    HepSymMatrix dummy(5,1); 
    dummy(1,1)=1.; dummy(2,2)=0.1*0.1;dummy(3,3)=1e-2*1e-2;
    dummy(4,4)=1.; dummy(5,5)=0.1*0.1;
    seed=HelixTraj(startPar,dummy);
    startpar=startPar;
    double dt0atz0=seed.zFlight(0.0)/Constants::c;
    if(_debugLvl>1) {
      std::cout<<"s0 "<<s0<<" time0 "<<time0<<std::endl;
      seed.printAll(std::cout);
    }
    double mommin=mom0;
    

    unsigned int nn=cellhits->size();
    strawhits.clear();

    if(_debugLvl>1)std::cout<<"create hotlist"<<std::endl;
    std::vector<std::pair<double,unsigned int> > sorthits;
    hitflt.clear();
    for(unsigned int i=0;i<nn;++i){
        const mu2e::StepPointMC& stepmc=*cellptr->at(i).at(0);

	double time=stepmc.time();
	double ll=(time-time0)*Constants::c;
	sorthits.push_back(std::pair<double,unsigned int>(ll,i));
	// if(stepmc.trackId().asInt()==1){
	//   strawhits.push_back(hitIndex(i));
	//   hitflt.push_back(ll);
	// }
    }

    std::sort(sorthits.begin(),sorthits.end());
    //sorted array
    double cosphidirprev=0.;
    int nhitstotal=0,nhitsperturn=0,nturn=1,nadd=0;

    for(unsigned int i=0;i<nn;++i){
      int istr=sorthits[i].second;
      const mu2e::StepPointMC& stepmc=*cellptr->at(istr).at(0);
	if(stepmc.trackId().asInt()==1){
	  double mm=stepmc.momentum().mag();
	  double phidir=stepmc.momentum().phi()-stepmc.position().phi();
	  double cosphidir=cos(phidir);

	  if((_oneturn&&cosphidirprev<0&&cosphidir>0)){
	    break;
	  }

	  if(cosphidirprev<0&&cosphidir>0){
	    hnhits->Fill(nhitsperturn);
	    nhitsperturn=0;
	    nturn++;
	  }

	  //	  std::cout<<nhitstotal<<" "<<_oneturn<<" "<<cosphidirprev<<" "<<cosphidir<<std::endl;
	  if(mm<mommin) mommin=mm;
	  nhitsperturn++;
	  nhitstotal++;
	  nadd++;
	  strawhits.push_back(hitIndex(istr));
	  hitflt.push_back(sorthits[i].first);
	  
	  cosphidirprev=cosphidir;
	  
	}
    }
    hnhits->Fill(nhitsperturn);
    hnhitstotal->Fill(nhitstotal);
    hnturns->Fill(nturn);
    heloss->Fill(mom0-mommin);    

    for(unsigned int ii=0;ii<nn;++ii){
      int istr=sorthits[ii].second;
      const mu2e::StepPointMC& stepmc=*cellptr->at(istr).at(0);
      double time=stepmc.time();
      double ll=(time-time0)*Constants::c;
      double dist=cellmchits->at(istr).driftDistance();
      unsigned int dcell_id=cellhits->at(istr).strawIndex().asUint();
      itwp->SelectCellDet(dcell_id);
      const StrawHit& sh = cellhits->at(istr);
      //const Straw& straw = tracker.getStraw(sh.strawIndex());

      int layer=itwp->GetSuperLayer();//dcell_id/10000;
      int cell=itwp->GetWire();//dcell_id%10000;
      if(_debugLvl>1)
	std::cout<<iev<<" "<<istr<<" lay,side,cell "<<layer<<" "<<itwp->isUpStream()<<" "<<cell
		 <<" time "<<time<<" hittime "<<sh.time()<<" dist "<<dist
		 <<" flt "<<ll<<" id "<<stepmc.trackId().asInt()
		 <<" xyz "<<stepmc.position().x()<<" "<<stepmc.position().y()<<" "<<stepmc.position().z()<<std::endl;

      T2D t2d;
      tcal->TimeToDistance(cellhits->at(istr).strawIndex(),(sh.time()-time),stepmc.momentum(),t2d);
      hdistres->Fill(t2d._rdrift-fabs(dist));
      //StrawHit sh2(sh.strawIndex(),distup/0.035+time/*+1.37459e-01/0.035 (wire signal prop)*/,sh.dt(),sh.energyDep());

    }

    recoinfo.clear();
    recoinfo.nhits=-1;
    recoinfo.nhitstot=nhitstotal;
    recoinfo.nturn=nturn;
    recoinfo.momin=mom0;
    recoinfo.momout=mommin;    

    recoinfo.t0=time0-s0/Constants::c+dt0atz0;
    
    return true;
  }

  
  void kalmanFitv2::FillHistos(TrkDef& trkdef,TrkKalFit& myfit){

    KalRep *kalrep = myfit._krep;
    TrkErrCode err= myfit._fit;

    recoinfo.nhits=myfit._krep?myfit._krep->hotList()->nActive():-1;

    if(!err.success()) {
      treefit->Fill();
      //      delete kalrep;
      return;
    }
    
    if(_debugLvl>2){
      std::cout<<"fit kalrep"<<std::endl;
      kalrep->printAll(std::cout);
    }

    double chi2=kalrep->chisquared(trkIn);
    double chi2out=kalrep->chisquared(trkOut);
    if(_debugLvl){
      cout<<"chi2 in "<<chi2<<" out "<<chi2out<<endl;
      //      cout<<"hello err="<<err<<endl;
      // seed.printAll(std::cout);
      cout<<"true t0 helix apprx at z=0 "<<recoinfo.t0
	  <<" t0 at first hit "<<time0<<" fitted at z=0 of helix niters "<<myfit._nt0iter
	  <<" t0 "<<myfit._t00._t0<<" +- "<<myfit._t00._t0err
	  <<" at end "<<myfit._t0._t0<<" +- "<<myfit._t0._t0err<<endl;
    }
    recoinfo.t0fit=myfit._t0._t0;
    recoinfo.errt0=myfit._t0._t0err;

    double Bval = _kfit.bField()->bFieldNominal();

    double fltlen=kalrep->firstHit()->globalLength();
    double fltlenEnd=kalrep->lastHit()->globalLength();
    recopar = kalrep->helix(fltlen);
    HepVector recoparams = recopar.params();
    HepSymMatrix recocovar = recopar.covariance();

    CLHEP::Hep3Vector fitmom = kalrep->momentum(fltlen);
    BbrVectorErr momerr = kalrep->momentumErr(fltlen);
    recoinfo.fitmom = fitmom.mag();
    Hep3Vector momdir = fitmom.unit();
    HepVector momvec(3);
    for(int icor=0;icor<3;icor++)
      momvec[icor] = momdir[icor];
    recoinfo.fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));

    HepVector hvpos(3),hvmom(3);
    HepSymMatrix ehpos(3),ehmom(3);
    HepMatrix ehposmom(3,3);
    kalrep->getAllWeights(fltlen,hvpos,hvmom,ehpos,ehmom,ehposmom);

    const TrkDifPieceTraj& recotraj = kalrep->pieceTraj();
    TrkLineTraj zaxis(HepPoint(0, 0, -10000), Hep3Vector(0, 0, 1), 9000);
    double helixturn=2*TMath::Pi()/seed.omega()/seed.cosDip();
    TrkPoca recpoca(recotraj,fltlen>0?0:-helixturn,zaxis,1000, 1e-12);
    double fltlenbeam=0.;
    if(recpoca.status().success())
      fltlenbeam = recpoca.flt1();
    TrkExchangePar recobeam = kalrep->helix(fltlenbeam);
    HepVector recoparams_beam = recobeam.params();
    HepSymMatrix recocovar_beam = recobeam.covariance();
    double momem_beam=Constants::c*1.0e-3*Bval/recoparams_beam[2]*sqrt(1.+recoparams_beam[4]*recoparams_beam[4]);

    double helixdz=2*TMath::Pi()/recobeam.omega()*recobeam.tanDip();
    kalrep->getAllWeights(fltlenbeam,hvpos,hvmom,ehpos,ehmom,ehposmom);
    if(_debugLvl)
      cout<<"gen pos "<<gen_pos<<" gen mom "<<gen_mom
	  <<" reco ("<<hvpos[0]<<","<<hvpos[1]<<","<<hvpos[2]<<") ("<<hvmom[0]<<","<<hvmom[1]<<","<<hvmom[2]<<")"
	  <<" helix dz "<<helixdz<<endl;



    HepVector recoparamsEnd =  kalrep->helix(fltlenEnd).params();
    HepVector startPar=startpar.params();

    double momem=Constants::c*1.0e-3*Bval/recoparams[2]*sqrt(1.+recoparams[4]*recoparams[4]);
    double momemEnd=Constants::c*1.0e-3*Bval/recoparamsEnd[2]*sqrt(1.+recoparamsEnd[4]*recoparamsEnd[4]);
    double momem0=Constants::c*1.0e-3*Bval/startPar[2]*sqrt(1.+startPar[4]*startPar[4]);
    hmom->Fill(momem-momem0);
    if(_debugLvl){
      cout<<"Fit mom "<<momem<<" min "<<momemEnd<<" beam "<<momem_beam
	  <<" mom0 "<<momem0<<" min "<<recoinfo.momout<<" beam "<<gen_mom.mag()<<endl;
      recopar.printAll(std::cout);
    }
      
    recoinfo.fit=1;
    recoinfo.chi2=chi2;
    recoinfo.chi2out=chi2out;

    treefit->Fill();

    hchi2->Fill(chi2/(recoinfo.nhits-5));
    hprob->Fill(TMath::Prob(chi2,(recoinfo.nhits-5)));
    hchi2_0->Fill(chi2out/(recoinfo.nhits-5));

    double maxpull=0;
    for(int i=0;i<5;i++){
      hpars[i]->Fill(recoparams[i]-startPar[i]);
      hdev->Fill(i,(recoparams[i]-startPar[i]));
      double pull=(recoparams[i]-startPar[i])/TMath::Sqrt(recocovar[i][i]);
      hpull->Fill(i,pull);
      if(fabs(pull)>maxpull) maxpull=fabs(pull);
    }
    
    if(maxpull>100){
      cout<<"problem, big pull value "<<maxpull<<endl;
      nbadfit++;
    }
    
    //delete kalrep;
    //delete hotlist;
    //delete retval;

  }
  
  void kalmanFitv2::endJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    cout<<"nbad "<<nbadfit<<std::endl;

    const char *titlenames[]={"pull d0;#delta parameter/#sigma",
			      "pull #phi0;#delta parameter/#sigma",
			      "pull #omega;#delta parameter/#sigma",
			      "pull z0;#delta parameter/#sigma",
			      "pull tan;#delta parameter/#sigma",
			      "chi2/ndf;#chi^{2}/ndf"};

    TCanvas *c1=tfs->make<TCanvas>("c1pull");
    c1->Divide(2,3,0.001,0.001);
    for(int i=0;i<5;i++){
      c1->cd(i+1);
      TH1D *hp=hpull->ProjectionY(Form("_py%i",i),i+1,i+1);
      hp->SetTitle(titlenames[i]);
      hp->SetAxisRange(-10,10);
      hp->Fit("gaus");
    }
    c1->cd(6);
    hchi2->SetAxisRange(0,10);
    hchi2->Draw();
    c1->Write();
    TCanvas *c2=tfs->make<TCanvas>("c1par");
    c2->Divide(2,3,0.001,0.001);

    for(int i=0;i<5;i++){
      c2->cd(i+1);
      hpars[i]->Fit("gaus");
    }
    c2->Write();
    hmom->SetAxisRange(-4,4);
    hmom->Fit("gaus","","",-0.5,0.75);
  }
}

using mu2e::kalmanFitv2;
DEFINE_ART_MODULE(kalmanFitv2);
