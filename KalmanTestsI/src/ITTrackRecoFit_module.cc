//
// Fit of the reconstructed tracks for the ITracker
//
// $Id: ITTrackRecoFit_module.cc,v 1.2 2012/09/17 14:44:30 ignatov Exp $
// $Author: ignatov $
// $Date: 2012/09/17 14:44:30 $
//
// Original author F. Ignatov and G. Tassielli
//

// C++ includes.
#include <iostream>
//#include <string>
#include <new>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <cstring>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cmath>

#include <boost/shared_ptr.hpp>

// Framework includes.
//#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP includes.
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"
//#include "DataProducts/inc/DPIndexVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "FastPatternReco/inc/FastPatRecoUtilsAndDataDef.hh"
#include "FastPatternReco/inc/TrkHelixFitIT.hh"
#include "FastPatternReco/inc/TrackSeedUtils.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

//temp El Visual
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"


//For track fit
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
//#include "TrkPatRec/inc/TrkHelixFit.hh"
//#include "TrkBase/TrkPoca.hh"
//#include "TrkPatRec/src/TrkPatRec_module.cc"
#include "TrkBase/HelixParams.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkPoca.hh"

#include "DchGeom/DchDetector.hh"

#include "KalmanTestsI/inc/DchGDchP.hh"
#include "KalmanTestsI/inc/KalFitI.hh"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
//#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TLatex.h"
#include "TTree.h"

namespace mu2e {

  class ITTrackRecoFit : public art::EDAnalyzer {
  public:

    explicit ITTrackRecoFit(fhicl::ParameterSet const& pset);
    virtual ~ITTrackRecoFit() {
    }
  
    virtual void beginJob();
    void beginRun(/*art::Run&*/);
    int firstEvent;
  
    void endJob();
  
    // This is called for each event.
    //void produce(art::Event & e);
    void analyze(art::Event const& e);
  
  private:
  
    // Start: run time parameters
  
    // The module label of this module.
    std::string _moduleLabel;
  
    // Label of the module that made the hits.
    std::string _makerModuleLabel;
  
    // Label of the module that made the hits.
    std::string _timeRejecterModuleLabel;
  
    // Label of the process for the remaping of the Tracker hit by Cell/Straw ID.
    std::string _mapTrackerHitByID;
  
    // Label to the Pattern Recognition module.
    std::string _patternRecoModuleLabel;
  
    // Label of the generator.
    std::string _generatorModuleLabel;
    // Label of the module that created the data products.
    std::string _g4ModuleLabel;
    // Instance names of data products
    std::string _targetStepPoints;
    string _hitMakerModuleLabel;

    // Diagnostics level.
    int _diagLevel;
    int _debugLvl;
  
    //kalman fit track
    DchDetector* dchdet;
    std::string _materialdb;
    KalFitI _kfit;
  
    // End: run time parameters
  
    // Pointers to histograms, ntuples, TGraphs.

    TH1I*         _hClockCicles;
    TH1F*         _hExecTime;

    // Some ugly but necessary ROOT related bookkeeping:
  
    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;
  
    // Save directory from beginJob so that we can go there in endJob. See note 3.
    //    TDirectory* _directory;
    bool FillMCInfo(const art::Event& event);

    std::vector<hitIndex> strawhits;
    HelixTraj seed;
    double time0;
    double s0;
    std::vector<double> hitflt;

    CLHEP::Hep3Vector gen_pos;
    CLHEP::Hep3Vector gen_mom;

    // Pointers to histograms, ntuples, TGraphs.
    void FillHistos(KalFitResult& kres); 

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
    HelixParams startpar;
    HelixParams recopar;



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
  
  ITTrackRecoFit::ITTrackRecoFit(fhicl::ParameterSet const& pset) :
    
    // Run time parameters
    _moduleLabel(pset.get<std::string>("module_label")),/*@module_label*/
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _timeRejecterModuleLabel(pset.get<std::string>("tRejecterModuleLabel")),
    _mapTrackerHitByID(pset.get<std::string>("mapTrackerHitByID")),
    _patternRecoModuleLabel(pset.get<std::string>("patternRecoModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
    _targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget")),
    _hitMakerModuleLabel(pset.get<std::string>("hitMakerModuleLabel", "makeDcH")),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _materialdb(pset.get<std::string>("materialdb", "KalmanFitTests/test/material_dumb.txt")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    // ROOT objects that are the main focus of this example.
    _hExecTime(0),
    // Some ugly but necessary ROOT related bookkeeping.
    _application(0),
    seed(TrkParams(HelixTraj::NHLXPRM)),
    nbadfit(0),
    startpar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1)),
    recopar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1))
  {
    std::cout<<"Constructed"<<std::endl;
    firstEvent=0;
    _debugLvl=_diagLevel;

  }

  void ITTrackRecoFit::beginJob(){
    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;
    _hClockCicles       = tfs->make<TH1I>( "hClockCicles",   "N clock cicles needed to analyze one Event by ITTrackRecoFit", 2000, 5000, 15000  );
    _hExecTime          = tfs->make<TH1F>( "hExecTime",   "Execution time to analyze one Event by BkgTrackRejecterByTime", 1000, 0.0, 1.0  );

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
    treefit->Branch("startpar","HelixParams",&startpar);
    treefit->Branch("recopar","HelixParams",&recopar);
    treefit->Branch("fitinfo",&recoinfo,"fit/I:chi2/F:chi2out/F:nhits/I:nhitstot/I:nturn/I:momin/F:momout/F:fitmom/F:fitmomerr/F:t0/F:t0fit/F:errt0/F");
    treefit->Branch("eventinfo",&eventinfo,"elosstgt/F:ntgt/I");


  }

  void ITTrackRecoFit::beginRun(/*art::Run&*/){
    std::cout<<"init dchdetector"<<std::endl;
    DchGDchP gdch(_materialdb);
    dchdet=new DchDetector(gdch,true);
  }
  
  //  void ITTrackRecoFit::produce(art::Event & event ) {
  void ITTrackRecoFit::analyze(art::Event const& event ) {
    if(!firstEvent){firstEvent=1;beginRun();}

    //--------------------------------------------

    //const Tracker& tracker = getTrackerOrThrow();
    //const ITracker &itr = static_cast<const ITracker&>( tracker );
    //CellGeometryHandle *itwp = itr.getCellGeometryHandle();

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    // art::Handle<TrackerHitByID> hitByIDHandle;
    // event.getByLabel(_mapTrackerHitByID,hitByIDHandle);
    // TrackerHitByID const* hitByID = hitByIDHandle.product();
    // TrackerHitByID::const_iterator hitByID_it;
    // std::pair<TrackerHitByID::const_iterator, TrackerHitByID::const_iterator> rangeHitByID_it;
    //std::set<size_t> hitLoked;

    art::Handle<TrackSeedCollection> tracksSeedHandle;
    event.getByLabel(_patternRecoModuleLabel,tracksSeedHandle);
    TrackSeedCollection const* tracksSeed = tracksSeedHandle.product();
  
    if(!FillMCInfo(event)) return;

    //clock_t startClock = clock();

    size_t nStrawPerEvent = hits->size();
    size_t nTimeClusPerEvent = tclusts->size();
    size_t nTracksSeeds = tracksSeed->size();
    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"event "<<event.id().event()<<" tot N hit "<<nStrawPerEvent<<" N tracks seed found "<<nTracksSeeds
	     <<" N time peaks "<<nTimeClusPerEvent<<std::endl;
    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
  
    for (size_t iTrackSeed=0; iTrackSeed<nTracksSeeds; ++iTrackSeed) {
      TrackSeed const& iTrkSeed(tracksSeed->at(iTrackSeed));
    
      //Fake conversion !!! Needed because I had to redefine hitIndex as different type to avoid compilation problem!! FIXME
      //all time peaks
      std::vector<hitIndex> tpeakhits;
      const TrackerHitTimeCluster&  tclust=*(iTrkSeed._relatedTimeCluster);
      for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); 
	   iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
	size_t iglbHit = iTCHit->key();
	tpeakhits.push_back(mu2e::hitIndex(iglbHit));
      }
      //all track selected
      std::vector<hitIndex> goodhits;
      const std::vector<HitIndex> &trkseedhits = iTrkSeed._fullTrkSeed._selectedTrackerHitsIdx;
      for (std::vector<HitIndex>::const_iterator loopPoints_it = trkseedhits.begin();
	   loopPoints_it != trkseedhits.end(); ++loopPoints_it) {
	goodhits.push_back( mu2e::hitIndex((*loopPoints_it)._index) );
      }

      HelixTraj seed(TrkParams(HelixTraj::NHLXPRM));
      HelixVal2HelixTraj(iTrkSeed._fullTrkSeed,seed);
      TrkDef seeddef(hits,goodhits);
      seeddef.setHelix(seed);
      seeddef.setT0(TrkT0(iTrkSeed._t0,iTrkSeed._errt0));
      KalFitResult seedfit(seeddef);
      seeddef.setIndices(goodhits);
      _kfit.makeTrack(seedfit);
      if(seedfit._fit.success()){
	TrkHotList* hots = seedfit._krep->hotList();
	std::cout<<"active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
	std::vector<hitIndex> strawhitsinactive;
      
	_kfit.reActivateHitsbyTurn(seedfit);
      
	std::cout<<"ReTurn:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
      
	_kfit.reActivateHitsbyChi2(seedfit);
      
	std::cout<<"ReChi2:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
      
	_kfit.addHitsUnique(seedfit,tpeakhits,false);
      
	std::cout<<"ReAddT:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;

	_kfit.reActivateHitsbyChi2(seedfit);
      
	std::cout<<"ReAddChi2:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl;
      
      
      }
      FillHistos(seedfit);

      KalRep *kalrep = seedfit._krep;
      TrkErrCode err= seedfit._fit;
      if(err.success()) {
	double chi2=kalrep->chisquared(trkIn);
	double chi2out=kalrep->chisquared(trkOut);
	if(_diagLevel){
	  std::cout<<"chi2 in "<<chi2<<" out "<<chi2out<<std::endl;
	  //      std::cout<<"hello err="<<err<<std::endl;
	  seed.printAll(std::cout);
	  std::cout<<"fitted at begin of helix iters "<<seedfit._nt0iter
		   <<" t0 "<<seeddef.t0()._t0<<" +- "<<seeddef.t0()._t0err
		   <<" at end "<<kalrep->t0()._t0<<" +- "<<kalrep->t0()._t0err<<std::endl;
	}
	double fltlen=kalrep->firstHit()->globalLength();
	HelixParams recopar = kalrep->helix(fltlen);
	HepVector recoparams = recopar.params();
	HepSymMatrix recocovar = recopar.covariance();
	double Bval = _kfit.bField()->bFieldNominal();
	double momem=Constants::c*1.0e-3*Bval/recoparams[2]*sqrt(1.+recoparams[4]*recoparams[4]);
	if(_diagLevel){
	  std::cout<<"Fit mom "<<momem<<std::endl;
	  recopar.printAll(std::cout);
	}
      }
    
    }

  } // end produce


  bool ITTrackRecoFit::FillMCInfo(const art::Event& event){
    int iev=event.id().event();

    GeomHandle<ITracker> itr;
    CellGeometryHandle *itwp = itr->getCellGeometryHandle();
    //const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    startpar=HelixParams(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
    recopar=HelixParams(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
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
    startpar=HelixParams(startPar,dummy);
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

	  // if((_oneturn&&cosphidirprev<0&&cosphidir>0)){
	  //   break;
	  // }

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

  
  void ITTrackRecoFit::FillHistos(KalFitResult& kres){

    KalRep *kalrep = kres._krep;
    TrkErrCode err= kres._fit;

    recoinfo.nhits=kres._krep?kres._krep->hotList()->nActive():-1;

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
	  <<" t0 at first hit "<<time0<<" fitted at z=0 of helix niters "<<kres._nt0iter
	  <<" t0 "<<kres._tdef.t0()._t0<<" +- "<<kres._tdef.t0()._t0err
	  <<" at end "<<kalrep->t0()._t0<<" +- "<<kalrep->t0()._t0err<<endl;
    }
    recoinfo.t0fit=kalrep->t0()._t0;
    recoinfo.errt0=kalrep->t0()._t0err;

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
    HelixParams recobeam = kalrep->helix(fltlenbeam);
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
  
  void ITTrackRecoFit::endJob()
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

}  // end namespace mu2e

using mu2e::ITTrackRecoFit;
DEFINE_ART_MODULE(ITTrackRecoFit);
