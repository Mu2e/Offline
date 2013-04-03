//
// Fit of the reconstructed tracks for the ITracker
//
// $Id: ITTrackRecoFit_module.cc,v 1.4 2013/04/03 22:08:21 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/04/03 22:08:21 $
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
#include "art/Framework/Core/EDProducer.h"
//#include "art/Framework/Core/EDAnalyzer.h"
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
#include "KalmanTests/inc/KalFitResult.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkPoca.hh"
// Data Output
#include "KalmanTests/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"

#include "DchGeom/DchDetector.hh"

#include "KalmanTestsI/inc/DchGDchP.hh"
#include "KalmanTestsI/inc/KalFitI.hh"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "KalmanTestsI/inc/kalFitOutUtils.hh"

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

  class ITTrackRecoFit : public art::EDProducer {
  public:

    explicit ITTrackRecoFit(fhicl::ParameterSet const& pset);
    virtual ~ITTrackRecoFit() {
    }
  
    virtual void beginJob();
    void beginRun(/*art::Run&*/);
  
    void endJob();
  
    // This is called for each event.
    void produce(art::Event & e);
    //void analyze(art::Event const& e);
  
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
    int _t0use;
    //kalman fit track
    DchDetector* dchdet;
    std::string _materialdb;
    bool _doEndCapExtrapol;
    TrkParticle _tpart; // particle type being searched for
    TrkFitDirection _fdir;  // fit direction in search

    KalFitI _kfit;
  
    std::string _iname; // data instance name
    //
    PayloadSaver _payloadSaver;
    // End: run time parameters
  
    // Pointers to histograms, ntuples, TGraphs.

    TH1I*         _hClockCicles;
    TH1F*         _hExecTime;

    // Some ugly but necessary ROOT related bookkeeping:
  
    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;
  
    std::vector<hitIndex> strawhits;
    HelixTraj seed;

    kalFitOutUtils _histoOut;

    bool _firstEv;
    Int_t _eventid;

    bool _doMinPrints;

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
    _t0use(pset.get<int>("t0use",0)),
    _materialdb(pset.get<std::string>("materialdb", "")),
    _doEndCapExtrapol(pset.get<bool>("doEndCapExtrapol",false)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    _payloadSaver(pset),
    // ROOT objects that are the main focus of this example.
    _hExecTime(0),
    // Some ugly but necessary ROOT related bookkeeping.
    _application(nullptr),
    seed(TrkParams(HelixTraj::NHLXPRM)),
    _histoOut (_g4ModuleLabel,_generatorModuleLabel,_targetStepPoints,_hitMakerModuleLabel,
               _diagLevel,0,pset.get<string>("vdStepPoints","virtualdetector"),
               GenId(pset.get<int>("genidselect",GenId::conversionGun))),
    _firstEv(true),
    _doMinPrints(false)
  {
    std::cout<<"Constructed"<<std::endl;
    _debugLvl=_diagLevel;
    _kfit.setMaxdist4Addhit(20.);
    // tag the data product instance by the direction and particle type found by this fitter
        _iname = _fdir.name() + _tpart.name();
        produces<KalRepCollection>(_iname);
        produces<KalRepPayloadCollection>();
    // set # bins for time spectrum plot
  }

  void ITTrackRecoFit::beginJob(){
          _histoOut.bookHitos();
          _eventid = 0;
  }

  void ITTrackRecoFit::beginRun(/*art::Run&*/){
    std::cout<<"init dchdetector"<<std::endl;
    DchGDchP gdch(_materialdb,_kfit.useDetailedWrSuppDescr(),_kfit.useSingleCellDescr());
    dchdet=new DchDetector(gdch,true);
  }
  
  void ITTrackRecoFit::produce(art::Event & event ) {
          //void ITTrackRecoFit::analyze(art::Event const& event ) {
            if(_firstEv){
                    beginRun();
                    _histoOut._bfield=&_kfit.bField();
                    _firstEv=false;
                    _doMinPrints = (_debugLvl>0);
            }
            ++_eventid;

            //_histoOut.FillMCInfo(event,strawhits,seed);
            if(!_histoOut.FillMCInfo(event,strawhits,seed)) return;

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

            // create output
            unique_ptr<KalRepCollection> tracks(new KalRepCollection );

            //clock_t startClock = clock();

            size_t nStrawPerEvent = hits->size();
            size_t nTimeClusPerEvent = tclusts->size();
            size_t nTracksSeeds = tracksSeed->size();

            if (_doMinPrints || event.id().event()%100==0) {
                    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
                    std::cout<<"event "<<event.id().event()<<" tot N hit "<<nStrawPerEvent<<" N tracks seed found "<<nTracksSeeds
                                    <<" N time peaks "<<nTimeClusPerEvent<<std::endl;
                    std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
            }

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
                unsigned int i=0;
                for(i=0;i<strawhits.size();i++){
                  if(strawhits[i]._index==(*loopPoints_it)._index) {
                    goodhits.push_back( mu2e::hitIndex((*loopPoints_it)._index,strawhits[i]._ambig) );
                    break;
                  }
                }
                if(i==strawhits.size()) {
                  goodhits.push_back( mu2e::hitIndex((*loopPoints_it)._index));
                }
              }

              std::vector<std::vector<hitIndex> > agoodhits;
              for(size_t i=0;i<iTrkSeed._loopSeeds.size();i++){
                //      HelixTraj recoseed(TrkParams(HelixTraj::NHLXPRM));
                //              HelixVal2HelixTraj(iTrkSeed._loopSeeds[i],recoseed);
                agoodhits.push_back(std::vector<hitIndex>());
                for(size_t j=0;j<iTrkSeed._loopSeeds[i]._selectedTrackerHitsIdx.size();j++){
                  agoodhits[i].push_back( mu2e::hitIndex((iTrkSeed._loopSeeds[i]._selectedTrackerHitsIdx[j])._index));
                }
              }


              if (_doMinPrints) {
                      std::cout<<tclust._selectedTrackerHits.size()<<" nhits in peak"<<iTrackSeed<<" from "<<hits->size()
                                       <<" t= "<<tclust._meanTime<<" min "<<tclust._minHitTime<<" max "<<tclust._maxHitTime
                                       <<" t0mc "<<_histoOut.time0<<" t0maxmc "<<_histoOut.time0_max<<" mc at z=0 "<<_histoOut.recoinfo.t0<<std::endl;
              }

              HelixTraj recoseed(TrkParams(HelixTraj::NHLXPRM));
              HelixVal2HelixTraj(iTrkSeed._fullTrkSeed,recoseed);


              //recoseed.parameters()->covariance()*=1000;
              //recoseed.parameters()->parameter()[HelixTraj::omegaIndex]*=0.99;
              double helixzturn=fabs(2*TMath::Pi()/recoseed.omega()*recoseed.tanDip());
              double z0seed=recoseed.parameters()->parameter()[HelixTraj::z0Index];
              int nzskip=ceil((z0seed+dchdet->zlen()*0.5)/helixzturn);
              recoseed.parameters()->parameter()[HelixTraj::z0Index]=z0seed-nzskip*helixzturn;
              if (_doMinPrints) { std::cout<<"fix z0 coor "<<z0seed<<" to "<<z0seed-nzskip*helixzturn<<" zturn "<<helixzturn<<std::endl; }

              TrkDef seeddef(hits, goodhits,recoseed,_tpart,_fdir);//TrkParticle(),TrkFitDirection());
              //TrkDef seeddef(hits, agoodhits[0],recoseed,TrkParticle(),TrkFitDirection());
              //TrkDef seeddef(hits, strawhits,recoseed,TrkParticle(),TrkFitDirection());
              seeddef.setEventId(_eventid);
              seeddef.setTrackId(iTrackSeed);

              //      TrkT0 t0(iTrkSeed._t0,iTrkSeed._errt0);
              TrkT0 t0(tclust._meanTime-7.32505e+01,7.);//-4.35 forr in wire speed
              if(_t0use==1) t0=TrkT0(_histoOut.recoinfo.t0,0.0);
              //if(_t0use==2) t0=TrkT0(_histoOut.recoinfo.t0,iTrkSeed._errt0);
              if(_t0use==2) t0=TrkT0(iTrkSeed._t0,iTrkSeed._errt0);

              seeddef.setT0(t0);

              _kfit._flt0=_histoOut.s0;
              _kfit._hitflt=_histoOut.hitflt;
              KalFitResult *seedfit=new KalFitResult(seeddef);
              _kfit.makeTrack(*seedfit);

              for(int itry=0;itry<2;itry++){
                if(!seedfit->_fit.success()){
                  for(unsigned int i=0;i<iTrkSeed._loopSeeds.size();i++){
                    if(itry==0&&i==0) continue;//already fitted by main seed
                    HelixTraj recoseed(TrkParams(HelixTraj::NHLXPRM));
                    HelixVal2HelixTraj(iTrkSeed._loopSeeds[i],recoseed);
                    double helixzturn=fabs(2*TMath::Pi()/recoseed.omega()*recoseed.tanDip());
                    double z0seed=recoseed.parameters()->parameter()[HelixTraj::z0Index];
                    int nzskip=ceil((z0seed+dchdet->zlen()*0.5)/helixzturn);
                    recoseed.parameters()->parameter()[HelixTraj::z0Index]=z0seed-nzskip*helixzturn;
                    if (_doMinPrints) {
                            std::cout<<"Try to fit from loop seed "<<i<<std::endl;
                            std::cout<<"fix z0 coor "<<z0seed<<" to "<<z0seed-nzskip*helixzturn<<" zturn "<<helixzturn<<std::endl;
                    }


                    seeddef=TrkDef(hits,itry==0?goodhits:agoodhits[i],recoseed,_tpart,_fdir);//TrkParticle(),TrkFitDirection());
                    seeddef.setT0(t0);
                    delete seedfit;
                    seedfit=new KalFitResult(seeddef);

                    _kfit._flt0=_histoOut.s0;
                    _kfit._hitflt=_histoOut.hitflt;
                    _kfit.makeTrack(*seedfit);
                    if(seedfit->_fit.success()) break;
                  }
                }
                if(seedfit->_fit.success()) break;
              }

              if(seedfit->_fit.success()){

                      TrkHotList* hots = seedfit->_krep->hotList();

                      if (_doMinPrints) { std::cout<<"active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }
                      std::vector<hitIndex> strawhitsinactive;

                      _kfit.reActivateHitsbyTurn(*seedfit);
                      if(seedfit->_fit.success()){
                              if (_doMinPrints) { std::cout<<"ReTurn:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }

                              if(seedfit->_fit.success()){
                                      _kfit.reActivateHitsbyChi2(*seedfit);
                                      if (_doMinPrints) { std::cout<<"ReChi2:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }

                                      if(seedfit->_fit.success()){
                                              _kfit.addHitsUnique(*seedfit,tpeakhits,false);
                                              if (_doMinPrints) { std::cout<<"ReAddT:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }
                                              if(seedfit->_fit.success()){
                                                      _kfit.reActivateHitsbyChi2(*seedfit);
                                                      if (_doMinPrints) { std::cout<<"ReAddChi2:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }
                                              }
                                      }
                              }
                      }
              }

              if (_doEndCapExtrapol && seedfit->_fit.success()) {
                      _kfit.makeExtrapolOutTrack(*seedfit);
              }

              _histoOut.FillHistos(*seedfit,seed,iTrackSeed);

              if(seedfit->_fit.success()){
                      // save successful kalman fits in the event
                      tracks->push_back( seedfit->stealTrack() );
              } else {
                      seedfit->deleteTrack();
              }

              delete seedfit;

            }

            if(nTracksSeeds==0){
                    _histoOut.recoinfo.clearrec();
                    _histoOut.treefit->Fill();
            }

            // put the tracks into the event
            art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
            _payloadSaver.put(*tracks, tracksID, event);
            event.put(std::move(tracks),_iname);
} // end produce

  
  void ITTrackRecoFit::endJob()
  {
          _histoOut.finalizeHistos();
  }
}  // end namespace mu2e

using mu2e::ITTrackRecoFit;
DEFINE_ART_MODULE(ITTrackRecoFit);
