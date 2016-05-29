//
// TTracker Kalman Fit launcher
//
// $Id: TrkRecFit_module.cc,v 1.4 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
// Mu2e
#include "TrkReco/inc/TrkDef.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/KalFit.hh"
#include "TrkDiag/inc/KalDiag.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std; 
using CLHEP::Hep3Vector;

namespace mu2e 
{
  class TrkRecFit : public art::EDProducer
  {
    public:
      explicit TrkRecFit(fhicl::ParameterSet const&);
      virtual ~TrkRecFit();
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
      void endJob();
    private:
      unsigned _iev;
      // MC tools
      KalDiag _kdiag;
      // configuration parameters
      int _diag,_debug;
      int _printfreq;
      bool _addhits; 
      art::Handle<mu2e::StrawHitCollection> _strawhitsH;
      art::Handle<TrackSeedCollection> _trksSeedH;
      // event object labels
      string _shLabel;
      string _shpLabel;
      string _shfLabel;
      string _seedFinderLabel; // Label to the Pattern Recognition module.
      StrawHitFlag _addsel;
      StrawHitFlag _addbkg;
      double _maxdtmiss;
      // outlier cuts
      double _maxseeddoca,_maxhelixdoca,_maxadddoca, _maxaddchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      double _helicity; // cache the value of helicity, which can be computed from the above with the BField
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      StrawHitFlagCollection* _flags;
      const StrawHitPositionCollection* _shpcol;
      const TrackSeedCollection * _tscol;
      // Kalman fitters.  Seed fit has a special configuration
      KalFit _seedfit, _kfit;
      string _iname; // data instance name
      //
      PayloadSaver _payloadSaver;
      // helper functions
      bool findData(const art::Event& e);
      void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec);
      void findMissingHits(KalRep* krep, vector<hitIndex>& indices);
      void createDiagnostics();
      void fillFitDiag(TrackSeed const& seed, TrkDef const& tdef,const KalRep* seedrep,const KalRep* krep);
      void HelixVal2HelixTraj (const HelixVal &helIn, HelixTraj &helOut);

      // fit tuple variables
      Int_t _eventid;
      Int_t _nadd;
      Int_t _helixfail,_seedfail,_kalfail;
      helixpar _hpar,_spar;
      helixpar _sparerr;
      Int_t _snhits, _snactive, _sniter, _sndof;
      Float_t _schisq, _st0;
      Int_t _nchit;
      Int_t _npeak;
      Float_t _tpeak;
      // hit filtering tuple variables
      vector<TrkHitFilter> _sfilt, _hfilt;
      // flow diagnostic
      TH1F* _cutflow, *_ccutflow;
      int _icepeak;
  };

  TrkRecFit::TrkRecFit(fhicl::ParameterSet const& pset) :
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel(pset.get<string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _seedFinderLabel(pset.get<string>("SeedCollectionLabel","RobustHelixFinder")),
    _addsel(pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _addbkg(pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _maxseeddoca(pset.get<double>("MaxSeedDoca",10.0)),
    _maxhelixdoca(pset.get<double>("MaxHelixDoca",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _helicity(0.0),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet()),_fdir),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet()),_fdir),
    _payloadSaver(pset)
  {
// tag the data product instance by the direction and particle type found by this fitter
    _iname = _fdir.name() + _tpart.name();
    produces<KalRepCollection>(_iname);
    produces<KalRepPtrCollection>(_iname);
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>(_iname);
  }

  TrkRecFit::~TrkRecFit(){}

  void TrkRecFit::beginJob(){
    // create diagnostics if requested
    if(_diag > 0)createDiagnostics();
    // create a histogram of throughput: this is a basic diagnostic that should ALWAYS be on
    art::ServiceHandle<art::TFileService> tfs;
    _cutflow=tfs->make<TH1F>("cutflow","Cutflow",10,-0.5,9.5);
    _cutflow->GetXaxis()->SetBinLabel(1,"All Events");
    _cutflow->GetXaxis()->SetBinLabel(2,"Time Peak");
    _cutflow->GetXaxis()->SetBinLabel(3,"Helix Fit");
    _cutflow->GetXaxis()->SetBinLabel(4,"Seed Fit");
    _cutflow->GetXaxis()->SetBinLabel(5,"Kalman Fit");

    if(_diag>1){
      _ccutflow=tfs->make<TH1F>("ccutflow","CE Cutflow",10,-0.5,9.5);
      _ccutflow->GetXaxis()->SetBinLabel(1,"All Events");
      _ccutflow->GetXaxis()->SetBinLabel(2,"CE hits tracker");
      _ccutflow->GetXaxis()->SetBinLabel(3,"CE hits in time window");
      _ccutflow->GetXaxis()->SetBinLabel(4,"CE time peak");
      _ccutflow->GetXaxis()->SetBinLabel(5,"CE Helix NHits");
      _ccutflow->GetXaxis()->SetBinLabel(6,"CE Helix Init");
      _ccutflow->GetXaxis()->SetBinLabel(7,"CE Helix XY Fit");
      _ccutflow->GetXaxis()->SetBinLabel(8,"CE Helix #phiZ Fit");
      _ccutflow->GetXaxis()->SetBinLabel(9,"CE Seed Fit");
      _ccutflow->GetXaxis()->SetBinLabel(10,"CE Kalman Fit");
    }
    _eventid = 0;
  }

  void TrkRecFit::beginRun(art::Run& ){
 // calculate the helicity
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    // positive helicity is clockwise rotation around the direction of axial motion (negative dphi/dz)  
    _helicity = copysign(1.0,-_fdir.dzdt()*_tpart.charge()*field.z());
  }

  void TrkRecFit::produce(art::Event& event ) {
    _eventid = event.event();
    _cutflow->Fill(0.0);
    if(_diag>1)_ccutflow->Fill(0.0);
    // create output
    unique_ptr<KalRepCollection>    tracks(new KalRepCollection );
    unique_ptr<KalRepPtrCollection> trackPtrs(new KalRepPtrCollection );
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event,_iname));
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"TrkRecFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkRecFit: data missing or incomplete"<< endl;
    }
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    unique_ptr<StrawHitFlagCollection> flags(_flags );
    // find mc truth if we're making diagnostics
    if(_diag > 0 && !_kdiag.findMCData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkRecFit: MC data missing or incomplete"<< endl;
      _nchit = _kdiag.countCEHits();
      if(_nchit>14)_ccutflow->Fill(1.0);
    }

    // find the time peaks in the time spectrum of selected hits.  Otherwise, take all
    // selected hits as a peak
    _icepeak = -1;

    if(_diag>0){
// fill primary particle MC truth information
      _kdiag.kalDiag(0,false);
    }
    // dummy objects
    static TrkDef dummydef;
    static TrackSeed dummyseed;
    static KalRep* dummyrep(0);
    // loop over the accepted time peaks
    if(_tscol->size()>0)_cutflow->Fill(1.0);
    if(_diag>1 && _icepeak >=0)_ccutflow->Fill(3.0);
    bool findhelix(false), findseed(false), findkal(false);

    for (size_t iTrackSeed=0; iTrackSeed<_tscol->size(); ++iTrackSeed) {
      TrackSeed const& iTrkSeed(_tscol->at(iTrackSeed));
      unsigned ipeak = iTrackSeed;
      findhelix = true;
      //all track selected
      std::vector<hitIndex> goodhits;

      const std::vector<hitIndex> &trkseedhits = iTrkSeed._timeCluster._strawHitIdxs;
      for (std::vector<hitIndex>::const_iterator loopPoints_it = trkseedhits.begin();
        loopPoints_it != trkseedhits.end(); ++loopPoints_it) {
        goodhits.push_back( mu2e::hitIndex(loopPoints_it->_index,loopPoints_it->_ambig) );
      }

      HelixTraj recoseed(TrkParams(HelixTraj::NHLXPRM));
      HelixVal2HelixTraj(iTrkSeed._helix,recoseed);
      if(_debug>1) {
// 	std::cout<<"Seed parameters:"<<std::endl<<iTrkSeed;
// 	std::cout<<"Converted into HelixTraj:"<<std::endl;
	recoseed.printAll(std::cout);
      }

      TrkDef seeddef(_shcol,goodhits,recoseed,_tpart,_fdir);
      TrkT0 t0(iTrkSeed._timeCluster._t0,iTrkSeed._timeCluster._errt0);
      seeddef.setT0(t0);
      TrkDef kaldef(seeddef);
      // initialize filters.  These are used only for diagnostics
      _hfilt.clear();
      _sfilt.clear();

      filterOutliers(seeddef,seeddef.helix(),_maxhelixdoca,_hfilt);
      // now, fit the seed helix from the filtered hits
      KalRep *seedrep(0), *krep(0);
      _seedfit.makeTrack(seeddef,seedrep);
      if(seedrep != 0 && seedrep->fitStatus().success()){
        findseed = true;
        // find the helix parameters from the helix fit, and initialize the full Kalman fit with this
        double locflt;
        const HelixTraj* shelix = dynamic_cast<const HelixTraj*>(seedrep->localTrajectory(seedrep->flt0(),locflt));
        kaldef.setHelix(*shelix);
        // filter the outliers
        filterOutliers(kaldef,seedrep->traj(),_maxseeddoca,_sfilt);
        _kfit.makeTrack(kaldef,krep);
        // if successfull, try to add missing hits
        if(krep != 0 && krep->fitStatus().success()){
          findkal = true;
          if(_addhits){
            // first, add back the hits on this track
            _kfit.unweedHits(krep,_maxaddchi);
            vector<hitIndex> misshits;
            findMissingHits(krep,misshits);
            if(misshits.size() > 0){
              _kfit.addHits(krep,_shcol,misshits,_maxaddchi);
            }
          }
        }
      }
      // fill fit diagnostics if requested
      if(_diag > 0)
	fillFitDiag(iTrkSeed,kaldef,seedrep,krep);
      if(_diag > 1 && (int)ipeak == _icepeak){
	if(seedrep != 0 && seedrep->fitStatus().success())_ccutflow->Fill(8.0);
	if(krep != 0 && krep->fitStatus().success())_ccutflow->Fill(9.0);
      }
      if(krep != 0 && krep->fitStatus().success()){	// flag the hits used in this track. 
	if(!krep->fitCurrent()){
	  cout << "Fit not current! " << endl;
	}
	  // flag the hits used in this track. 
	if(ipeak<16){
	  for(auto ihit=krep->hitVector().begin();ihit != krep->hitVector().end();++ihit){
	    if((*ihit)->isActive())_flags->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::trackBit(ipeak));
	  }
	}
	// save successful kalman fits in the event
	tracks->push_back(krep);
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
      } else
	delete krep;
      // cleanup the seed fit
      delete seedrep;
    }
    // cutflow diagnostics
    if(findhelix)_cutflow->Fill(2.0);
    if(findseed)_cutflow->Fill(3.0);
    if(findkal)_cutflow->Fill(4.0);
    // add a dummy entry in case there are no helices
    if(_diag > 0 && _tscol->size() == 0)
      fillFitDiag(dummyseed,dummydef,dummyrep,dummyrep);
    // put the tracks into the event
    art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*tracks, tracksID, event);
    event.put(move(tracks),_iname);
    event.put(move(trackPtrs),_iname);
    event.put(move(flags),_iname);
  }

  void TrkRecFit::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

  // find the input data objects 
  bool TrkRecFit::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _shpcol = 0;
    _tscol = 0;

    if(evt.getByLabel(_shLabel,_strawhitsH))
      _shcol = _strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
    if (evt.getByLabel(_seedFinderLabel,_trksSeedH))
      _tscol = _trksSeedH.product();
// don't require stereo hits: they are only used for diagnostics
    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && _tscol != 0;
  }

  void TrkRecFit::filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec){
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    traj.getInfo(0.0,tpos,tdir);
    // tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const vector<hitIndex>& indices = mytrk.strawHitIndices();
    vector<hitIndex> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = hits->at(indices[ihit]._index);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-20,20);
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tpos.z())/tdir.z();
      TrkPoca hitpoca(traj,fltlen,htraj,0.0);
      // flag hits with small residuals
      if(fabs(hitpoca.doca()) < maxdoca){
        goodhits.push_back(indices[ihit]);
      }
      // optional diagnostics
      if(_diag > 0){
        // summarize the MC truth for this strawhit
        TrkHitFilter thfilter;
        HepPoint tpos =  traj.position(hitpoca.flt1());
        thfilter._pos = CLHEP::Hep3Vector(tpos.x(),tpos.y(),tpos.z());
        thfilter._doca = hitpoca.doca();
        if(_kdiag.mcData()._mcdigis != 0){
          StrawDigiMC const& mcdigi = _kdiag.mcData()._mcdigis->at(indices[ihit]._index);
          // use TDC channel 0 to define the MC match
          StrawDigi::TDCChannel itdc = StrawDigi::zero;
          if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
          art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
          art::Ptr<SimParticle> const& spp = spmcp->simParticle();
          thfilter._mcpdg = spp->pdgId();
          thfilter._mcproc = spp->creationCode();
          thfilter._mcgen = -1;
          if(spp->genParticle().isNonnull())
            thfilter._mcgen = spp->genParticle()->generatorId().id();
        }
        thfvec.push_back(thfilter);
      }
    }
    // update track
    mytrk.setIndices(goodhits);
  }

  void TrkRecFit::findMissingHits(KalRep* krep,vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_addsel)&& !_flags->at(istr).hasAnyProperty(_addbkg)){
        StrawHit const& sh = _shcol->at(istr);
        if(fabs(_shcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
          // make sure we haven't already used this hit

          vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
          if(ifnd == tshv.end()){
            // good in-time hit.  Compute DOCA of the wire to the trajectory
            Straw const& straw = tracker.getStraw(sh.strawIndex());
            CLHEP::Hep3Vector hpos = straw.getMidPoint();
            CLHEP::Hep3Vector hdir = straw.getDirection();
            // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
            HepPoint spt(hpos.x(),hpos.y(),hpos.z());
            TrkLineTraj htraj(spt,hdir,-20,20);
            // estimate flightlength along track.  This assumes a constant BField!!!
            double fltlen = (hpos.z()-tpos.z())/tdir.z();
            TrkPoca hitpoca(krep->pieceTraj(),fltlen,htraj,0.0);
            // flag hits with small residuals
            if(fabs(hitpoca.doca()) < _maxadddoca){
              misshits.push_back(istr);
            }
          }
        }
      }
    }
  }

  void TrkRecFit::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
    // extend the KalDiag track diagnostic tuple
    TTree* trkdiag = _kdiag.createTrkDiag();
    trkdiag->Branch("eventid",&_eventid,"eventid/I");
    trkdiag->Branch("nadd",&_nadd,"nadd/I");
    trkdiag->Branch("helixfail",&_helixfail,"helixfail/I");
    trkdiag->Branch("seedfail",&_seedfail,"seedfail/I");
    trkdiag->Branch("kalfail",&_kalfail,"kalfail/I");
    trkdiag->Branch("hpar",&_hpar,"hd0/F:hp0/F:hom/F:hz0/F:htd/F");
    trkdiag->Branch("spar",&_spar,"sd0/F:sp0/F:som/F:sz0/F:std/F");
    trkdiag->Branch("serr",&_sparerr,"sd0err/F:sp0err/F:somerr/F:sz0err/F:stderr/F");
    trkdiag->Branch("st0",&_st0,"st0/F");
    trkdiag->Branch("snhits",&_snhits,"snhits/I");
    trkdiag->Branch("sndof",&_sndof,"sndof/I");
    trkdiag->Branch("sniter",&_sniter,"sniter/I");
    trkdiag->Branch("snactive",&_snactive,"snactive/I");
    trkdiag->Branch("schisq",&_schisq,"schisq/F");
    trkdiag->Branch("nchit",&_nchit,"nchit/I");
    trkdiag->Branch("npeak",&_npeak,"npeak/I");
    trkdiag->Branch("tpeak",&_tpeak,"tpeak/F");
    trkdiag->Branch("seedfilt",&_sfilt);
    trkdiag->Branch("helixfilt",&_hfilt);
  }

  void TrkRecFit::fillFitDiag(TrackSeed const& seed,TrkDef const& tdef,
  const KalRep* seedrep, const KalRep* krep) {
    // Peak information for this seed
    TimeCluster const*  tclust = &seed._timeCluster;
    if(tclust != 0){
      // time peak information
      _tpeak   = tclust->_t0;//tclust->_meanTime;
      _npeak   = tclust->_strawHitIdxs.size();
      _helixfail = 0;
    } else {
      _tpeak = -1.0;
      _npeak = -1;
      _helixfail = -1;
    }
    // helix information
    _hpar = helixpar(seed._helix);
    // seed fit information
    if(seedrep != 0 && seedrep->fitStatus().success()){
      _seedfail = seedrep->fitStatus().failure();
      _snhits = seedrep->nHits(); 
      _snactive = seedrep->nActive();
      _sniter = seedrep->iterations();
      _sndof = seedrep->nDof();
      _schisq = seedrep->chisq();
      _st0 = seedrep->t0()._t0;
      double loclen;
      const TrkSimpTraj* ltraj = seedrep->localTrajectory(0.0,loclen);
      _spar = helixpar(ltraj->parameters()->parameter());
      _sparerr = helixpar(ltraj->parameters()->covariance());
    } else {
      _seedfail = -1;
      _snhits = -1;
      _snactive = -1;
      _sniter = -1;
      _sndof = -1;
      _schisq = -1.0;
      _st0 = -1.0;
    }
    // count # of added hits
    _nadd = 0;
    TrkStrawHitVector tshv;
    // fill kalman fit info.  This needs to be last, as it calls TTree::Fill().
    //
    if(krep != 0) {
      _kalfail = krep->fitStatus().failure();
      _kdiag.kalDiag(krep);
      convert(krep->hitVector(),tshv);
      for(auto ihit=tshv.begin();ihit != tshv.end();++ihit){
	if((*ihit)->hitFlag()==TrkStrawHit::addedHit)++_nadd;
      }
    } else {
      _kalfail = -1;
    }
  }

  void TrkRecFit::HelixVal2HelixTraj (const HelixVal &helIn, HelixTraj &helOut) {
    CLHEP::HepVector helParams(5);
    // compare the input with this configuration's helicity: these must have the same sign
    if(_helicity*helIn.helicity() < 0.0){
      throw cet::exception("RECO")<<"mu2e::TrkRecFit: Input helix values have wrong helicity"<< endl;
    }
    // Resolve the physical ambiguity on the product of charge and direction
    if(_fdir.dzdt()*helIn._tanDip < 0.0 ){ // sign of dz/dt and tanDip disagree: invert
      // invert both charge and direction-dependent parameters (= same helicity)
      // d0 changes sign as the angular momentum changes sign
      helParams(1) = -helIn._d0;
      // phi0 flips by pi, but keep it in the range -pi -> pi
      helParams(2) = helIn._phi0 + M_PI;
      static double twopi=2*M_PI;
      if(helParams(2) > twopi)helParams(2) -= twopi;
      // omega changes sign as the charge changes sign
      helParams(3) = -helIn._omega;
      // z0 is unchanged
      helParams(4) = helIn._z0;
      // tandip changes sign as the axial velocity changes sign
      helParams(5) = -helIn._tanDip;
    } else {
      // no inversion necessary
      helParams(1) = helIn._d0;
      helParams(2) = helIn._phi0;
      helParams(3) = helIn._omega;
      helParams(4) = helIn._z0;
      helParams(5) = helIn._tanDip;
    }
    CLHEP::HepSymMatrix conv(5,1);

    HelixTraj tmpHelix(helParams,conv);
    helOut=tmpHelix;
  }

}
using mu2e::TrkRecFit;
DEFINE_ART_MODULE(TrkRecFit);
