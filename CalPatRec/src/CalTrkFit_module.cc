///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalTrkFit_module.hh"
#include "CalPatRec/inc/Ref.hh"
#include "CalPatRec/inc/AlgorithmIDCollection.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "CalPatRec/inc/KalFitResult.hh"

// Mu2e BaBar
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

//TrkReco
#include "TrkReco/inc/TrkUtilities.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "TVector2.h"

using namespace std;
using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {
  //-----------------------------------------------------------------------------
  // comparison functor for sorting by Z(wire)
  //-----------------------------------------------------------------------------
  struct straw_zcomp : public binary_function<StrawHitIndex,StrawHitIndex,bool> {
    bool operator()(StrawHitIndex const& h1, StrawHitIndex const& h2) {

      mu2e::GeomHandle<mu2e::TTracker> handle;
      const TTracker* t = handle.get();
      const Straw* s1 = &t->getStraw(StrawIndex(h1));
      const Straw* s2 = &t->getStraw(StrawIndex(h2));

      return s1->getMidPoint().z() < s2->getMidPoint().z();
    }
  }; // a semicolumn here is required

  //-----------------------------------------------------------------------------
  // module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
  //-----------------------------------------------------------------------------
  CalTrkFit::CalTrkFit(fhicl::ParameterSet const& pset) :
    _diagLevel       (pset.get<int>                ("diagLevel"                      )),
    _debugLevel      (pset.get<int>                ("debugLevel"                     )),
    _printfreq       (pset.get<int>                ("printFrequency"                 )),
    _useAsFilter     (pset.get<int>                ("useAsFilter"                    )),    
    _addhits         (pset.get<bool>               ("addhits"                        )),
    _zsave           (pset.get<vector<double> >("ZSavePositions",
						vector<double>{-1522.0,0.0,1522.0}   )), // front, middle and back of the tracker
    _shLabel         (pset.get<string>             ("StrawHitCollectionLabel"        )),
    _shDigiLabel     (pset.get<string>             ("StrawDigiCollectionLabel"       )),
    _shpLabel        (pset.get<string>             ("StrawHitPositionCollectionLabel")),
    _shfLabel        (pset.get<string>             ("StrawHitFlagCollectionLabel"    )),
    _trkseedLabel    (pset.get<string>             ("TrackSeedModuleLabel"           )),
    _maxdtmiss       (pset.get<double>             ("MaxDtMiss"                      )),
    _maxadddoca      (pset.get<double>             ("MaxAddDoca"                     )),
    _maxaddchi       (pset.get<double>             ("MaxAddChi"                      )),
    _goodseed        (pset.get<vector<string> >("GoodKalSeedFitBits",vector<string>{})),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle"               ))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _kfit            (pset.get<fhicl::ParameterSet>("KalFitHack",fhicl::ParameterSet())),
    _kfresult        (0),
    _payloadSaver    (pset)
  {
    //    fStopwatch = new TStopwatch();

    // tag the data product instance by the direction
    // and particle type found by this fitter

    produces<KalRepCollection       >();
    produces<KalRepPtrCollection    >();
    produces<KalRepPayloadCollection>();
    produces<AlgorithmIDCollection  >();
    produces<StrawHitFlagCollection >();
    produces<KalSeedCollection      >();

    _minNMCHits = 25;

    fHackData = new THackData("HackData","Hack Data");
    gROOT->GetRootFolder()->Add(fHackData);
    //-----------------------------------------------------------------------------
    // provide for interactive disanostics
    //-----------------------------------------------------------------------------
    _ref = new Ref("CalTrkFitRef","Ref to CalTrkFit", &_kfit);

    TFolder* f_mu2e;

    f_mu2e = (TFolder*) gROOT->GetRootFolder()->FindObject("Mu2e");
    if (f_mu2e == NULL) f_mu2e = gROOT->GetRootFolder()->AddFolder("Mu2e","Mu2e Folder");

    if (f_mu2e) {
      _folder = f_mu2e->AddFolder("CalTrkFit","CalTrkFit Folder");
      _folder->Add(_ref);
    }

    fgTimeOffsets     = new SimParticleTimeOffset(pset.get<fhicl::ParameterSet>("TimeOffsets"));

    _helTraj = 0;
  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  CalTrkFit::~CalTrkFit() {
    delete _ref;
    if (_helTraj) delete _helTraj;
    //    delete fStopwatch;
  }

  //-----------------------------------------------------------------------------
  void CalTrkFit::beginJob(){

    if(_diagLevel > 0) bookHistograms();

    _eventid = 0;
  }

  //-----------------------------------------------------------------------------
  bool CalTrkFit::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    // calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();

    _kfit.setTracker(_tracker);
    _kfit.setTrackerCalib(_trackerCalib);

    return true;
  }

  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalTrkFit::bookHistograms() {
    if (_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      
      art::TFileDirectory sf_dir = tfs->mkdir("KalFit");
      _hist.nhits       = sf_dir.make<TH1F>("nhits" , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
      _hist.chi2[0]     = sf_dir.make<TH1F>("chi20"  , "chi2 distribution: all tacks", 100, 0., 10.);
      _hist.chi2[1]     = sf_dir.make<TH1F>("chi21"  , "chi2 distribution:  tacks with nhits>15", 100, 0., 10.);
      _hist.p   [0]     = sf_dir.make<TH1F>("p0"  , "p distribution: all tacks", 400, 0., 200.);
      _hist.p   [1]     = sf_dir.make<TH1F>("p1"  , "p distribution:  tacks with nhits > 15", 400, 0., 200.);
      
      _hist.kaldoca    [0]   = tfs->make<TH1F>("hkaldoca_0","doca kalfit active hits; doca [mm]", 1000, -20., 20);
      _hist.kaldoca    [1]   = tfs->make<TH1F>("hkaldoca_1","doca kalfit non active hits; doca [mm]", 1000, -20., 20);

      _hist.ntracks      [0]     = tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
      _hist.ntracks      [1]     = tfs->make<TH1F>("nseeds1"  , "number of track candidates with nhits > 15 ", 21, -0.5, 20.5);
    }
  }

  //-----------------------------------------------------------------------------
  // find the input data objects
  //-----------------------------------------------------------------------------
  bool CalTrkFit::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalTrkFit::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalTrkFit::findData: StrawHitPositionCollection with label=%s not found.\n",
             _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalTrkFit::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }

    //-----------------------------------------------------------------------------
    // find list of MC hits - for debugging only
    //-----------------------------------------------------------------------------
    art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_shDigiLabel,mcptrHandle);
    if (mcptrHandle.isValid()) {
      _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      _listOfMCStrawHits = NULL;
    }

    art::Handle<KalSeedCollection> trkseedsHandle;
    if (evt.getByLabel(_trkseedLabel, trkseedsHandle)){
      _trkseeds = trkseedsHandle.product();
    }else {
      _trkseeds = 0;
    }
    
    //-----------------------------------------------------------------------------
    // done
    //-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_trkseeds != 0);
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  bool CalTrkFit::filter(art::Event& event ) {
    const char*               oname = "CalTrkFit::produce";
    bool                      findkal (false);
    int                       nhits;
    int                       ntrkcandidates;
    int                       gen_index, sim_id;

    ::KalRep*                 krep;

    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    fgTimeOffsets->updateMap(event);

    _ntracks[0]       = 0;
    _ntracks[1]       = 0;
    _nhits_from_gen   = 0;

    // reset the fit iteration counter
    _kfit.setNIter(0);
    //     t1 = fStopwatch->RealTime();
    //     fStopwatch->Continue();
    // event printout
    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);


    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n",oname);
      return false;
    }

    unique_ptr<KalRepCollection>       tracks   (new KalRepCollection                );
    unique_ptr<KalRepPtrCollection>    trackPtrs(new KalRepPtrCollection             );
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection           );
    unique_ptr<StrawHitFlagCollection> shfcol   (new StrawHitFlagCollection(*_shfcol));
    unique_ptr<KalSeedCollection>      kscol    (new KalSeedCollection()             );

    // find the data
 
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event));

    double pEntrance(.0), step_time(-9999.);
    double time_threshold(500.);

    const KalSeed*   trkSeed(0);
    
 
    //-----------------------------------------------------------------------------
    // count the number of MC straw hits generated by the CE
    //-----------------------------------------------------------------------------
    nhits = _shcol->size();

    for (int i=0; i<nhits; i++) {
      mu2e::PtrStepPointMCVector const& mcptr( _listOfMCStrawHits->at(i ) );
      const mu2e::StepPointMC* Step = mcptr[0].operator ->();

      gen_index = -1;
      if (Step) {
        art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle();

        if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
        else                         gen_index = -1;

        sim_id        = simptr->id().asInt();
      }
      if (gen_index >0 && sim_id == 1) {
        step_time = fgTimeOffsets->timeWithOffsetsApplied(*Step);
        step_time = fmod(step_time,_mbtime);
        if (step_time > time_threshold) {
          ++_nhits_from_gen;
          if (Step->momentum().mag() > pEntrance) {
            pEntrance = Step->momentum().mag();
          }
        }
      }
    }
    if (pEntrance < 100. ) _nhits_from_gen = 0;

    _kfit.setStepPointMCVectorCollection(_listOfMCStrawHits);

    //-----------------------------------------------------------------------------
    // loop over found time peaks - for us, - "eligible" calorimeter clusters
    //-----------------------------------------------------------------------------
    ntrkcandidates = _trkseeds->size();
    
    for (int ipeak=0; ipeak<ntrkcandidates; ipeak++) {
      trkSeed = &_trkseeds->at(ipeak);

      //create the list of std::vector<StrawHitIndex> from  TrkStrawHitSeed
      _nindex = trkSeed->hits().size();
      const mu2e::TrkStrawHitSeed*           seedStraw;
      std::vector<StrawHitIndex> hits_kf;
      for (int i=0; i< _nindex; ++i){
	seedStraw  = &trkSeed->hits().at(i);
	hits_kf.push_back(seedStraw->index());
      }
      
      TrkDefHack             kaldef (_shcol, hits_kf, _tpart, _fdir);
      
      // track fitting objects for this track candidate
      init(_kfresult, &kaldef );

// find the segment at the 0 flight
      double flt0 = trkSeed->flt0();
      auto kseg = trkSeed->segments().end();
      for(auto iseg= trkSeed->segments().begin(); iseg != trkSeed->segments().end(); ++iseg){
	if(iseg->fmin() <= flt0 && iseg->fmax() > flt0){
	  kseg = iseg;
	  break;
	}
      }
      if(kseg == trkSeed->segments().end()){
	std::cout << "Helix segment range doesn't cover flt0" << std::endl;
	kseg = trkSeed->segments().begin();
      }
      // create a trajectory from the seed. This shoudl be a general utility function that
      // can work with multi-segment seeds FIXME!
      // create CLHEP objects from seed native members.  This will
      // go away when we switch to SMatrix FIXME!!!
      HepVector    pvec(5,0);
      HepSymMatrix pcov(5,0);
      kseg->helix().hepVector(pvec);
      kseg->covar().symMatrix(pcov);
      // Create the traj from these
      if (_helTraj == 0)  _helTraj = new HelixTraj(pvec, pcov);
      else               *_helTraj = HelixTraj(pvec, pcov);
      kaldef.setHelix(*_helTraj);
      
      //--------------------------------------------------------------------------------
      // 2015-03-23 G. Pezzu: fill info about the doca
      //--------------------------------------------------------------------------------
      if (trkSeed->status().hasAllProperties(_goodseed)) {
	// check the seed has the same basic parameters as this module expects
	if(trkSeed->particle() != _tpart || trkSeed->fitDirection() != _fdir ) {
	  throw cet::exception("RECO")<<"mu2e::KalFinalFit: wrong particle or direction"<< endl;
	}
	// seed should have at least 1 segment
	if(trkSeed->segments().size() < 1){
	  throw cet::exception("RECO")<<"mu2e::KalFinalFit: no segments"<< endl;
	}

	if (_debugLevel > 0) printf("CalTrkFit::produce] calling _kfit.makeTrack\n");

	_kfit.makeTrack(*_kfresult, trkSeed->caloCluster().get());

	if (_debugLevel > 0) {
	  printf("[CalTrkFit::produce] kalfit status = %i\n", _kfresult->_fit.success());
	  _kfit.printHits(*_kfresult,"CalTrkFit::produce kalfit_001");
	}

	if (_kfresult->_fit.success()) {
	  findkal = true;

	  if (_addhits) {
	    //-----------------------------------------------------------------------------
	    // this is the default. First, add back the hits on this track
	    // if successfull, try to add missing hits, at this point external errors were
	    // set to zero
	    // assume this is the last iteration
	    //-----------------------------------------------------------------------------
	    //	      int last_iteration = _kfit.maxIteration();
	    int last_iteration = -1;
	      
	    _kfit.unweedHits(*_kfresult,_maxaddchi);
	    if (_debugLevel > 0) _kfit.printHits(*_kfresult,"CalTrkFit::produce after unweedHits");
	      
	    std::vector<StrawHitIndex> misshits;
	    findMissingHits(*_kfresult,misshits);
	    //-----------------------------------------------------------------------------
	    // if new hits have been added, add then and refit the track.
	    // Otherwise - just refit the track one last time
	    // in both cases
	    //-----------------------------------------------------------------------------
	    if (misshits.size() > 0) {
	      _kfit.addHits(*_kfresult,_shcol,misshits, _maxaddchi, trkSeed->caloCluster().get());
	    }
	    else {
	      _kfit.fitIteration(*_kfresult,last_iteration, trkSeed->caloCluster().get());
	    }

	    if (_debugLevel > 0) _kfit.printHits(*_kfresult,"CalTrkFit::produce after addHits");
	    //-----------------------------------------------------------------------------
	    // and weed hits again to insure that addHits doesn't add junk
	    //-----------------------------------------------------------------------------
	    _kfit.weedHits(*_kfresult,last_iteration);
	  }
	  //-----------------------------------------------------------------------------
	  // now evaluate the T0 and its error using the straw hits
	  //-----------------------------------------------------------------------------
	  _kfit.updateT0(*_kfresult);
	    
	  if (_debugLevel > 0) {
	    _kfit.printHits(*_kfresult,"CalTrkFit::produce : final, after weedHits");
	  }
	  //-----------------------------------------------------------------------------
	  // done, fill debug histograms
	  //-----------------------------------------------------------------------------
	  if (_diagLevel > 0) {
	    TrkHitVector const& hot_l = _kfresult->_krep->hitVector();
	    const mu2e::TrkStrawHit* hit;
	    int                      hit_index;
	    Hep3Vector               tdir;
	    HepPoint                 tpos;

	    _kfresult->_krep->traj().getInfo(0.0,tpos,tdir);

	    for (int i=0; i< _nindex; ++i){
	      int               hIndex= hits_kf[i];
	      StrawHit const*   sh    = & _shcol->at(hIndex);
	      Straw const&      straw = _tracker->getStraw(sh->strawIndex());
	      CLHEP::Hep3Vector hpos  = straw.getMidPoint();
	      CLHEP::Hep3Vector hdir  = straw.getDirection();
	      bool              found = false;

	      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	      HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
	      TrkLineTraj       htraj(spt,hdir,-20,20);
	      // estimate flightlength along track.  This assumes a constant BField!!!
	      double           fltlen = (hpos.z()-tpos.z())/tdir.z();
	      TrkPoca          hitpoca(_kfresult->_krep->traj(),fltlen,htraj,0.0);

	      double           doca   = hitpoca.doca();
	      for(auto it=hot_l.begin(); it<hot_l.end(); it++) {
		hit = static_cast<const mu2e::TrkStrawHit*> (*it);
		if (!hit->isActive()) continue;
		hit_index = hit->index();
		if (int(hits_kf[i]) == hit_index){
		  found = true;
		  break;
		}
	      }

	  
	      if (found) _hist.kaldoca[0]->Fill(doca);
	      else       _hist.kaldoca[1]->Fill(doca);
	    
	    }
	  }
	}
      }
	  
      

      if (_kfresult->_fit.success()) {
	// save successful kalman fits in the event.
	// start from _kfresult, as stealTrack clears the hit pointers
	// _kfresult doesn't own anything
	//-----------------------------------------------------------------------------
        krep = _kfresult->stealTrack();

	// flg all hits as belonging to this track
	// // this is probably obsolete, FIXME!
	if(size_t(ipeak)<StrawHitFlag::_maxTrkId){
	  for(auto ihit=krep->hitVector().begin();ihit != krep->hitVector().end();++ihit){
	    if((*ihit)->isActive())shfcol->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::trackBit(ipeak));
	  }
	}
	

        tracks->push_back(krep);
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	//        tp->SetCprIndex(tracks->size());

	int best = AlgorithmID::CalPatRecBit;
	int mask = 1 << AlgorithmID::CalPatRecBit;

	algs->push_back(AlgorithmID(best,mask));

	// convert successful fits into 'seeds' for persistence.  Start with the input
	KalSeed fseed(*trkSeed);
	// reference the seed fit in this fit
	art::Handle<KalSeedCollection> ksH;
	event.getByLabel(_trkseedLabel, ksH);
	fseed._kal = art::Ptr<KalSeed>(ksH, ipeak);
	// fill other information
	fseed._t0 = krep->t0();
	fseed._flt0 = krep->flt0();
	fseed._status.merge(TrkFitFlag::kalmanOK);
	if(krep->fitStatus().success()==1) fseed._status.merge(TrkFitFlag::kalmanConverged);
	TrkUtilities::fillHitSeeds(krep, fseed._hits);
	// sample the fit at the requested z positions.  This should
	// be in terms of known positions (front of tracker, ...) FIXME!
	for(auto zpos : _zsave) {
	  // compute the flightlength for this z
	  // double fltlen = TrkUtilities::zFlight(krep->pieceTraj(),zpos);
	  double fltlen = krep->pieceTraj().zFlight(zpos);
	  // sample the momentum at this flight.  This belongs in a separate utility FIXME
	  BbrVectorErr momerr = krep->momentumErr(fltlen);
	    // sample the helix
	  double locflt(0.0);
	  const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(fltlen,locflt));
	  // fill the segment
	  KalSegment kseg;
	  TrkUtilities::fillSegment(*htraj,momerr,kseg);
	  fseed._segments.push_back(kseg);
	}
	kscol->push_back(fseed);
      }
      else {
	//-----------------------------------------------------------------------------
	// fit failed, just delete the track
	//-----------------------------------------------------------------------------
        _kfresult->deleteTrack();
      }
      //-----------------------------------------------------------------------------
      // cleanup the seed fit - why it is not being done ?
      //-----------------------------------------------------------------------------
      //       _sfresult->deleteTrack();

      if (_debugLevel > 0) {
	if (_nhits_from_gen >= _minNMCHits) {
	  //	  if (tp->_tmin > 400.){
	    if (!findkal) {
	      printf("[CalTrkFit::produce] LOOK AT: more than 25 MC hits and findKal not converged! event = %i\n", _iev);
	    }
	    //	  }
	}
      }
    }

    if (_diagLevel > 0){
      _hist.ntracks[0]->Fill(_ntracks[0]);
      _hist.ntracks[1]->Fill(_ntracks[1]);
    }
    //-----------------------------------------------------------------------------
    // put reconstructed tracks into the event record
    //-----------------------------------------------------------------------------
    //  END:;
    int     ntracks = tracks->size();

    art::ProductID krcolID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*tracks, krcolID, event);
    event.put(std::move(tracks)   );
    event.put(std::move(trackPtrs));
    event.put(std::move(algs     ));
    event.put(std::move(shfcol   ));
    event.put(std::move(kscol    ));

    if (_useAsFilter == 1) {
      if (ntracks > 0){
	return true;
      } else{
	return false;
      }
    } else {
      return true;
    }
  }

  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalTrkFit::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }



  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalTrkFit::findMissingHits(KalFitResult& kalfit,std::vector<StrawHitIndex>& misshits) {
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    int        radius_ok;
    double     dt;

    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    if (_debugLevel > 0) {
      printf("[CalTrkFit::findMissingHits]      shId    sec     panel       doca   \n");
    }

    for (unsigned istr=0; istr<nstrs;++istr) {
      //----------------------------------------------------------------------
      // 2015-02-11 gianipez and P. Murat changed the selection bit
      //            for searching for missed hits
      //----------------------------------------------------------------------
      StrawHit const& sh = _shcol->at(istr);
      //-----------------------------------------------------------------------------
      // I think, we want to check the radial bit: if it is set, than at least one of
      // the two measured times is wrong...
      //-----------------------------------------------------------------------------
      radius_ok = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      dt        = _shcol->at(istr).time()-kalfit._krep->t0()._t0;

      if (radius_ok && (fabs(dt) < _maxdtmiss)) {
        // make sure we haven't already used this hit
        std::vector<TrkStrawHit*>::iterator ifnd = find_if(kalfit._hits.begin(),kalfit._hits.end(),FindTrkStrawHit(sh));
        if(ifnd == kalfit._hits.end()){
          // good in-time hit.  Compute DOCA of the wire to the trajectory
          Straw const& straw = _tracker->getStraw(sh.strawIndex());
          CLHEP::Hep3Vector hpos = straw.getMidPoint();
          CLHEP::Hep3Vector hdir = straw.getDirection();
          // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
          HepPoint spt(hpos.x(),hpos.y(),hpos.z());
          TrkLineTraj htraj(spt,hdir,-20,20);
          // estimate flightlength along track.  This assumes a constant BField!!!
          double fltlen = (hpos.z()-tpos.z())/tdir.z();
          TrkPoca hitpoca(kalfit._krep->pieceTraj(),fltlen,htraj,0.0);

          if (_debugLevel > 0) {
            printf("[CalTrkFit::findMissingHits] %8i  %6i  %8i  %10.3f \n",
                   straw.index().asInt(),
                   straw.id().getPlane(),
                   straw.id().getPanel(),
                   hitpoca.doca());
          }
	  //-----------------------------------------------------------------------------
	  // flag hits with small residuals
	  //-----------------------------------------------------------------------------
          if (fabs(hitpoca.doca()) < _maxadddoca) {
            misshits.push_back(istr);
          }
        }
      }
    }
  }


  //-----------------------------------------------------------------------------
  void CalTrkFit::init(KalFitResult*& KRes, TrkDefHack* TDef) {

    if (KRes != 0) {
      KRes->deleteTrack();
      delete KRes;
    }
    KRes = new KalFitResult();
    KRes->_tdef = TDef;
  }
}

using mu2e::CalTrkFit;
DEFINE_ART_MODULE(CalTrkFit);
