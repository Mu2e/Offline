///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalSeedFit_module.hh"
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
#include "TrkReco/inc/TrkUtilities.hh"

//BaBar
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

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
  CalSeedFit::CalSeedFit(fhicl::ParameterSet const& pset) :
    _diagLevel       (pset.get<int>                ("diagLevel"                      )),
    _debugLevel      (pset.get<int>                ("debugLevel"                     )),
    _printfreq       (pset.get<int>                ("printFrequency"                 )),
    _useAsFilter     (pset.get<int>                ("useAsFitler"                    )),    
    _rescueHits      (pset.get<int>                ("RescueHits"                     )),    
    _addhits         (pset.get<bool>               ("addhits"                        )),
    _shLabel         (pset.get<string>             ("StrawHitCollectionLabel"        )),
    _shpLabel        (pset.get<string>             ("StrawHitPositionCollectionLabel")),
    _shfLabel        (pset.get<string>             ("StrawHitFlagCollectionLabel"    )),
    _helixSeedLabel    (pset.get<string>           ("HelixSeedModuleLabel"           )),
    _tpeaksLabel     (pset.get<string>             ("CalTimePeakModuleLabel"         )),
    _maxdtmiss       (pset.get<double>             ("MaxDtMiss"                      )),
    _maxadddoca      (pset.get<double>             ("MaxAddDoca"                     )),
    _maxaddchi       (pset.get<double>             ("MaxAddChi"                      )),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle"               ))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _seedfit         (pset.get<fhicl::ParameterSet>("SeedFitHack",fhicl::ParameterSet())),
    _sfresult(0)
  {
    //    fStopwatch = new TStopwatch();

    // tag the data product instance by the direction
    // and particle type found by this fitter

    produces<KalSeedCollection>     ();
    produces<AlgorithmIDCollection> ();

    //    produces<StrawHitFlagCollection>(_iname);

    _minNMCHits = 25;

    fHackData = new THackData("HackData","Hack Data");
    gROOT->GetRootFolder()->Add(fHackData);
    //-----------------------------------------------------------------------------
    // provide for interactive disanostics
    //-----------------------------------------------------------------------------
    _ref = new Ref("CalSeedFitRef","Ref to CalSeedFit", &_seedfit);

    TFolder* f_mu2e;

    f_mu2e = (TFolder*) gROOT->GetRootFolder()->FindObject("Mu2e");
    if (f_mu2e == NULL) f_mu2e = gROOT->GetRootFolder()->AddFolder("Mu2e","Mu2e Folder");

    if (f_mu2e) {
      _folder = f_mu2e->AddFolder("CalSeedFit","CalSeedFit Folder");
      _folder->Add(_ref);
    }

    fgTimeOffsets     = new SimParticleTimeOffset(pset.get<fhicl::ParameterSet>("TimeOffsets"));

    _helTraj = 0;
  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  CalSeedFit::~CalSeedFit() {
    delete _ref;
    if (_helTraj) delete _helTraj;
    //    delete fStopwatch;
  }

  //-----------------------------------------------------------------------------
  void CalSeedFit::beginJob(){

    if(_diagLevel > 0) bookHistograms();

    _eventid = 0;
  }

  //-----------------------------------------------------------------------------
  bool CalSeedFit::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    // calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();

    _seedfit.setTracker(_tracker);
    _seedfit.setTrackerCalib(_trackerCalib);
    
    return true;
  }

  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalSeedFit::bookHistograms() {
    if (_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      
      art::TFileDirectory sf_dir = tfs->mkdir("SeedFit");
      _hist.seedFit.nhits       = sf_dir.make<TH1F>("nhits" , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
      _hist.seedFit.seeddoca    [0]     = tfs->make<TH1F>("hseeddoca_0","doca seedfit active hits; doca [mm]", 1000, -20., 20);
      _hist.seedFit.seeddoca    [1]     = tfs->make<TH1F>("hseeddoca_1","doca seedfit non active hits; doca [mm]", 1000, -20., 20);
      _hist.seedFit.seeddoca    [2]     = tfs->make<TH1F>("hseeddoca_2","doca seedfit all hits; doca [mm]", 1000, -20., 20);

      _hist.seedFit.NpointsSeed [0]     = tfs->make<TH1F>("hnpseed_0","# of points de-activatedby seedfit; N points de-activated [#]", 50, 0., 50);
      _hist.seedFit.NpointsSeed [1]     = tfs->make<TH1F>("hnpseed_1","# of points not found in the pattern-reco and rescued ; N points rescued [#]", 50, 0., 50);
      _hist.seedFit.chi2[0]     = sf_dir.make<TH1F>("chi20"  , "chi2 distribution: all tacks", 100, 0., 10.);
      _hist.seedFit.chi2[1]     = sf_dir.make<TH1F>("chi21"  , "chi2 distribution:  tacks with nhits>15", 100, 0., 10.);
      _hist.seedFit.p   [0]     = sf_dir.make<TH1F>("p0"  , "p distribution: all tacks", 400, 0., 200.);
      _hist.seedFit.p   [1]     = sf_dir.make<TH1F>("p1"  , "p distribution:  tacks with nhits > 15", 400, 0., 200.);
      
      
      _hist.ntracks      [0]    = tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
      _hist.ntracks      [1]    = tfs->make<TH1F>("nseeds1"  , "number of track candidates with nhits > 15 ", 21, -0.5, 20.5);
    }
  }

  //-----------------------------------------------------------------------------
  // find the input data objects
  //-----------------------------------------------------------------------------
  bool CalSeedFit::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitPositionCollection with label=%s not found.\n",
             _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalSeedFit::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }

    //-----------------------------------------------------------------------------
    // find list of MC hits - for debugging only
    //-----------------------------------------------------------------------------
    art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_shLabel,"StrawHitMCPtr",mcptrHandle);
    if (mcptrHandle.isValid()) {
      _listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      _listOfMCStrawHits = NULL;
    }

    if (evt.getByLabel(_helixSeedLabel, _helixSeedsHandle)){
      _helixSeeds = _helixSeedsHandle.product();
    }else {
      _helixSeeds = 0;
    }
    
    art::Handle<CalTimePeakCollection> tpeaksHandle;
    if (evt.getByLabel(_tpeaksLabel, tpeaksHandle)){
      _tpeaks = tpeaksHandle.product();
    }else {
      _tpeaks = 0;
    }

    //-----------------------------------------------------------------------------
    // done
    //-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_helixSeeds != 0) && (_tpeaks != 0);
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  bool CalSeedFit::filter(art::Event& event ) {
    const char*               oname = "CalSeedFit::produce";
    char                      message[200];
    bool                      findhelix (false), findseed (false), findkal (false);
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

    //     t1 = fStopwatch->RealTime();
    //     fStopwatch->Continue();
    // event printout
    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    unique_ptr<KalSeedCollection>      tracks   (new KalSeedCollection     );
    unique_ptr<AlgorithmIDCollection>  algs     (new AlgorithmIDCollection);

    double pEntrance(.0), step_time(-9999.);
    double time_threshold(500.);

    const HelixSeed*   helixSeed(0);
    CalTimePeak* tp(0);

    // find the data
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n",oname);
      goto END;
    }
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

    _seedfit.setStepPointMCVectorCollection(_listOfMCStrawHits);

    //-----------------------------------------------------------------------------
    // loop over found time peaks - for us, - "eligible" calorimeter clusters
    //-----------------------------------------------------------------------------
    ntrkcandidates = _helixSeeds->size();
    
    for (int ipeak=0; ipeak<ntrkcandidates; ipeak++) {
      helixSeed = &_helixSeeds->at(ipeak);
      tp      = (CalTimePeak*)&_tpeaks  ->at(ipeak);

      TrkDefHack             seeddef(_shcol, helixSeed->timeCluster()->hits(), _tpart, _fdir);

      // track fitting objects for this track candidate
      init(_sfresult, &seeddef);

      findhelix = true;

      //-----------------------------------------------------------------------------
      // P.Murat: here hits are ordered by index - WHY?
      // the Kalman fitter needs them ordered in Z(straw)
      //-----------------------------------------------------------------------------
      _hitIndices = helixSeed->timeCluster()->hits();
      
      std::sort(_hitIndices.begin(), _hitIndices.end(), [ ]( const StrawHitIndex& lhs,
							     const StrawHitIndex& rhs )
		{
		  return lhs < rhs;
		} );
      
      _nindex = _hitIndices.size();
      
      seeddef.setIndices (_hitIndices);

      //-----------------------------------------------------------------------------
      // seed fit - fit through the wires of found hits, not using the drift times
      //-----------------------------------------------------------------------------
      _seedfit.makeTrack(*_sfresult, tp);

      //--------------------------------------------------------------------------------
      // 2014-11-24 gianipez added the following diagnnostic
      //--------------------------------------------------------------------------------
      if (_debugLevel > 0) {
	sprintf(message,
		"CalSeedFit::produce seedfit::makeTrack : fit_success = %i\n",
		_sfresult->_fit.success());
	_seedfit.printHits(*_sfresult,message);
      }
      
      

      if (_sfresult->_fit.success()){
	_ntracks[0] += 1;


	if (_rescueHits == 1){
// //-----------------------------------------------------------------------------
// // use helix parameters by the seed fit to initialize the full Kalman fit
// //-----------------------------------------------------------------------------
//           double           locflt;
//           const HelixTraj* shelix;
//           shelix = (const HelixTraj*) _sfresult->_krep->localTrajectory(_sfresult->_krep->flt0(),locflt);
//           seeddef.setHelix(*shelix);
//----------------------------------------------------------------------
//2015-02-07 G. Pezzu added new selection using seedFit results
//----------------------------------------------------------------------
	  _goodhits.clear();

          const mu2e::TrkStrawHit* hit;
          int                      hit_index;
	  TrkHitVector const&  hot_l = _sfresult->_krep->hitVector();
          const StrawHit*          sh;
          const Straw*             straw;
          const CLHEP::Hep3Vector  *wpos, *wdir;

                                        //  Trajectory info
          Hep3Vector               tdir;
          HepPoint                 tpos;
          double                   doca, rdrift, fltlen;
          bool                     found(false), active;
          int                      banner_11_printed(0);

          _sfresult->_krep->traj().getInfo(0.0,tpos,tdir);

          _nrescued = 0;

          for (int i=0; i< _nindex; ++i) {
            hit_index = _hitIndices[i];
            sh        = &_shcol->at(hit_index);
            straw     = &_tracker->getStraw(sh->strawIndex());
            wpos      = &straw->getMidPoint();
            wdir      = &straw->getDirection();
            rdrift    = -9990;
            found     = false;
            active    = false;
	    
            HepPoint      wpt  (wpos->x(),wpos->y(),wpos->z());
            TrkLineTraj   wtraj(wpt,*wdir,-20,20);
//-----------------------------------------------------------------------------
// estimate flightlength along track. This assumes a constant BField!!!
// in principle, this should work well enough, however, may want to check
//-----------------------------------------------------------------------------
            fltlen = (wpos->z()-tpos.z())/tdir.z();

            TrkPoca   wpoca(_sfresult->_krep->traj(),fltlen,wtraj,0.0);

            doca = wpoca.doca();

            for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
              hit    = static_cast<const mu2e::TrkStrawHit*> (*it);
              rdrift = hit->driftRadius();
              int shIndex = int(hit->index());
              if (hit_index == shIndex) {
                found = true;
                break;
              }
            }

            if ( std::fabs(doca) < _maxadddoca) {
              active = true;
              _goodhits.push_back(hit_index);
            }

            if (_debugLevel > 0) {
              if (banner_11_printed == 0) {
                banner_11_printed = 1;
                printf("[CalPatRec::produce] -------------------------------------\n");
                printf("[CalPatRec::produce]  ih  A   Sind      Rdrift        doca\n");
                printf("[CalPatRec::produce] -------------------------------------\n");
              }

              printf("[CalPatRec::produce]  %2i  %1i  %5i  %10.3f  %10.3f \n",
                     i, active? 1:0, straw->index().asInt(), rdrift, doca );
            }

            if (_diagLevel > 0) {
              _hist.seedFit.seeddoca[2]->Fill(doca);

              if (found)  _hist.seedFit.seeddoca[0]->Fill(doca);
              else        _hist.seedFit.seeddoca[1]->Fill(doca);

              if (active) _hist.seedFit.seeddoca[0]->Fill(doca);
              else        _hist.seedFit.seeddoca[1]->Fill(doca);
            }

            if (!found && active) ++_nrescued;
          }

       
//-----------------------------------------------------------------------------
// at this point the full kalman fit starts
//-----------------------------------------------------------------------------
          seeddef.setIndices(_goodhits);
          if (_debugLevel > 0) printf("CalPatRec::produce] calling _kfit.makeTrack\n");

          if (_nrescued > 0) _seedfit.makeTrack(*_sfresult, tp);//repeat the fit only if new hits have been found
	  
	  if ( !(_sfresult->_fit.success()) ){//if adding the new hits the fit does not converge, return to the previous case
	    seeddef.setIndices (_hitIndices);
	    _seedfit.makeTrack(*_sfresult, tp);
	  }
	  
	  if (_diagLevel > 0) fillSeedFitHistograms(*_sfresult);
	}



	if (_diagLevel > 0){
	  TrkHitVector const& hot_l = _sfresult->_krep->hitVector();
	  int                 nactivated(0);

	  const mu2e::TrkStrawHit* hit;
	  for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
	    hit = static_cast<const mu2e::TrkStrawHit*> (*it);
	    if (hit->isActive()) ++nactivated;
	  }
	    
	  double  chi2           = _sfresult->_krep->chisq()/nactivated;

	  double  h1_fltlen      = _sfresult->_krep->firstHit()->kalHit()->hit()->fltLen();
	  double  hn_fltlen      = _sfresult->_krep->lastHit ()->kalHit()->hit()->fltLen();
	  double  entlen         = std::min(h1_fltlen, hn_fltlen);

	  CLHEP::Hep3Vector fitmom = _sfresult->_krep->momentum(entlen);

	  _hist.seedFit.nhits  ->Fill(nactivated);
	  _hist.seedFit.chi2[0]->Fill(chi2);     
	  _hist.seedFit.p   [0]->Fill(fitmom.mag());     

	  if (nactivated > 15){
	    _ntracks[1] += 1;
	    _hist.seedFit.chi2[1]->Fill(chi2);     
	    _hist.seedFit.p   [1]->Fill(fitmom.mag());     
	  }
	}

	krep      = _sfresult->stealTrack();
	// convert the status into a FitFlag
	TrkFitFlag seedok(TrkFitFlag::seedOK);
	KalSeed    kseed (_tpart, _fdir, krep->t0(), krep->flt0(), seedok);
	// fill ptr to the helix seed
	kseed._helix = art::Ptr<HelixSeed>(_helixSeedsHandle, ipeak);
	// calo cluser ptr from Helix Seed
	kseed._caloCluster = helixSeed->caloCluster(); 
	// extract the hits from the rep and put the hitseeds into the KalSeed
	TrkUtilities::fillHitSeeds(krep,kseed._hits);
	kseed._status.merge(TrkFitFlag::hitsOK);

	// extract the helix trajectory from the fit (there is just 1)
	double locflt;
	const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(krep->flt0(),locflt));
	// use this to create segment.  This will be the only segment in this track
	if(htraj != 0){
	  KalSegment kseg;
	  // sample the momentum at this point
	  BbrVectorErr momerr = krep->momentumErr(krep->flt0());
	  TrkUtilities::fillSegment(*htraj, momerr, kseg);
	  kseed._segments.push_back(kseg);
	  // push this seed into the collection
	  tracks->push_back(kseed);
	  
	  if(_debugLevel > 1){
	    printf("Seed fit segment parameters \n");
	    for(size_t ipar=0;ipar<5;++ipar) cout << kseg.helix()._pars[ipar] << " ";
	    printf(" covariance \n");
	    for(size_t ipar=0;ipar<15;++ipar)
	      cout << kseg.covar()._cov[ipar] << " ";
	    cout << endl;
	  }
	} else {
	  throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
	}
	
	tp->SetCprIndex(tracks->size());
	    
	int best = AlgorithmID::CalPatRecBit;
	int mask = 1 << AlgorithmID::CalPatRecBit;
	    
	algs->push_back(AlgorithmID(best,mask));
	    
      }else {
	_sfresult->deleteTrack();
      }
	  
//-----------------------------------------------------------------------------
// cleanup the seed fit - why it is not being done ?
//-----------------------------------------------------------------------------
//       _sfresult->deleteTrack();

      if (_debugLevel > 0) {
	if (_nhits_from_gen >= _minNMCHits) {
	  if (tp->_tmin > 400.){
	    if (!findhelix) {
	      printf("[CalSeedFit::produce] LOOK AT: more than 25 MC hits and findHelix not converged! event = %i\n", _iev);
	    }
	    if (findhelix && !findseed){
	      printf("[CalSeedFit::produce] LOOK AT: findhelix converged and findseed not! event = %i\n", _iev);
	    }
	    if (findseed && !findkal){
	      printf("[CalSeedFit::produce] LOOK AT: findseed converged and findkal not! event = %i\n", _iev);
	    }
	  }
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
  END:;
    int     ntracks = tracks->size();

    event.put(std::move(tracks));
    event.put(std::move(algs  ));

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
  void CalSeedFit::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }



//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalSeedFit::findMissingHits(KalFitResult& kalfit,std::vector<StrawHitIndex>& misshits) {
                                        //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    int        radius_ok;
    double     dt;

    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    if (_debugLevel > 0) {
      printf("[CalSeedFit::findMissingHits]      shId    sec     panel       doca   \n");
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
            printf("[CalSeedFit::findMissingHits] %8i  %6i  %8i  %10.3f \n",
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
// seed fit diagnostic histograms, DOCA is the distance to the wire
//-----------------------------------------------------------------------------
  void CalSeedFit::fillSeedFitHistograms(KalFitResult& SFResult) {

    int                           ndeactivated(0);
    TrkHitVector const&           hot_l = SFResult._krep->hitVector();
    const mu2e::TrkStrawHit*      hit;
 
    for (auto it=hot_l.begin(); it<hot_l.end(); it++) {
      hit = static_cast<const mu2e::TrkStrawHit*> (*it);
      if (!hit->isActive()) ++ndeactivated;
    }

    _hist.seedFit.NpointsSeed[0]->Fill(ndeactivated);
    _hist.seedFit.NpointsSeed[1]->Fill(_nrescued);
  }

//-----------------------------------------------------------------------------
  void CalSeedFit::init(KalFitResult*& KRes, TrkDefHack* TDef) {

    if (KRes != 0) {
      KRes->deleteTrack();
      delete KRes;
    }
    KRes = new KalFitResult();
    KRes->_tdef = TDef;
  }
}

using mu2e::CalSeedFit;
DEFINE_ART_MODULE(CalSeedFit);
