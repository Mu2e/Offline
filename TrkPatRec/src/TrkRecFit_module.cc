//
// TTracker Kalman Fit launcher
//
// $Id: TrkRecFit_module.cc,v 1.1 2014/08/22 16:49:03 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:49:03 $
//
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// data
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPoca.hh"
// Mu2e
#include "TrkPatRec/inc/TrkRecBase.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "FastPatternReco/inc/TrackSeedUtils.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
// C++
#include <iostream>

using namespace std;

namespace mu2e 
{
  class TrkRecFit : public TrkRecBase, public art::EDProducer
  {
    public:
      explicit TrkRecFit(fhicl::ParameterSet const&);
      virtual ~TrkRecFit();
      virtual void beginJob();
//      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event );
//    virtual  void endJob();
    private:
      // configuration parameters
      art::Handle<mu2e::StrawHitCollection> _strawhitsH;
      art::Handle<TrackerHitTimeClusterCollection> _tclusthitH;
      art::Handle<TrackSeedCollection> _trksSeedH;
      // event object labels
      string _tpkfLabel;
      string _ptrnRecLabel; // Label to the Pattern Recognition module.
      // cache of event objects
      const TrackerHitTimeClusterCollection* _tccol;
      const TrackSeedCollection * _tscol;
      bool findData(const art::Event& e);
  };

  TrkRecFit::TrkRecFit(fhicl::ParameterSet const& pset) :
     TrkRecBase(pset),
    _tpkfLabel(pset.get<string>("TrackerHitTimeClusterCollection","TimePeakFinder")),
    _ptrnRecLabel(pset.get<string>("PatternRecoModuleLabel","RobustHelixFinder"))
  {
    produces<KalRepCollection>(_iname);
    produces<KalRepPtrCollection>(_iname);
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>(_iname);
  }

  TrkRecFit::~TrkRecFit(){}

  void TrkRecFit::beginJob(){
    bgnJob();
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
    if((_iev%_printfreq)==0)cout<<"TrkRecFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    unique_ptr<StrawHitFlagCollection> flags(_flags );
    // find mc truth if we're making diagnostics
    if(_diag > 0 && !_kfitmc.findMCData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkRecFit: MC information missing "<< endl;
    }
    if(_diag > 1){
      fillStrawDiag();
      if(_nchit>14)_ccutflow->Fill(1.0);
      if(_nchit>14&&_ctime>_tmin)_ccutflow->Fill(2.0);
    }

    // find the time peaks in the time spectrum of selected hits.  Otherwise, take all
    // selected hits as a peak
    _tpeaks.clear();
    _icepeak = -1;

    //if(_findtpeak){
      loadTimePeaks(_tpeaks,_tccol);
      if(_diag>1 && _kfitmc.mcData()._mcsteps != 0){
        for(size_t ip=0;ip<_tpeaks.size();++ip){
          TrkTimePeak const& tp = _tpeaks[ip];
          fillPeakDiag(ip,tp);
        }
      }
    //}

    if(_diag>0){
// fill primary particle MC truth information
      _kfitmc.mcTrkInfo(_kfitmc.mcData()._simparts->begin()->second);
    }
    // fill diagnostics if requested
    if(_diag > 2 && _nchit>0 && _ctime > _tmin)fillTimeDiag();
    // dummy objects
    static TrkDef dummydef;
    static HelixDef dummyhdef;
    static HelixFitResult dummyhfit(dummyhdef);
    static KalFitResult dummykfit(dummydef);
    // loop over the accepted time peaks
    if(_tpeaks.size()>0)_cutflow->Fill(1.0);
    if(_diag>1 && _icepeak >=0)_ccutflow->Fill(3.0);
    bool findhelix(false), findseed(false), findkal(false);

    if (_debug>0 || (_iev%_printfreq)==0) {
            std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"event "<<_iev<<" tot N hit "<<_shcol->size()<<" N tracks seed found "<<_tscol->size()
                            <<" N time peaks "<<_tccol->size()<<std::endl;
            std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
    }

    for (size_t iTrackSeed=0; iTrackSeed<_tscol->size(); ++iTrackSeed) {
      TrackSeed const& iTrkSeed(_tscol->at(iTrackSeed));
      unsigned ipeak = iTrkSeed._relatedTimeCluster.key();
      findhelix = true;

      //Fake conversion !!! Needed because I had to redefine hitIndex as different type (HitIndex) to avoid compilation problem!! FIXME
      //all time peaks
//      std::vector<hitIndex> tpeakhits;
//      const TrackerHitTimeCluster&  tclust=*(iTrkSeed._relatedTimeCluster);
//      for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin();
//           iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
//        size_t iglbHit = iTCHit->key();
//        tpeakhits.push_back(mu2e::hitIndex(iglbHit));
//      }
      //all track selected
      std::vector<hitIndex> goodhits;

      const std::vector<HitIndex> &trkseedhits = iTrkSeed._fullTrkSeed._selectedTrackerHitsIdx;
      for (std::vector<HitIndex>::const_iterator loopPoints_it = trkseedhits.begin();
        loopPoints_it != trkseedhits.end(); ++loopPoints_it) {
        goodhits.push_back( mu2e::hitIndex(loopPoints_it->_index,loopPoints_it->_ambig) );
        /*unsigned int i=0;
        for(i=0;i<strawhits.size();i++){
          if(strawhits[i]._index==(*loopPoints_it)._index) {
            goodhits.push_back( mu2e::hitIndex((*loopPoints_it)._index,strawhits[i]._ambig) );
            break;
          }
        }
        if(i==strawhits.size()) {
          goodhits.push_back( mu2e::hitIndex((*loopPoints_it)._index));
        }*/
      }

      HelixTraj recoseed(TrkParams(HelixTraj::NHLXPRM));
      HelixVal2HelixTraj(iTrkSeed._fullTrkSeed,recoseed);


      TrkDef seeddef(_shcol,goodhits,recoseed,_tpart,_fdir);

//      TrkDef seeddef(_shcol,_tpart,_fdir);
//      seeddef.setHelix(recoseed);
//      const std::vector<HitIndex> &trkseedhits = iTrkSeed._fullTrkSeed._selectedTrackerHitsIdx;
//      for (std::vector<HitIndex>::const_iterator loopPoints_it = trkseedhits.begin();
//        loopPoints_it != trkseedhits.end(); ++loopPoints_it) {
//        seeddef.appendHit(loopPoints_it->_index,loopPoints_it->_ambig);
//      }

      seeddef.setEventId(_iev);
      seeddef.setTrackId(iTrackSeed);

      TrkT0 t0(iTrkSeed._t0,iTrkSeed._errt0);
      //TrkT0 t0(tclust._meanTime,tclust._sigma);
      //if(_t0use==1) t0=TrkT0(_histoOut.recoinfo.t0,0.0);
      ////if(_t0use==2) t0=TrkT0(_histoOut.recoinfo.t0,iTrkSeed._errt0);
      //if(_t0use==2) t0=TrkT0(iTrkSeed._t0,iTrkSeed._errt0);

      seeddef.setT0(t0);
      TrkDef kaldef(seeddef);
      HelixFitResult helixfit(seeddef);
      KalFitResult seedfit(seeddef);
      KalFitResult kalfit(kaldef);
      // initialize filters.  These are used only for diagnostics
      _hfilt.clear();
      _sfilt.clear();
      // now, fit the seed helix from the filtered hits
      //seedfit._tdef.helix().printAll(std::cout);
      _seedfit.makeTrack(seedfit);
      if(seedfit._fit.success()){
        findseed = true;
        cout<<"Seed found"<<endl;
        // find the helix parameters from the helix fit, and initialize the full Kalman fit with this
        double locflt;
        const HelixTraj* shelix = dynamic_cast<const HelixTraj*>(seedfit._krep->localTrajectory(seedfit._krep->flt0(),locflt));
        kaldef.setHelix(*shelix);
        // filter the outliers
        filterOutliers(kaldef,seedfit._krep->traj(),_maxseeddoca,_diag,&_sfilt,&_kfitmc);
        _kfit.makeTrack(kalfit);
        // if successfull, try to add missing hits
        if(kalfit._fit.success()){
          findkal = true;
          cout<<"Fit found"<<endl;
          if(_addhits){
            // first, add back the hits on this track
            _kfit.unweedHits(kalfit,_maxaddchi);
            vector<hitIndex> misshits;
            findMissingHits(kalfit,misshits);
            if(misshits.size() > 0){
              _kfit.addHits(kalfit,_shcol,misshits,_maxaddchi);
            }
          }
        }
      }
      // fill fit diagnostics if requested
      if(_diag > 0)
	fillFitDiag(ipeak,helixfit,seedfit,kalfit);
      if(_diag > 1 && (int)ipeak == _icepeak){
	if(helixfit._fit.success()){
	  _ccutflow->Fill(4.0);
	  _ccutflow->Fill(5.0);
	  _ccutflow->Fill(6.0);
	  _ccutflow->Fill(7.0);
	} else {
	  if(helixfit._fit.failure()>1)_ccutflow->Fill(4.0);
	  if(helixfit._fit.failure()>2)_ccutflow->Fill(5.0);
	  if(helixfit._fit.failure()>3)_ccutflow->Fill(6.0);
	}
	if(seedfit._fit.success())_ccutflow->Fill(8.0);
	if(kalfit._fit.success())_ccutflow->Fill(9.0);
      }
      if(kalfit._fit.success()){
	// flag the hits used in this track.  This should use the track id, FIXME!!! (in the BaBar code)
	if(ipeak<16){
	  for(size_t ihit=0;ihit<kalfit._hits.size();++ihit){
	    const TrkStrawHit* tsh = kalfit._hits[ihit];
	    if(tsh->isActive())_flags->at(tsh->index()).merge(StrawHitFlag::trackBit(ipeak));
	  }
	}
	// save successful kalman fits in the event
	tracks->push_back( kalfit.stealTrack() );
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
      } else
	kalfit.deleteTrack();
      // cleanup the seed fit
      seedfit.deleteTrack();
    }
    if(findhelix)_cutflow->Fill(2.0);
    if(findseed)_cutflow->Fill(3.0);
    if(findkal)_cutflow->Fill(4.0);
    // add a dummy entry in case there are no peaks
    if(_diag > 0 && _tpeaks.size() == 0)
      fillFitDiag(-1,dummyhfit,dummykfit,dummykfit);
    // put the tracks into the event
    art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*tracks, tracksID, event);
    event.put(move(tracks),_iname);
    event.put(move(trackPtrs),_iname);
    event.put(move(flags),_iname);
  }

  // find the input data objects 
  bool TrkRecFit::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _shpcol = 0;
    _stcol = 0;
    _tccol = 0;
    _tscol = 0;

    if(evt.getByLabel(_shLabel,_strawhitsH))
      _shcol = _strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StereoHitCollection> stH;
    if(evt.getByLabel(_stLabel,stH))
      _stcol = stH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
    if (evt.getByLabel(_tpkfLabel,_tclusthitH))
      _tccol = _tclusthitH.product();
    if (evt.getByLabel(_ptrnRecLabel,_trksSeedH))
      _tscol = _trksSeedH.product();
// don't require stereo hits: they are only used for diagnostics
    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && _tccol != 0 && _tscol != 0;
  }

}
using mu2e::TrkRecFit;
DEFINE_ART_MODULE(TrkRecFit);
