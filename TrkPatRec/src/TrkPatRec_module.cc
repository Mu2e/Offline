// $Id: TrkPatRec_module.cc,v 1.83 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/22 16:10:41 $
//

// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
// BaBar
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
// Mu2e
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
// C++
#include <iostream>

using namespace std;

// Mu2e
#include "TrkPatRec/inc/TrkRecBase.hh"

namespace mu2e 
{
  class TrkPatRec : public TrkRecBase, public art::EDProducer
  {
    public:
      explicit TrkPatRec(fhicl::ParameterSet const&);
      virtual ~TrkPatRec();
      virtual void beginJob();
//      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event );
//    virtual  void endJob();
  };

  TrkPatRec::TrkPatRec(fhicl::ParameterSet const& pset) :
     TrkRecBase(pset)
  {
    produces<KalRepCollection>(_iname);
    produces<KalRepPtrCollection>(_iname);
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>(_iname);
  }

  TrkPatRec::~TrkPatRec() {}

  void TrkPatRec::beginJob(){
    bgnJob();
  }

  void TrkPatRec::produce(art::Event& event ) {
    _eventid = event.event();
    _cutflow->Fill(0.0);
    if(_diag>1)_ccutflow->Fill(0.0);
    // create output
    unique_ptr<KalRepCollection>    tracks(new KalRepCollection );
    unique_ptr<KalRepPtrCollection> trackPtrs(new KalRepPtrCollection );
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event,_iname));
    // event printout
    _iev=event.id().event();
    if((_iev%_printfreq)==0)cout<<"TrkPatRec: event="<<_iev<<endl;
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
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: MC information missing "<< endl;
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
//    if(_findtpeak){
//      findTimePeaks();
//    } else {
//      createTimePeak();
//    }
    if (_debug>0) { cout<<"TrkPatRec::produce - starting TPFinder"<<endl; }
    if(_findtpeak){
     if (_debug>0) { cout<<"TrkPatRec::produce - finding TP"<<endl; }
     findTimePeaks( _shcol, _shpcol, _flags,
                     _tsel, _hsel, _tbkg, _hbkg,
                     _tpeaks, _maxnpeak,
                     _nbins, _tmin, _tmax,
                     _1dthresh, _ymin, _maxdt, _minnhits,
                     _cleanpeaks, _pmva, _peakMVA, _PMVAType, 
                     _minpeakmva, _maxpeakdt, _maxpeakdphi );
      // if requested, fill diagnostics
      if(_diag>1 && _kfitmc.mcData()._mcsteps != 0){
        for(size_t ip=0;ip<_tpeaks.size();++ip){
          TrkTimePeak const& tp = _tpeaks[ip];
          fillPeakDiag(ip,tp);
        }
      }

    } else {
      if (_debug>0) { cout<<"TrkPatRec::produce - creating TP"<<endl; }
      createTimePeak( _shcol, _flags,
                      _tsel, _hsel, _tbkg, _hbkg,
                      _tpeaks, _minnhits );

    }
    if (_debug>0) { cout<<"TrkPatRec::produce - TPFinder done"<<endl; }

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
    for(unsigned ipeak=0;ipeak<_tpeaks.size();++ipeak){
      // create track definitions for the helix fit from this initial information 
      HelixDef helixdef(_shcol,_shpcol,_tpeaks[ipeak]._trkptrs,_tpart,_fdir);
      // set some identifiers
      helixdef.setEventId(_eventid);
      helixdef.setTrackId(ipeak);
      // copy this for the other fits
      TrkDef seeddef(helixdef);
      TrkDef kaldef(helixdef);
      // track fitting objects for this peak
      HelixFitResult helixfit(helixdef);
      KalFitResult seedfit(seeddef);
      KalFitResult kalfit(kaldef);
      // initialize filters.  These are used only for diagnostics
      _hfilt.clear();
      _sfilt.clear();
      // robust helix fit
      if(_hfit.findHelix(helixfit,_icepeak==(int)ipeak)){
	findhelix = true;
	// convert the result to standard helix parameters, and initialize the seed definition helix
	HepVector hpar;
	HepVector hparerr;
	_hfit.helixParams(helixfit,hpar,hparerr);
	HepSymMatrix hcov = vT_times_v(hparerr);
	seeddef.setHelix(HelixTraj(hpar,hcov));
	// Filter outliers using this helix
	filterOutliers(seeddef,seeddef.helix(),_maxhelixdoca,_diag,&_hfilt);
	// now, fit the seed helix from the filtered hits
	_seedfit.makeTrack(seedfit);
	if(seedfit._fit.success()){
	  findseed = true;
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

}
using mu2e::TrkPatRec;
DEFINE_ART_MODULE(TrkPatRec);
