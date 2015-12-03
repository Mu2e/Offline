// $Id: TrkPatRec_module.cc,v 1.86 2014/09/18 09:34:08 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/18 09:34:08 $
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
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalDiag.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "TrkPatRec/inc/HelixFit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
// Mu2e
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include "TMarker.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
using namespace std; 
using namespace boost::accumulators;

namespace mu2e 
{
  class TrkPatRec : public art::EDProducer
  {
    public:
      enum fitType {helixFit=0,seedFit,kalFit};
      explicit TrkPatRec(fhicl::ParameterSet const&);
      virtual ~TrkPatRec();
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
      void endJob();
    private:
      unsigned _iev;
  // diagnostics
      KalDiag* _kdiag;
      // configuration parameters
      int _diag,_debug;
      int _printfreq;
      bool _addhits; 
      bool _useobsolete;
      // event object labels
      string _shLabel;
      string _shpLabel;
      string _stLabel;
      string _shfLabel;
      StrawHitFlag _tsel, _hsel, _addsel;
      StrawHitFlag _tbkg, _hbkg, _addbkg;
      double _maxdt, _maxdtmiss;
      // time spectrum parameters
      bool _findtpeak;
      unsigned _maxnpeak;
      unsigned _minnhits;
      bool _cleanpeaks;
      double _minpeakmva, _maxpeakdt, _maxpeakdphi;
      std::string _PMVAType; // type of MVA
      std::string _PMVAWeights; // file of MVA weights
      double _maxphirange;
      double _tmin;
      double _tmax;
      double _tbin;
      unsigned _nbins;
      double _ymin;
      double _1dthresh,_tssigma;
      // outlier cuts
      double _maxseeddoca,_maxhelixdoca,_maxadddoca, _maxaddchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      StrawHitFlagCollection* _flags;
      const StrawHitPositionCollection* _shpcol;
      const StereoHitCollection* _stcol;
      // Kalman fitters.  Seed fit has a special configuration
      KalFit _seedfit, _kfit;
      // robust helix fitter
      HelixFit _hfit;
      // cache of time peaks
      vector<TrkTimePeak> _tpeaks;
      string _iname; // data instance name
      //
      PayloadSaver _payloadSaver;
      // helper functions
      bool findData(const art::Event& e);
      void findTimePeaks();
      void findTimePeaks(TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak);
      void createTimePeak();
      void cleanTimePeak(TrkTimePeak& tp);
      void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec);
      void findMissingHits(KalFitResult& kalfit, vector<hitIndex>& indices);
      void createDiagnostics();
      void fillStrawDiag();
      void fillTimeDiag();
      void fillPeakDiag(size_t ip, TrkTimePeak const& tp);
      void fillFitDiag(int ipeak, HelixFitResult const& helixfit,
	  KalFitResult const& seedfit,KalFitResult const& kalfit);
      void fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const;
      void initializeReaders();

      TMVA::Reader *_peakMVA; // MVA for peak cleaning
      TimePeakMVA _pmva; // input variables to TMVA for peak cleaning

      // strawhit tuple variables
      TTree *_shdiag;
      Int_t _eventid;
      threevec _shp;
      Float_t _edep;
      Float_t _time, _deltat, _rho;
      Int_t _nmcsteps;
      Int_t _mcnunique,_mcnmax;
      Int_t _mcpdg,_mcgen,_mcproc;
      Int_t _mcppdg,_mcpproc;
      Int_t _mcgid, _mcgpdg;
      Float_t _mcge, _mcgt;
      threevec _mcshp, _mcop, _mcpop, _mcgpos;
      Float_t _mcoe, _mcpoe, _mcom, _mcpom;
      Float_t _mcshlen,_mcshd;
      Float_t _mcedep;
      Float_t _pdist,_pperp,_pmom;
      Float_t _mctime, _mcptime;
      Int_t _esel,_rsel, _timesel,  _delta, _stereo, _isolated;
      Int_t _plane, _panel, _layer, _straw;
      Float_t _shpres, _shrres, _shchisq, _shdt, _shdist;
      Bool_t _xtalk;
      // time peak diag variables
      TTree* _tpdiag;
      Int_t _tpeventid, _peakid, _pmax, _nphits, _ncphits, _nchits;
      Float_t _ptime, _pdtimemax, _ctime, _cdtimemax;;
      Float_t _pphi, _cphi, _cphirange, _pdphimax, _cdphimax;
      vector<TimePeakHitInfo> _tphinfo;

      // fit tuple variables
      Int_t _nadd,_ipeak;
      Float_t _hcx, _hcy, _hr, _hdfdz, _hfz0;
      Float_t _mccx, _mccy, _mcr, _mcdfdz, _mcfz0;
      Int_t _helixfail,_seedfail,_kalfail;
      helixpar _hpar,_spar;
      helixpar _hparerr,_sparerr;
      Int_t _snhits, _snactive, _sniter, _sndof, _snweediter;
      Float_t _schisq, _st0;
      Int_t _nchit;
      Int_t _npeak, _nmc;
      Float_t _peakmax, _tpeak;
      // hit filtering tuple variables
      vector<TrkHitFilter> _sfilt, _hfilt;
      // flow diagnostic
      TH1F* _cutflow, *_ccutflow;
      int _icepeak;
  };

  TrkPatRec::TrkPatRec(fhicl::ParameterSet const& pset) :
    _kdiag(new KalDiag(pset.get<fhicl::ParameterSet>("KalDiag"))),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
    _useobsolete(pset.get<bool>("UseObsolete",false)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _stLabel(pset.get<string>("StereoHitCollectionLabel","MakeStereoHits")),
    _shfLabel(pset.get<string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _tsel(pset.get<vector<string> >("TimeSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"} )),
    _hsel(pset.get<vector<string> >("HelixFitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"} )),
    _addsel(pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _tbkg(pset.get<vector<string> >("TimeBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _hbkg(pset.get<vector<string> >("HelixFitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _addbkg(pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _maxdt(pset.get<double>("DtMax",30.0)),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _findtpeak(pset.get<bool>("FindTimePeaks",true)),
    _maxnpeak(pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits(pset.get<unsigned>("MinNHits",10)),
    _cleanpeaks(pset.get<bool>("CleanTimePeaks",true)),
    _minpeakmva(pset.get<double>("MinTimePeakMVA",0.5)),
    _maxpeakdt(pset.get<double>("MaxTimePeakDeltat",25.0)),
    _maxpeakdphi(pset.get<double>("MaxTimePeakDeltaPhi",1.0)),
    _PMVAType(pset.get<std::string>("TimePeakMVAType","MLP method")),
    _tmin(pset.get<double>("tmin",500.0)),
    _tmax(pset.get<double>("tmax",1700.0)),
    _tbin(pset.get<double>("tbin",20.0)),
    _ymin(pset.get<double>("ymin",8.0)),
    _1dthresh(pset.get<double>("OneDPeakThreshold",5.0)),
    _maxseeddoca(pset.get<double>("MaxSeedDoca",10.0)),
    _maxhelixdoca(pset.get<double>("MaxHelixDoca",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet())),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet()),_kdiag),
    _hfit(pset.get<fhicl::ParameterSet>("HelixFit",fhicl::ParameterSet())),
    _payloadSaver(pset)
    {
    if(!_useobsolete){
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: This module is obsolete.  Use one of the TrkPatRec sequences defined in TrkPatRec/fcl/prolog.fcl"<< endl;
    }

    // tag the data product instance by the direction and particle type found by this fitter
    _iname = _fdir.name() + _tpart.name();
    produces<KalRepCollection>(_iname);
    produces<KalRepPtrCollection>(_iname);
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>(_iname);

    //    produces<KalFitResultCollection>(_iname);

    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    // location-independent files
    ConfigFileLookupPolicy configFile;
    std::string weights = pset.get<std::string>("PeakMVAWeights","TrkPatRec/test/TimePeak.weights.xml");
    _PMVAWeights = configFile(weights);
  }

  TrkPatRec::~TrkPatRec(){}

  void TrkPatRec::beginJob(){
    initializeReaders();
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

  void TrkPatRec::beginRun(art::Run& ){}

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
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: data missing or incomplete"<< endl;
    }
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    unique_ptr<StrawHitFlagCollection> flags(_flags );
    //    unique_ptr<KalFitResultCollection> kfresults(new KalFitResultCollection);
    // find mc truth if we're making diagnostics
    if(_diag > 0 ){
      bool goodmc = _kdiag->findMCData(event);
      if(!goodmc)
	throw cet::exception("RECO")<<"mu2e::TrkPatRec: MC data missing or incomplete"<< endl;
// fill primary particle MC truth information
      _kdiag->kalDiag(0,false);
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
    if(_findtpeak){
      findTimePeaks();
    } else {
      createTimePeak();
    }
    // fill diagnostics if requested
    if(_diag > 2 && _nchit>0 && _ctime > _tmin)fillTimeDiag();
    // dummy objects
    static TrkDef dummydef;
    static HelixDef dummyhdef;
    static HelixFitResult dummyhfit(dummyhdef);
    static KalFitResult dummykfit(&dummydef);
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
      KalFitResult seedfit(&seeddef);
      KalFitResult kalfit(&kaldef);
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
	filterOutliers(seeddef,seeddef.helix(),_maxhelixdoca,_hfilt);
	// now, fit the seed helix from the filtered hits
	_seedfit.makeTrack(seedfit);
	if(seedfit._fit.success()){
	  findseed = true;
	  // find the helix parameters from the helix fit, and initialize the full Kalman fit with this
	  double locflt;
	  const HelixTraj* shelix = dynamic_cast<const HelixTraj*>(seedfit._krep->localTrajectory(seedfit._krep->flt0(),locflt));
	  kaldef.setHelix(*shelix);
	  // filter the outliers
	  filterOutliers(kaldef,seedfit._krep->traj(),_maxseeddoca,_sfilt);
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
		// 2015-04-12 P.Murat: assume this is a call corresponding to final=1
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
	//	kfresults->push_back(kalfit);
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
    //    event.put(move(kfresults),_iname);
  }

  void TrkPatRec::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

  // find the input data objects 
  bool TrkPatRec::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _shpcol = 0;
    _stcol = 0;
    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(evt.getByLabel(_shLabel,strawhitsH))
      _shcol = strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StereoHitCollection> stH;
    if(evt.getByLabel(_stLabel,stH))
      _stcol = stH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
// don't require stereo hits: they are only used for diagnostics
    return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
  }

  void TrkPatRec::findTimePeaks() {
    TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
    // loop over straws hits and fill time spectrum plot for tight hits
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
	double time = _shcol->at(istr).time();
	timespec.Fill(time);
      }
    }
    std::vector<Float_t> xpeaks,ypeaks;
    findTimePeaks(timespec,xpeaks,ypeaks);
    unsigned np = xpeaks.size();
    // Loop over peaks, looking only at those with a minimum peak value
    for (unsigned ip=0; ip<np; ++ip) {
      Float_t xp = xpeaks[ip];
      Float_t yp = ypeaks[ip];
      TrkTimePeak tpeak(xp,yp);
      if(yp > _ymin){
    // associate all hits in the time window with this peak
	for(size_t istr=0; istr<nstrs;++istr){
	  if(_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_hbkg)){
	    double time = _shcol->at(istr).time();
	    if(fabs(time-xp) < _maxdt){
	      tpeak._trkptrs.push_back(istr);
	    }
	  }
	}
	// if requested, clean the peaks
	if(_cleanpeaks)cleanTimePeak(tpeak);
	if(tpeak._trkptrs.size() >= _minnhits)_tpeaks.push_back(tpeak);
      }
    }
    // sort the peaks so that the largest comes first
    sort(_tpeaks.begin(),_tpeaks.end(),greater<TrkTimePeak>());
    // if requested, fill diagnostics
    if(_diag>1 && _kdiag->mcData()._mcdigis != 0){
      for(size_t ip=0;ip<_tpeaks.size();++ip){
	TrkTimePeak const& tp = _tpeaks[ip];
	fillPeakDiag(ip,tp);
      }
    }
  }

  void TrkPatRec::createTimePeak() {
    // find the median time
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
	double time = _shcol->at(istr).time();
	tacc(time);
      }
    }
    unsigned np = boost::accumulators::extract::count(tacc);  
    if(np >= _minnhits){
      double mtime  = median(tacc);
// create a time peak from the full subset of selected hits
      TrkTimePeak tpeak(mtime,(double)nstrs);
      for(unsigned istr=0; istr<nstrs;++istr){
	if(_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_hbkg)){
	  tpeak._trkptrs.push_back(istr);
	}
      }
      _tpeaks.push_back(tpeak);
    }
  }

  void TrkPatRec::filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec){
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
	if(_kdiag->mcData()._mcdigis != 0){
	  StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(indices[ihit]._index);
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

  void TrkPatRec::findMissingHits(KalFitResult& kalfit,vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_addsel)&& !_flags->at(istr).hasAnyProperty(_addbkg)){
	StrawHit const& sh = _shcol->at(istr);
	if(fabs(_shcol->at(istr).time()-kalfit._krep->t0()._t0) < _maxdtmiss) {
	  // make sure we haven't already used this hit
	  vector<TrkStrawHit*>::iterator ifnd = find_if(kalfit._hits.begin(),kalfit._hits.end(),FindTrkStrawHit(sh));
	  if(ifnd == kalfit._hits.end()){
	    // good in-time hit.  Compute DOCA of the wire to the trajectory
	    Straw const& straw = tracker.getStraw(sh.strawIndex());
	    CLHEP::Hep3Vector hpos = straw.getMidPoint();
	    CLHEP::Hep3Vector hdir = straw.getDirection();
	    // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
	    HepPoint spt(hpos.x(),hpos.y(),hpos.z());
	    TrkLineTraj htraj(spt,hdir,-20,20);
	    // estimate flightlength along track.  This assumes a constant BField!!!
	    double fltlen = (hpos.z()-tpos.z())/tdir.z();
	    TrkPoca hitpoca(kalfit._krep->pieceTraj(),fltlen,htraj,0.0);
	    // flag hits with small residuals
	    if(fabs(hitpoca.doca()) < _maxadddoca){
	      misshits.push_back(istr);
	    }
	  }
	}
      }
    }
  }


  void TrkPatRec::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
    // straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("shpos",&_shp,"x/F:y/F:z/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("deltat",&_deltat,"deltat/F");
    _shdiag->Branch("rho",&_rho,"rho/F");
    _shdiag->Branch("plane",&_plane,"plane/I");
    _shdiag->Branch("panel",&_panel,"panel/I");
    _shdiag->Branch("layer",&_layer,"layer/I");
    _shdiag->Branch("straw",&_straw,"straw/I");
    _shdiag->Branch("mcshpos",&_mcshp,"x/F:y/F:z/F");
    _shdiag->Branch("mcopos",&_mcop,"x/F:y/F:z/F");
    _shdiag->Branch("mcpopos",&_mcpop,"x/F:y/F:z/F");
    _shdiag->Branch("mcoe",&_mcoe,"F");
    _shdiag->Branch("mcom",&_mcom,"F");
    _shdiag->Branch("mcpoe",&_mcpoe,"F");
    _shdiag->Branch("mcpom",&_mcpom,"F");
    _shdiag->Branch("mcshlen",&_mcshlen,"mcshlen/F");
    _shdiag->Branch("mcshd",&_mcshd,"mcshd/F");
    _shdiag->Branch("mcedep",&_mcedep,"mcedep/F");
    _shdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/I");
    _shdiag->Branch("mcnunique",&_mcnunique,"mcnunique/I");
    _shdiag->Branch("mcnmax",&_mcnmax,"mcnmax/I");
    _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
    _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
    _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
    _shdiag->Branch("mctime",&_mctime,"mctime/F");
    _shdiag->Branch("mcppdg",&_mcppdg,"mcpdg/I");
    _shdiag->Branch("mcpproc",&_mcpproc,"mcpproc/I");
    _shdiag->Branch("mcptime",&_mcptime,"mcptime/F");
    _shdiag->Branch("mcgid",&_mcgid,"mcgid/I");
    _shdiag->Branch("mcgpdg",&_mcgpdg,"mcgpdg/I");
    _shdiag->Branch("mcge",&_mcge,"mcge/F");
    _shdiag->Branch("mcgt",&_mcgt,"mcgt/F");
    _shdiag->Branch("mcgpos",&_mcgpos,"x/F:y/F:z/F");
    _shdiag->Branch("esel",&_esel,"esel/I");
    _shdiag->Branch("rsel",&_rsel,"rsel/I");
    _shdiag->Branch("tsel",&_timesel,"tsel/I");
    _shdiag->Branch("delta",&_delta,"delta/I");
    _shdiag->Branch("stereo",&_stereo,"stereo/I");
    _shdiag->Branch("isolated",&_isolated,"isolated/I");
    _shdiag->Branch("pdist",&_pdist,"pdist/F");
    _shdiag->Branch("pperp",&_pperp,"pperp/F");
    _shdiag->Branch("pmom",&_pmom,"pmom/F");
    _shdiag->Branch("pres",&_shpres,"pres/F");
    _shdiag->Branch("rres",&_shrres,"rres/F");
    _shdiag->Branch("chisq",&_shchisq,"chisq/F");
    _shdiag->Branch("dt",&_shdt,"dt/F");
    _shdiag->Branch("dist",&_shdist,"dist/F");
    _shdiag->Branch("xtalk",&_xtalk,"xtalk/B");
    // time peak diagnostics
    _tpdiag=tfs->make<TTree>("tpdiag","time peak diagnostics");
    _tpdiag->Branch("eventid",&_tpeventid,"eventid/I");
    _tpdiag->Branch("peakid",&_peakid,"peakid/I");
    _tpdiag->Branch("pmax",&_pmax,"pmax/I");
    _tpdiag->Branch("nphits",&_nphits,"nphits/I");
    _tpdiag->Branch("ncphits",&_ncphits,"ncphits/I");
    _tpdiag->Branch("nchits",&_nchits,"nchits/I");
    _tpdiag->Branch("ptime",&_ptime,"ptime/F");
    _tpdiag->Branch("ctime",&_ctime,"ctime/F");
    _tpdiag->Branch("pdtimemax",&_pdtimemax,"pdtimemax/F");
    _tpdiag->Branch("cdtimemax",&_cdtimemax,"cdtimemax/F");
    _tpdiag->Branch("pphi",&_pphi,"pphi/F");
    _tpdiag->Branch("cphi",&_cphi,"cphi/F");
    _tpdiag->Branch("cphirange",&_cphirange,"cphirange/F");
    _tpdiag->Branch("pdphimax",&_pdphimax,"pdphimax/F");
    _tpdiag->Branch("cdphimax",&_cdphimax,"cdphimax/F");
    _tpdiag->Branch("tphinfo",&_tphinfo);
 
    // extend the KalDiag track diagnostic tuple
    TTree* trkdiag = _kdiag->createTrkDiag();
    trkdiag->Branch("eventid",&_eventid,"eventid/I");
    trkdiag->Branch("nadd",&_nadd,"nadd/I");
    trkdiag->Branch("ipeak",&_ipeak,"ipeak/I");
    trkdiag->Branch("hcx",&_hcx,"hcx/F");
    trkdiag->Branch("hcy",&_hcy,"hcy/F");
    trkdiag->Branch("hr",&_hr,"hr/F");
    trkdiag->Branch("hdfdz",&_hdfdz,"hdfdz/F");
    trkdiag->Branch("hfz0",&_hfz0,"hfz0/F");
    trkdiag->Branch("mccx",&_mccx,"mccx/F");
    trkdiag->Branch("mccy",&_mccy,"mccy/F");
    trkdiag->Branch("mcr",&_mcr,"mcr/F");
    trkdiag->Branch("mcdfdz",&_mcdfdz,"mcdfdz/F");
    trkdiag->Branch("mcfz0",&_mcfz0,"mcfz0/F");
    trkdiag->Branch("helixfail",&_helixfail,"helixfail/I");
    trkdiag->Branch("seedfail",&_seedfail,"seedfail/I");
    trkdiag->Branch("kalfail",&_kalfail,"kalfail/I");
    trkdiag->Branch("hpar",&_hpar,"hd0/F:hp0/F:hom/F:hz0/F:htd/F");
    trkdiag->Branch("herr",&_hparerr,"hd0err/F:hp0err/F:homerr/F:hz0err/F:htderr/F");
    trkdiag->Branch("spar",&_spar,"sd0/F:sp0/F:som/F:sz0/F:std/F");
    trkdiag->Branch("serr",&_sparerr,"sd0err/F:sp0err/F:somerr/F:sz0err/F:stderr/F");
    trkdiag->Branch("st0",&_st0,"st0/F");
    trkdiag->Branch("snhits",&_snhits,"snhits/I");
    trkdiag->Branch("sndof",&_sndof,"sndof/I");
    trkdiag->Branch("sniter",&_sniter,"sniter/I");
    trkdiag->Branch("snweediter",&_snweediter,"snweediter/I");
    trkdiag->Branch("snactive",&_snactive,"snactive/I");
    trkdiag->Branch("schisq",&_schisq,"schisq/F");
    trkdiag->Branch("nchit",&_nchit,"nchit/I");
    trkdiag->Branch("npeak",&_npeak,"npeak/I");
    trkdiag->Branch("tpeak",&_tpeak,"tpeak/F");
    trkdiag->Branch("nmc",&_nmc,"nmc/I");
    trkdiag->Branch("seedfilt",&_sfilt);
    trkdiag->Branch("helixfilt",&_hfilt);
  }

  void TrkPatRec::fillStrawDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = getTrackerOrThrow();
    _nchit = 0;
    double cfmin = 2*M_PI;
    double cfmax = -2*M_PI;
    accumulator_set<double, stats<tag::mean> > tacc;
    accumulator_set<double, stats<tag::mean> > facc;
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _shcol->at(istr);
      StrawHitPosition const& shp = _shpcol->at(istr);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      _plane = straw.id().getPlane();
      _panel = straw.id().getPanel();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();

      _shp = shp.pos();
      _stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);
      _edep = sh.energyDep();
      _time = sh.time();
      _deltat = sh.dt();
      _rho = shp.pos().perp();
      // summarize the MC truth for this strawhit.  Preset the values in case MC is missing/incomplete
      _mcgid = -1;
      _mcgpdg = -1;
      _mcge = -1.0;
      _mcgt = -1.0;
      _mcgpos = threevec();
      _mcppdg=0;
      _mcpproc=-1;
      _mcptime=0.0;
      _mcpop = threevec();
      _mcpoe = _mcpom = -1.0;
      _mcnmax = -1;
      _mcpdg = -1;
      _mcgen = -1;
      _mcproc = -1;
      _mctime = -1;
      _mcshp = threevec();
      _mcop = threevec();
      _mcoe = -1;
      _mcom = -1;
      _mcshlen = -1;
      _mcshd = -1;
      _mcppdg=0;
      _mcpproc=-1;
      _mcptime=0.0;
      _mcpop = threevec(); 
      _mcpoe = _mcpom = -1.0;
      _xtalk = false;
      if(_kdiag->mcData()._mcdigis != 0){
	StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(istr);
	// use TDC channel 0 to define the MC match
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	CLHEP::Hep3Vector dprod = spmcp->position()-det->toDetector(spp->startPosition());
	static Hep3Vector zdir(0.0,0.0,1.0);
	_pdist = dprod.mag();
	_pperp = dprod.perp(zdir);
	_pmom = spmcp->momentum().mag();
	_mcnunique = mcdigi.stepPointMCs().size();
	// compute energy sum
	_mcedep = mcdigi.energySum();
	_mcnmax = mcdigi.stepPointMCs().size();
	_mcpdg = spp->pdgId();
	_mcproc = spp->creationCode();
	_mcgen = -1;
	if(spp->genParticle().isNonnull())
	  _mcgen = spp->genParticle()->generatorId().id();
	_mctime = spmcp->time();
	_mcshp = spmcp->position();
	_mcop = det->toDetector(spp->startPosition());
	_mcoe = spp->startMomentum().e();
	_mcom = spp->startMomentum().vect().mag();
	_mcshlen = (spmcp->position()-straw.getMidPoint()).dot(straw.getDirection());
	_mcshd = (spmcp->position()-straw.getMidPoint()).dot(straw.getDirection().cross(spmcp->momentum().unit()));
	bool conversion = (_mcpdg == 11 && _mcgen == 2 && spmcp->momentum().mag()>90.0);
	if(conversion){
	  ++_nchit;
  // compute the average time and average phi of conversion hits
	  tacc(sh.time());
	  double phi = shp.pos().phi();
	  if(extract_result<tag::count>(facc) > 0){
	    double dphi = phi - extract_result<tag::mean>(facc);
	    if(dphi > M_PI){
	      phi -= 2*M_PI;
	    } else if(dphi < -M_PI){
	      phi += 2*M_PI;
	    }
	  }
	  facc(phi);
	  if(phi>cfmax)cfmax = phi;
	  if(phi<cfmin)cfmin = phi;
	}
  // immediate parent information
	if(spp.isNonnull() && spp->parent().isNonnull()){
	  const art::Ptr<SimParticle>& psp = spp->parent();
	  _mcppdg = psp->pdgId();
	  _mcpproc = psp->creationCode();
	  _mcptime = psp->startGlobalTime();
	  _mcpop = det->toDetector(psp->startPosition());
	  _mcpoe = psp->startMomentum().e();
	  _mcpom = psp->startMomentum().vect().mag();
	}
// generator information
	if(spp.isNonnull()){
	art::Ptr<SimParticle> sp = spp;
	// find the first parent which comes from a generator
	  while(sp->genParticle().isNull() && sp->parent().isNonnull()){
	    sp = sp->parent();
	  }
	  if(sp->genParticle().isNonnull()){
	    _mcgid = sp->genParticle()->generatorId().id();
	    _mcgpdg = sp->genParticle()->pdgId();
	    _mcge = sp->genParticle()->momentum().e();
	    _mcgt = sp->genParticle()->time();
	    _mcgpos = det->toDetector(sp->genParticle()->position());
	  }
	}
	_xtalk = spmcp->strawIndex() != sh.strawIndex();
      }
      _esel = _flags->at(istr).hasAllProperties(StrawHitFlag::energysel);
      _rsel = _flags->at(istr).hasAllProperties(StrawHitFlag::radsel);
      _timesel = _flags->at(istr).hasAllProperties(StrawHitFlag::timesel);
      _stereo = _flags->at(istr).hasAllProperties(StrawHitFlag::stereo);
      _isolated = _flags->at(istr).hasAllProperties(StrawHitFlag::isolated);
      _delta = _flags->at(istr).hasAllProperties(StrawHitFlag::delta);
      _shpres = _shpcol->at(istr).posRes(StrawHitPosition::phi);
      _shrres = _shpcol->at(istr).posRes(StrawHitPosition::rho);
//  Info depending on stereo hits
      if(_stcol != 0 && _shpcol->at(istr).stereoHitIndex() >= 0){
	_shchisq = _stcol->at(_shpcol->at(istr).stereoHitIndex()).chisq();
	_shdt = _stcol->at(_shpcol->at(istr).stereoHitIndex()).dt();
	_shdist = _stcol->at(_shpcol->at(istr).stereoHitIndex()).dist();
      } else {
	_shchisq = -1.0;
	_shdt = 0.0;
	_shdist = -1.0;
      }
      _shdiag->Fill();
    }
    _nchits = _nchit;
    _cphirange = cfmax-cfmin;
    _cphi = extract_result<tag::mean>(facc);
    _ctime = extract_result<tag::mean>(tacc);
  }

  void TrkPatRec::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp, *ptsp;

    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    char tdsname[100];
    char tpsname[100];
    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(csname,100,"convtspectrum%i",_iev);
    snprintf(tsname,100,"tighttspectrum%i",_iev);
    snprintf(lsname,100,"loosetspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltatspectrum%i",_iev);
    snprintf(tpsname,100,"protontspectrum%i",_iev);
    ttsp = tfs->make<TH1F>(tsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ttsp->SetLineColor(kCyan);
    ltsp = tfs->make<TH1F>(lsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ltsp->SetLineColor(kGreen);
    rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kBlue);
    ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kRed);
    ctsp->SetFillColor(kRed);
    ptsp = tfs->make<TH1F>(tpsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ptsp->SetLineColor(kBlack);
    ptsp->SetFillColor(kBlack);
    tdtsp = tfs->make<TH1F>(tdsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    tdtsp->SetLineColor(kOrange);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _shcol->at(istr).time();
      bool conversion(false);
      bool proton(false);
      // summarize the MC truth for this strawhit
      if(_kdiag->mcData()._mcsteps != 0) {
	StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(istr);
	// use TDC channel 0 to define the MC match
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(StrawDigi::zero)) itdc = StrawDigi::one;
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	int gid(-1);
	if(spp->genParticle().isNonnull())
	  gid = spp->genParticle()->generatorId().id();

	conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	proton = spp->pdgId()==2212;
      }
      // fill plots
      rtsp->Fill(time);
      if(_flags->at(istr).hasAllProperties(_tsel)){
	ttsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
	tdtsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_hbkg)){
	ltsp->Fill(time);
	if(proton)
	  ptsp->Fill(time);
      }
      if(conversion)
	ctsp->Fill(time);

    }
    // plot time peaks
    TList* flist = tdtsp->GetListOfFunctions();
    for(auto ipeak=_tpeaks.begin();ipeak!=_tpeaks.end();++ipeak){
      TMarker* smark = new TMarker(ipeak->_tpeak,ipeak->_peakmax,23);
      smark->SetMarkerColor(kRed);
      smark->SetMarkerSize(1.5);
      flist->Add(smark);
    }
  }

  void TrkPatRec::fillFitDiag(int ipeak,HelixFitResult const& helixfit,
      KalFitResult const& seedfit, KalFitResult const& kalfit) {
    // convenience numbers
    static const double pi(M_PI);
    static const double twopi(2*pi);
    static const double halfpi(0.5*pi);
    // initialize some variables
    _ipeak = ipeak;
    _nmc = 0;
    if(ipeak >= 0){
      const TrkTimePeak& tpeak = _tpeaks[ipeak];
      // time peak information
      _peakmax = tpeak._peakmax;
      _tpeak = tpeak._tpeak;
      _npeak = tpeak._trkptrs.size();
      for(vector<hitIndex>::const_iterator istr= tpeak._trkptrs.begin(); istr != tpeak._trkptrs.end(); ++istr){
	// summarize the MC truth for this strawhit
	if(_kdiag->mcData()._mcsteps != 0) {

	  StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(istr->_index);
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  int gid(-1);
	  if(spp->genParticle().isNonnull())
	    gid = spp->genParticle()->generatorId().id();
	  bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	  if(conversion)
	    ++_nmc;
	}
      } 
    } else {
      _peakmax = -1.0;
      _tpeak = -1.0;
      _npeak = -1;
    }
    // fit status 
    _helixfail = helixfit._fit.failure();
    _seedfail = seedfit._fit.failure();
    _kalfail = kalfit._fit.failure();
    // helix information
    HepVector hpar;
    HepVector hparerr;
    _hfit.helixParams(helixfit,hpar,hparerr);
    _hpar = helixpar(hpar);
    _hparerr = helixpar(hparerr);
    _hcx = helixfit._center.x(); _hcy = helixfit._center.y(); _hr = helixfit._radius;
    _hdfdz = helixfit._dfdz; _hfz0 = helixfit._fz0;
    // seed fit information
    if(seedfit._fit.success()){
      _snhits = seedfit._tdef->strawHitIndices().size();
      _snactive = seedfit._krep->nActive();
      _sniter = seedfit._krep->iterations();
      _sndof = seedfit._krep->nDof();
      _schisq = seedfit._krep->chisq();
      _st0 = seedfit._krep->t0()._t0;
      _snweediter = seedfit._nweediter;
      double loclen;
      const TrkSimpTraj* ltraj = seedfit._krep->localTrajectory(0.0,loclen);
      _spar = helixpar(ltraj->parameters()->parameter());
      _sparerr = helixpar(ltraj->parameters()->covariance());
    } else {
      _snhits = -1;
      _snactive = -1;
      _sniter = -1;
      _sndof = -1;
      _schisq = -1.0;
      _st0 = -1.0;
      _snweediter = -1;
    }
    // use MC truth to define hits and seed helix
    TrkDef mctrk(_shcol,_tpart,_fdir);
    // should be chosing the track ID for conversion a better way, FIXME!!!
    cet::map_vector_key itrk(1);
    if(_kdiag->trkFromMC(itrk,mctrk)){
      // find true center, radius
      double rtrue = fabs(1.0/mctrk.helix().omega());
      double rad = 1.0/mctrk.helix().omega() + mctrk.helix().d0();
      double cx = -rad*sin(mctrk.helix().phi0());
      double cy = rad*cos(mctrk.helix().phi0());
      _mccx = cx; _mccy = cy; _mcr = rtrue;
      _mcdfdz = mctrk.helix().omega()/mctrk.helix().tanDip();
      // fix loop for MC values
      _mcfz0 = -mctrk.helix().z0()*mctrk.helix().omega()/mctrk.helix().tanDip() + mctrk.helix().phi0() - copysign(halfpi,mctrk.helix().omega());
      int nloop = (int)rint((helixfit._fz0 - _mcfz0)/twopi);
      _mcfz0 += nloop*twopi;
    }
    // count # of added hits
    _nadd = 0;
    for(vector<TrkStrawHit*>::const_iterator ish=kalfit._hits.begin();ish!=kalfit._hits.end();++ish){
      if((*ish)->usability()==3)++_nadd;
    }
    // fill kalman fit info.  This needs to be last, as it calls TTree::Fill().
    _kdiag->kalDiag(kalfit._krep);
  }


  void TrkPatRec::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHit const& sh = _shcol->at(ish);
    StrawHitPosition const& shp = _shpcol->at(ish);
    shinfo._pos = shp.pos();
    shinfo._time = sh.time();
    shinfo._rho = shp.pos().perp();
    shinfo._pres = shp.posRes(StrawHitPosition::phi);
    shinfo._rres = shp.posRes(StrawHitPosition::rho);
// info depending on stereo hits
    if(_stcol != 0 && shp.stereoHitIndex() >= 0){
      shinfo._chisq = _stcol->at(shp.stereoHitIndex()).chisq();
      shinfo._stdt = _stcol->at(shp.stereoHitIndex()).dt();
      shinfo._dist = _stcol->at(shp.stereoHitIndex()).dist();
    } else {
      shinfo._chisq = -1.0;
      shinfo._stdt = 0.0;
      shinfo._dist = -1.0;
    }
    shinfo._edep = sh.energyDep();
    const Straw& straw = tracker.getStraw( sh.strawIndex() );
    shinfo._plane = straw.id().getPlane();
    shinfo._panel = straw.id().getPanel();
    shinfo._layer = straw.id().getLayer();
    shinfo._straw = straw.id().getStraw();
    shinfo._esel = shp.flag().hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel = shp.flag().hasAllProperties(StrawHitFlag::radsel);
    shinfo._delta = shp.flag().hasAllProperties(StrawHitFlag::delta);
    shinfo._stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);

    if(_kdiag->mcData()._mcdigis != 0) {

      StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(ish);
      // use TDC channel 0 to define the MC match
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      // create MC info and fill
      TrkStrawHitInfoMC tshinfomc;
      // hit t0 needs correction for event offset FIXME!!!
      shinfo._mct0 = spmcp->time();
      shinfo._mcht = mcdigi.wireEndTime(itdc);
      shinfo._mcpdg = spp->pdgId();
      shinfo._mcproc = spp->creationCode();
      shinfo._mcedep = mcdigi.energySum();
      shinfo._mcgen = -1;
      if(spp->genParticle().isNonnull())
	shinfo._mcgen = spp->genParticle()->generatorId().id();

      shinfo._mcpos = spmcp->position();
      shinfo._mctime = spmcp->time();
      shinfo._mcedep = mcdigi.energySum();;
      shinfo._mcmom = spmcp->momentum().mag();
      double cosd = spmcp->momentum().cosTheta();
      shinfo._mctd = cosd/sqrt(1.0-cosd*cosd);
    }
  }

  void TrkPatRec::fillPeakDiag(size_t ip,TrkTimePeak const& tp) {
    _tpeventid = _eventid;
    _peakid = ip;
    _pmax = tp._peakmax;
    _nphits = tp._trkptrs.size();
    _ncphits = 0;
    _pdtimemax = 0.0;
    _pdphimax = 0.0;
    _cdtimemax = 0.0;
    _cdphimax = 0.0;
    _tphinfo.clear();
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::mean > > tacc;
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHit const& sh = _shcol->at(ish);
      StrawHitPosition const& shp = _shpcol->at(ish);
      double time = sh.time();
      tacc(time);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      if(extract_result<tag::count>(facc) > 0){
	double dphi = phi - extract_result<tag::mean>(facc);
	if(dphi > M_PI){
	  phi -= 2*M_PI;
	} else if(dphi < -M_PI){
	  phi += 2*M_PI;
	}
      }
      facc(phi);
    }
    _pphi =extract_result<tag::mean>(facc);
    _ptime =extract_result<tag::mean>(tacc);
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHitPosition const& shp = _shpcol->at(ish);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      double dphi = phi - _pphi;
      if(dphi > M_PI)
	dphi -= 2*M_PI;
      else if(dphi < -M_PI)
	dphi += 2*M_PI;
      if(fabs(dphi) > _pdphimax)_pdphimax = fabs(dphi);

      double dt = _shcol->at(ish).time() - _ptime;
      if(fabs(dt) > _pdtimemax)_pdtimemax=fabs(dt);

      StrawDigiMC const& mcdigi = _kdiag->mcData()._mcdigis->at(ish);
      // use TDC channel 0 to define the MC match
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      int gid(-1);
      if(spp->genParticle().isNonnull())
	gid = spp->genParticle()->generatorId().id();

      bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);

      if(conversion){
	_ncphits++;
	if(fabs(dt) > _cdtimemax)_cdtimemax = fabs(dt);
      }
      if(conversion && dphi > _cdphimax)_cdphimax = dphi;

      _pmva._dt = dt;
      _pmva._dphi = dphi;
      _pmva._rho = pos.perp();
      double mvaout = _peakMVA->EvaluateMVA(_PMVAType);

      TimePeakHitInfo tph;
      tph._dt = dt;
      tph._dphi = dphi;
      tph._rho = pos.perp();
      tph._mva = mvaout;
      tph._mcpdg = spp->pdgId();
      tph._mcgen = gid;
      tph._mcproc = spp->creationCode();
      _tphinfo.push_back(tph);
    }
    _tpdiag->Fill();
    if(_ncphits > 0.5*_nchits)_icepeak = ip;
  }


  void TrkPatRec::cleanTimePeak(TrkTimePeak& tpeak) {
    static const double twopi(2*M_PI);
    // iteratively filter outliers
    double worstmva(100.0);
    double pphi(1000.0);
    double ptime(-1000.0);
    do{
  // first, compute the average phi and range.  Take care for wrapping
      accumulator_set<double, stats<tag::mean > > facc;
      accumulator_set<double, stats<tag::mean > > tacc;
      for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
	unsigned ish = tpeak._trkptrs[ips]._index;
	double time = _shcol->at(ish).time();
	tacc(time);
	double phi = _shpcol->at(ish).pos().phi();
	if(extract_result<tag::count>(facc) > 0){
	  double dphi = phi - extract_result<tag::mean>(facc);
	  if(dphi > M_PI){
	    phi -= twopi;
	  } else if(dphi < -M_PI){
	    phi += twopi;
	  }
	}
	facc(phi);
      }
      pphi = extract_result<tag::mean>(facc);
      ptime = extract_result<tag::mean>(tacc);
// find the least signal-like hit
      size_t iworst =0;
      worstmva = 100.0;
      for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
	unsigned ish = tpeak._trkptrs[ips]._index;
	double dt = _shcol->at(ish).time() - ptime;
	double rho = _shpcol->at(ish).pos().perp();
	double phi = _shpcol->at(ish).pos().phi();
	double dphi = phi - pphi;
	if(dphi > M_PI){
	  dphi -= twopi;
	} else if(dphi < -M_PI){
	  dphi += twopi;
	}
	// compute MVA
	_pmva._dt = dt;
	_pmva._dphi = dphi;
	_pmva._rho = rho; 
	double mvaout = _peakMVA->EvaluateMVA(_PMVAType);
	if(mvaout < worstmva){
	  worstmva = mvaout;
	  iworst = ips;
	}
      }
      // remove the worst hit
      if(worstmva < _minpeakmva){
	std::swap(tpeak._trkptrs[iworst],tpeak._trkptrs.back());
	tpeak._trkptrs.pop_back();
      }
    } while(tpeak._trkptrs.size() >= _minnhits && worstmva < _minpeakmva);
// final pass: hard cut on dt and dphi
    vector<size_t> toremove;
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::mean > > tacc;
    for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
      unsigned ish = tpeak._trkptrs[ips]._index;
      double dt = _shcol->at(ish).time() - ptime;
      double phi = _shpcol->at(ish).pos().phi();
      double dphi = phi - pphi;
      if(dphi > M_PI){
	dphi -= twopi;
	phi -= twopi;
      } else if(dphi < -M_PI){
	dphi += twopi;
	phi += twopi;
      }
      if(fabs(dt) < _maxpeakdt && fabs(dphi) < _maxpeakdphi){
	tacc(_shcol->at(ish).time());
	facc(phi);
      } else {
	toremove.push_back(ips);
      }
    }
    // actually remove the hits; must start from the back
    std::sort(toremove.begin(),toremove.end(),std::greater<size_t>());
    for(auto irm=toremove.begin();irm!=toremove.end();++irm){
      std::swap(tpeak._trkptrs[*irm],tpeak._trkptrs.back());
      tpeak._trkptrs.pop_back();
    }
    // update peak properties
    pphi = extract_result<tag::mean>(facc);
    ptime = extract_result<tag::mean>(tacc);
    tpeak._tpeak = ptime;
    tpeak._phi = pphi;
  }

  void TrkPatRec::initializeReaders() {
    _peakMVA = new TMVA::Reader();
    _peakMVA->AddVariable("_dt",&_pmva._dt);
    _peakMVA->AddVariable("_dphi",&_pmva._dphi);
    _peakMVA->AddVariable("_rho",&_pmva._rho);
    _peakMVA->BookMVA(_PMVAType,_PMVAWeights);
  }

  void TrkPatRec::findTimePeaks( TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak) { 
    int ibin(1); // root starts bin numbers at 1
    int nbins = tspect.GetNbinsX();
    while(ibin < nbins+1){
      double y = tspect.GetBinContent(ibin);
      if(y >= _1dthresh){
	double xmax = tspect.GetBinCenter(ibin);
	double ymax = y;
	double yprev = y;
	// find contiguous bins above threshold
	int jbin = ibin+1;
	bool descending(false);
	while(jbin < nbins+1 && tspect.GetBinContent(jbin) >= _1dthresh){
	  y =tspect.GetBinContent(jbin);
	  descending |= yprev-y > sqrt(yprev);
	  // don't follow next maximum
	  if(descending && y-yprev > sqrt(y)){
	    break;
	  } else {
	    if(y > ymax){
	      ymax = y;
	      xmax = tspect.GetBinCenter(jbin);
	    }
	    yprev = y;
	    ibin = jbin;
	    ++jbin;
	  }
	}
	xpeak.push_back(xmax);
	ypeak.push_back(ymax);
      }
      ++ibin;
    }
  }
}
using mu2e::TrkPatRec;
DEFINE_ART_MODULE(TrkPatRec);
