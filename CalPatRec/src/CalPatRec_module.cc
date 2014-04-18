// $Id: CalPatRec_module.cc,v 1.8 2014/04/18 16:41:56 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2014/04/18 16:41:56 $
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
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
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
#include "BaBar/BaBar.hh"
#include "BaBar/BbrStringUtils.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/KalFitHack.hh"
// #include "CalPatRec/inc/THackData.hh"

#include "TrkBase/TrkPoca.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
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
// #include "TSpectrum.h"
// #include "TSpectrum2.h"
// #include "TSpectrum3.h"
// #include "TMVA/Reader.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>
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

namespace mu2e {
  class CalPatRec : public art::EDProducer {
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalPatRec(fhicl::ParameterSet const&);
    virtual ~CalPatRec();
    
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce(art::Event& event ); 
    virtual void endJob();

  private:
    unsigned     _iev;
					// configuration parameters
    int          _diag; 
    int          _debug;
    int          _printfreq;
    bool         _addhits; 
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string  _shLabel;
    std::string  _shpLabel;
    std::string  _shfLabel;
    std::string  _ccmLabel; // caloClusterModuleLabel
    std::string  _ccAlgorithm;
    std::string  _ccSeeding;

    std::string  _dtspecpar;

    StrawHitFlag     _tsel, _hsel, _ksel;
    StrawHitFlag     _bkgsel;
    double           _maxedep;
    double           _mindt;
    double           _maxdt;
    double           _maxdtmiss;
    double           _fbf;
					// time spectrum parameters
    bool             _findtpeak;
    unsigned         _maxnpeak;
    unsigned         _minnhits;
    double           _tmin;
    double           _tmax;
    double           _tbin;
    unsigned         _nbins;
    double           _ymin;
    double           _1dthresh;
    double           _minClusterEnergy;	// min seed energy
    double           _pitchAngle;
					// outlier cuts
    double           _maxseeddoca;
    double           _maxhelixdoca;
    double           _maxadddoca;
    double           _maxaddchi;
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search

					// cache of event objects

    const StrawHitCollection*         _shcol;
    const StrawHitFlagCollection*     _shfcol;
    StrawHitFlagCollection*           _flags;
    const StrawHitPositionCollection* _shpcol;
    const CaloClusterCollection*      _ccCollection;

					// Kalman fitters.  Seed fit has a special configuration
    KalFitHack               _seedfit;
    KalFitHack               _kfit;
					// robust helix fitter
    HelixFitHack             _hfit;
					// cache of time peaks
    CalTimePeakCollection*   _tpeaks;
    std::string              _iname;	// data instance name

    PayloadSaver             _payloadSaver;

    std::string              _ccAlgLabel;   // calo cluster algorithm label
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
    void findTimePeaks    (CalTimePeakCollection* TimePeakColl);
    void createTimePeak   (CalTimePeakCollection* TimePeakColl);
    void filterOutliers   (TrkDef& mytrk,Trajectory const& traj,double maxdoca,std::vector<TrkHitFilter>& thfvec);
    void findMissingHits  (KalFitResult& kalfit, std::vector<hitIndex>& indices);
    void createDiagnostics();
    void fillStrawDiag    ();
    void fillTimeDiag     ();
    void fillFitDiag      (int                       ipeak   , 
			   HelixFitHackResult const& helixfit,
			   KalFitResult       const& seedfit ,
			   KalFitResult       const& kalfit  );

    void fillStrawHitInfo (size_t ish, StrawHitInfo& shinfo) const;

    //      void initializeReaders();
					// MC tools
    KalFitMC _kfitmc;
//-----------------------------------------------------------------------------
// strawhit tuple variables
//-----------------------------------------------------------------------------
    TTree*   _shdiag;
    Int_t    _eventid;
    threevec _shp;
    Float_t  _edep;
    Float_t  _time, _rho;
    Int_t    _nmcsteps;
    Int_t    _mcnunique,_mcnmax;
    Int_t    _mcpdg,_mcgen,_mcproc;
    threevec _mcshp, _mcop;
    Float_t  _mcshlen;
    Float_t  _mcedep,_mcemax;
    Float_t  _pdist,_pperp,_pmom;
    Float_t  _mctime;
    Int_t    _esel,_rsel, _timesel,  _delta, _stereo, _isolated;
    Int_t    _device, _sector, _layer, _straw;
    Int_t    _ishpeak, _ntpeak, _nshtpeak;
    //    Float_t  _shtpeak;
    Float_t  _shpres, _shrres, _shchisq, _shdt, _shdist;
    Float_t  _shmct0, _shmcmom, _shmctd;
					// fit tuple variables
    Int_t    _nadd,_ipeak;
    Float_t  _hcx, _hcy, _hr, _hdfdz, _hfz0;
    Float_t  _mccx, _mccy, _mcr, _mcdfdz, _mcfz0;
    Int_t    _helixfail,_seedfail,_kalfail;
    helixpar _hpar,_spar;
    helixpar _hparerr,_sparerr;
    Int_t    _snhits, _snactive, _sniter, _sndof, _snweediter;
    Float_t  _schisq, _st0;
    Int_t    _nchit;
    Int_t    _npeak, _nmc;
    Float_t  _tpeak;
					// hit filtering tuple variables

    std::vector<TrkHitFilter> _sfilt, _hfilt;

					// flow diagnostic
    TH1F* _cutflow;

    //    THackData* fHackData;
  };

  CalPatRec::CalPatRec(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>        ("diagLevel",0)),
    _debug       (pset.get<int>        ("debugLevel",0)),
    _printfreq   (pset.get<int>        ("printFrequency",100)),
    _addhits     (pset.get<bool>       ("addhits",true)),
    _shLabel     (pset.get<std::string>("StrawHitCollectionLabel"        ,"makeSH")),
    _shpLabel    (pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStrawHitPositions")),
    _shfLabel    (pset.get<std::string>("StrawHitFlagCollectionLabel"    ,"FlagBkgHits"   )),
    _ccmLabel    (pset.get<std::string>("caloClusterModuleLabel"         ,"MakeCaloCluster")),
    _ccAlgorithm (pset.get<std::string>("caloClusterAlgorithm"           , "closest")),
    _ccSeeding   (pset.get<std::string>("caloClusterSeeding"             , "energy")),

    _dtspecpar   (pset.get<std::string>("DeltaTSpectrumParams","nobackgroundnomarkovgoff")),
    _tsel        (pset.get<std::vector<std::string> >("TimeSelectionBits")),
    _hsel        (pset.get<std::vector<std::string> >("HelixFitSelectionBits")),
    _ksel        (pset.get<std::vector<std::string> >("KalmanFitSelectionBits")),
    _bkgsel      (pset.get<std::vector<std::string> >("BackgroundSelectionBits")),
    _maxedep     (pset.get<double>("MaxStrawEDep",0.005)),
    _mindt       (pset.get<double>("DtMin",-70.0)),
    _maxdt       (pset.get<double>("DtMax", 20.0)),
    _maxdtmiss   (pset.get<double>("DtMaxMiss",55.0)),
    _fbf         (pset.get<double>("PhiEdgeBuffer",1.1)),
    _findtpeak   (pset.get<bool>("FindTimePeaks",true)),
    _maxnpeak    (pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits    (pset.get<unsigned>("MinNHits" ,20)),
    _tmin            (pset.get<double>("tmin"             ,0.0)),
    _tmax            (pset.get<double>("tmax"             ,2000.0)),
    _tbin            (pset.get<double>("tbin"             ,20.0)),
    _ymin            (pset.get<double>("ymin"             ,4)),
    _1dthresh        (pset.get<double>("OneDPeakThreshold",4.0)),
    _minClusterEnergy(pset.get<double>("minClusterEnergy" ,60.)),
    _pitchAngle      (pset.get<double>("_pitchAngle"      ,0.67)),
    _maxseeddoca     (pset.get<double>("MaxSeedDoca"      ,10.0)),
    _maxhelixdoca    (pset.get<double>("MaxHelixDoca"     ,40.0)),
    _maxadddoca      (pset.get<double>("MaxAddDoca"       ,2.75)),
    _maxaddchi       (pset.get<double>("MaxAddChi"        ,4.0)),
    _tpart           ((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir            ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _seedfit     (pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet())),
    _kfit        (pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("HelixFit",fhicl::ParameterSet())),
    _payloadSaver(pset),
    _kfitmc      (pset.get<fhicl::ParameterSet>("KalFitMC",fhicl::ParameterSet()))
  {
					// tag the data product instance by the direction 
					// and particle type found by this fitter
    _iname = _fdir.name() + _tpart.name();
    produces<KalRepCollection>      (_iname);
    produces<KalRepPtrCollection>   (_iname);
    produces<StrawHitFlagCollection>(_iname);
    produces<CalTimePeakCollection> (_iname);

    produces<KalRepPayloadCollection>();

					// set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);

    _ccAlgLabel = "Algo"+babar::String::toUpper(_ccAlgorithm)+"SeededBy"+
                         babar::String::toUpper(_ccSeeding);
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalPatRec::~CalPatRec() {
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginJob(){
					// create diagnostics ntuple if requested
    if(_diag > 0)createDiagnostics();
					// create a histogram of throughput: this is 
					// a basic diagnostic that should ALWAYS be on

    art::ServiceHandle<art::TFileService> tfs;
    _cutflow = tfs->make<TH1F>("cutflow","Cutflow",10,-0.5,9.5);
    _eventid = 0;
  }

//-----------------------------------------------------------------------------
  void CalPatRec::beginRun(art::Run& ){
  }

//-----------------------------------------------------------------------------
// find the input data objects 
//-----------------------------------------------------------------------------
  bool CalPatRec::findData(const art::Event& evt) {

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitCollection with label=%s not found.\n",
	     _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (evt.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitPositionCollection with label=%s not found.\n",
	     _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalPatRec::findData: StrawHitFlagCollection with label=%s not found.\n",
	     _shfLabel.data());
    }

    art::Handle<CaloClusterCollection> ccH;
    if (evt.getByLabel(_ccmLabel,_ccAlgLabel, ccH)) {
      _ccCollection = ccH.product();
    }
    else {
      _ccCollection = 0;
      printf(" >>> ERROR in CalPatRec::findData: CaloClusterCollection with label=%s:%s not found.\n",
	     _ccmLabel.data(),_ccAlgLabel.data());
    }
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_shpcol != 0) && (_ccCollection != 0);
  }


//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  void CalPatRec::produce(art::Event& event ) {

    bool                      findhelix (false), findseed (false), findkal (false);
    int                       nhits;
    int                       npeaks;
    ::KalRep*                 krep;
					// dummy objects
    static TrkDef             dummydef;
    static HelixDefHack       dummyhdef;

    static HelixFitHackResult dummyhfit(dummyhdef);
    static KalFitResult       dummykfit(dummydef);

    static StrawHitFlag       esel(StrawHitFlag::energysel), flag;

//     static int first_call(1.);
//     if (first_call) {
//       fHackData = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
//     }
					// event printout
    _eventid = event.event();
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) cout<<"CalPatRec: event="<<_iev<<endl;

    _cutflow->Fill(0.0);
					// create output
    _tpeaks = new CalTimePeakCollection;
    unique_ptr<KalRepCollection>       tracks(new KalRepCollection      );
    unique_ptr<KalRepPtrCollection>    trackPtrs(new KalRepPtrCollection);
    unique_ptr<CalTimePeakCollection>  tpeaks(_tpeaks);
    
    _flags = new StrawHitFlagCollection();
    unique_ptr<StrawHitFlagCollection> flags (_flags);

    art::ProductID kalRepsID(getProductID<KalRepCollection>(event,_iname));

					// find the data
    if (!findData(event)) {
      printf("CalPatRec::produce ERROR: No straw hits found, RETURN\n");
                                                            goto END;
    }
//-----------------------------------------------------------------------------
// all needed pieces of data have been found, tighten the energy cut and 
// copy flags
//-----------------------------------------------------------------------------
    nhits = _shcol->size();
    for (int i=0; i<nhits; i++) {
      flag = _shfcol->at(i);
      if (_shcol->at(i).energyDep() > _maxedep && flag.hasAllProperties(esel))
	flag.clear(esel);
      _flags->push_back(flag);
    }
//-----------------------------------------------------------------------------
// find mc truth if we're making diagnostics
//-----------------------------------------------------------------------------
    if (_diag > 0){
      if(!_kfitmc.findMCData(event)){
	cout<<"MC information missing "<< endl;
      }
    }
//-----------------------------------------------------------------------------
// find the time peaks in the time spectrum of selected hits.  
// Otherwise, take all selected hits as a peak
//-----------------------------------------------------------------------------
    if(_findtpeak) {
      findTimePeaks(_tpeaks);
    } else {
      createTimePeak(_tpeaks);
    }
//-----------------------------------------------------------------------------
// diagnostics, MC truth
//-----------------------------------------------------------------------------
    if(_diag>0){
      _kfitmc.mcTrkInfo(_kfitmc.mcData()._simparts->begin()->second);
    }
    if (_diag > 2)fillTimeDiag();
    if (_diag > 1)fillStrawDiag();
    if (_tpeaks->size()>0)_cutflow->Fill(1.0);
//-----------------------------------------------------------------------------
// loop over found time peaks - for us, - "eligible" calorimeter clusters 
//-----------------------------------------------------------------------------
    npeaks = _tpeaks->size();
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      CalTimePeak* tp = &_tpeaks->at(ipeak);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
      if (_debug > 0) {
	const StrawHit*     hit;
	int nh = tp->_trkptrs.size();
	printf(" peak # ipeak = %2i; nhits = %5lu\n",ipeak,tp->_trkptrs.size());
	if (_debug > 1) {
	  for (int ih=0; ih<nh; ih++) {
	    hitIndex ind = tp->_trkptrs[ih];
	    hit = &_shcol->at(ind._index);
	    printf("index = %5i time=%10.3f energy = %10.3f\n",
		   hit->strawIndex().asInt(),hit->time(),hit->energyDep());
	  }
	}
      }
//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information 
//-----------------------------------------------------------------------------
      HelixDefHack helixdef(_shcol,_shpcol,tp->_trkptrs,_tpart,_fdir);

					// set some identifiers
      helixdef.setEventId(_eventid);
      helixdef.setTrackId(ipeak);
					// copy this for the other fits

      TrkDef             seeddef(helixdef);
      TrkDef             kaldef (helixdef);

					// track fitting objects for this peak

      HelixFitHackResult helixfit(helixdef);
      KalFitResult       seedfit (seeddef);
      KalFitResult       kalfit  (kaldef);

					// initialize filters. These are used only for diagnostics
      _hfilt.clear();
      _sfilt.clear();
//-----------------------------------------------------------------------------
// pattern recognition step - find initial approximation for Kalman fitter
//-----------------------------------------------------------------------------
      int rc = _hfit.findHelix(helixfit,tp);
      if (rc) {
//-----------------------------------------------------------------------------
// pattern recognition succeeded
// convert the result to standard helix parameters, and initialize the seed definition helix
//-----------------------------------------------------------------------------
	findhelix = true;
	HepVector hpar;
	HepVector hparerr;
	_hfit.helixParams(helixfit,hpar,hparerr);
	HepSymMatrix hcov = vT_times_v(hparerr);
	seeddef.setHelix(HelixTraj(hpar,hcov));

					// Filter outliers using this helix

	filterOutliers(seeddef,seeddef.helix(),_maxhelixdoca,_hfilt);

					// now, fit the seed helix from the filtered hits

	_seedfit.makeTrack(seedfit, tp);

	if (seedfit._fit.success()) {
	  findseed = true;
//-----------------------------------------------------------------------------
// find the helix parameters from the helix fit, and initialize 
// the full Kalman fit with this
//-----------------------------------------------------------------------------
	  double locflt;
	  const HelixTraj* shelix;
	  shelix = dynamic_cast<const HelixTraj*>(seedfit._krep->localTrajectory(seedfit._krep->flt0(),locflt));
	  kaldef.setHelix(*shelix);
					// filter the outliers
	  filterOutliers(kaldef,seedfit._krep->traj(),_maxseeddoca,_sfilt);

					// 2013-09-16: use HackData - not sure, what for
	  //	  fHackData->fSHelix = shelix;

	  _kfit.makeTrack(kalfit,tp);
//-----------------------------------------------------------------------------
// if successfull, try to add missing hits
//-----------------------------------------------------------------------------
	  if (kalfit._fit.success()) {
	    findkal = true;
	    if (_addhits) {
					// first, add back the hits on this track
	      _kfit.unweedHits(kalfit,_maxaddchi);
	      std::vector<hitIndex> misshits;
	      findMissingHits(kalfit,misshits);
	      if (misshits.size() > 0) {
		_kfit.addHits(kalfit,_shcol,misshits,_maxaddchi);
	      }
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
//  fill fit diagnostics histograms if requested
//-----------------------------------------------------------------------------
      if (_diag > 0) fillFitDiag(ipeak,helixfit,seedfit,kalfit);

      if (kalfit._fit.success()) {
					// flag hits used in this track.  
					// This should use the track id, FIXME!!! (in the BaBar code)
	if (ipeak<16) {
	  for (size_t ihit=0; ihit<kalfit._hits.size(); ++ihit){
	    const TrkStrawHit* tsh = kalfit._hits[ihit];
	    if(tsh->isActive())_flags->at(tsh->index()).merge(StrawHitFlag::trackBit(ipeak));
	  }
	}
//-----------------------------------------------------------------------------
// save successful kalman fits in the event
//-----------------------------------------------------------------------------
	krep = kalfit.stealTrack();
	tracks->push_back(krep);
        int index = tracks->size()-1;
        trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	tp->SetCprIndex(tracks->size());
      } 
      else {
	kalfit.deleteTrack();
      }
//-----------------------------------------------------------------------------
// cleanup the seed fit
//-----------------------------------------------------------------------------
      seedfit.deleteTrack();
    }
//-----------------------------------------------------------------------------
// diagnostics in the end
//-----------------------------------------------------------------------------
    if (findhelix) _cutflow->Fill(2.0);
    if (findseed ) _cutflow->Fill(3.0);
    if (findkal  ) _cutflow->Fill(4.0);

    // add a dummy entry in case there are no peaks
    if(_diag > 0 && _tpeaks->size() == 0)
      fillFitDiag(-1,dummyhfit,dummykfit,dummykfit);
//-----------------------------------------------------------------------------
// put tracks into the event
//-----------------------------------------------------------------------------
  END:;
    art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*tracks, tracksID, event);
    event.put(std::move(tracks),   _iname);
    event.put(std::move(trackPtrs),_iname);
    event.put(std::move(flags ),   _iname);
    event.put(std::move(tpeaks),   _iname);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  void CalPatRec::findTimePeaks(CalTimePeakCollection* TimePeakColl) {

    int                 ncl, nsh;
    double              time, dt, tof, zstraw, cl_time, stime; 
    double              xcl, ycl, zcl, dz_cl;
    const CaloCluster*  cl;
    const StrawHit*     hit;
    const Straw*        straw;
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
//     mu2e::GeomHandle<mu2e::DiskCalorimeter> ch;
//     const mu2e::DiskCalorimeter* cal = ch.get();

    mu2e::GeomHandle<mu2e::TTracker> ttH;
    const mu2e::TTracker* tracker = ttH.get();
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nsh = _shcol->size();
    ncl = _ccCollection->size();

    for (int ic=0; ic<ncl; ic++) {
      cl      = &_ccCollection->at(ic);

      if (cl->energyDep() > _minClusterEnergy) {
	cl_time = cl->time();
	xcl     = cl->cog3Vector().x()+3904.;
	ycl     = cl->cog3Vector().y();
	zcl     = cl->cog3Vector().z(); // cal->disk(cl->vaneId()).origin().z();

	dz_cl   = zcl-tracker->z0();
	CalTimePeak tpeak(cl,xcl,ycl,dz_cl);
	tpeak._tmin = cl_time+_mindt;
	tpeak._tmax = cl_time+_maxdt;
//-----------------------------------------------------------------------------
// record hits in time with each peak, and accept them if they have a minimum # of hits
//-----------------------------------------------------------------------------
	stime = 0;
	for(int istr=0; istr<nsh;++istr) {
	  if (_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)) {
	    hit    = &_shcol->at(istr);
	    time   = hit->time();
	    straw  = &tracker->getStraw(hit->strawIndex());
	    zstraw = straw->getMidPoint().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
//-----------------------------------------------------------------------------
	    tof = (dz_cl-zstraw)/sin(_pitchAngle)/CLHEP::c_light;
	    dt  = cl_time-(time+tof);

	    if ((dt < _maxdt) && (dt >= _mindt)) {
	      tpeak._trkptrs.push_back(istr);
	      stime += time;
	    }
	  }
	}

	tpeak._tpeak = stime/(tpeak._trkptrs.size()+1.e-12);
	if (tpeak._trkptrs.size() > _minnhits) TimePeakColl->push_back(tpeak);
      }
    }
  }

//-----------------------------------------------------------------------------
// 2014-04-08 P.M.: I don't think this function is called any more
//-----------------------------------------------------------------------------
  void CalPatRec::createTimePeak(CalTimePeakCollection* TimePeakColl) {
// find the median time
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
    unsigned nstrs = _shcol->size();
    double   time;

    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)) {
	time = _shcol->at(istr).time();
	tacc(time);
      }
    }

    unsigned np = boost::accumulators::extract::count(tacc);  
    if(np >= _minnhits){
      double mtime  = median(tacc);
      // create a time peak from the full subset of selected hits
      CalTimePeak tpeak(0, 0., 0., 0.);
      for(unsigned istr=0; istr<nstrs;++istr){
	if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	  tpeak._trkptrs.push_back(istr);
	}
      }
      tpeak._tpeak = mtime;
      TimePeakColl->push_back(tpeak);
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::filterOutliers(TrkDef&                    mytrk  , 
				 Trajectory const&          traj   , 
				 double                     maxdoca, 
				 std::vector<TrkHitFilter>& thfvec ) {
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;
    traj.getInfo(0.0,tpos,tdir);

    // tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const std::vector<hitIndex>& indices = mytrk.strawHitIndices();
    std::vector<hitIndex> goodhits;

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
	if(_kfitmc.mcData()._mcsteps != 0){
	  const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ihit);
	  thfilter._mcpdg = mcsum[0]._pdgid;
	  thfilter._mcgen = mcsum[0]._gid;
	  thfilter._mcproc = mcsum[0]._pid;
	}
	thfvec.push_back(thfilter);
      }
    }
    // update track
    mytrk.setIndices(goodhits);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::findMissingHits(KalFitResult& kalfit,std::vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
					//  Trajectory info
    Hep3Vector tdir;
    HepPoint   tpos;

    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_ksel)&& !_flags->at(istr).hasAnyProperty(_bkgsel)){
	StrawHit const& sh = _shcol->at(istr);
	if(fabs(_shcol->at(istr).time()-kalfit._krep->t0()._t0) < _maxdtmiss) {
	  // make sure we haven't already used this hit
	  std::vector<TrkStrawHit*>::iterator ifnd = find_if(kalfit._hits.begin(),kalfit._hits.end(),FindTrkStrawHit(sh));
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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
//-----------------------------------------------------------------------------
// straw hit tuple
//-----------------------------------------------------------------------------
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("shpos",&_shp,"x/F:y/F:z/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("rho",&_rho,"rho/F");
    _shdiag->Branch("device",&_device,"device/I");
    _shdiag->Branch("sector",&_sector,"sector/I");
    _shdiag->Branch("layer",&_layer,"layer/I");
    _shdiag->Branch("straw",&_straw,"straw/I");
    _shdiag->Branch("ishpeak",&_ishpeak,"ishpeak/I");
    _shdiag->Branch("ntpeak",&_ntpeak,"ntpeak/I");
    //    _shdiag->Branch("tpeak",&_shtpeak,"tpeak/F");
    _shdiag->Branch("nshtpeak",&_nshtpeak,"nshtpeak/I");
    _shdiag->Branch("mcshpos",&_mcshp,"x/F:y/F:z/F");
    _shdiag->Branch("mcopos",&_mcop,"x/F:y/F:z/F");
    _shdiag->Branch("mcshlen",&_mcshlen,"mcshlen/F");
    _shdiag->Branch("mcedep",&_mcedep,"mcedep/F");
    _shdiag->Branch("mcemax",&_mcemax,"mcemax/F");
    _shdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/I");
    _shdiag->Branch("mcnunique",&_mcnunique,"mcnunique/I");
    _shdiag->Branch("mcnmax",&_mcnmax,"mcnmax/I");
    _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
    _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
    _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
    _shdiag->Branch("mctime",&_mctime,"mctime/F");
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
    _shdiag->Branch("mct0",&_shmct0,"mct0/F");
    _shdiag->Branch("mcmom",&_shmcmom,"mcmom/F");
    _shdiag->Branch("mctd",&_shmctd,"mctd/F");
//-----------------------------------------------------------------------------
// extend the KalFitMC track diagnostic tuple
//-----------------------------------------------------------------------------
    TTree* trkdiag = _kfitmc.createTrkDiag();

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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::fillStrawDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = getTrackerOrThrow();
    _nchit = 0;
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _shcol->at(istr);
      StrawHitPosition const& shp = _shpcol->at(istr);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      _device = straw.id().getDevice();
      _sector = straw.id().getSector();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();

      _shp = shp.pos();
      _stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);
      _edep = sh.energyDep();
      _time = sh.time();
      _rho = shp.pos().perp();
      // find proximity for different radii
      double esum(0.0);
      // MC information
      //      StrawHitMCTruth const& mcstrawhit = (_kfitmc.mcData()._mcstrawhits->at(istr));
      PtrStepPointMCVector const& mcptr(_kfitmc.mcData()._mchitptr->at(istr));
      // compute weighted distance from particle production
      _pdist = 0.0;
      _pperp = 0.0;
      _pmom = 0.0;
      _nmcsteps = mcptr.size();
      for( size_t imc=0; imc< mcptr.size(); ++imc ) {
	StepPointMC const& mchit = *mcptr[imc];
	// distance from production
	double edep = mchit.eDep();
	esum += edep;
	CLHEP::Hep3Vector dprod = mchit.position()-det->toDetector(mchit.simParticle()->startPosition());
	_pdist += dprod.mag()*edep;
	static Hep3Vector zdir(0.0,0.0,1.0);
	_pperp += dprod.perp(zdir)*edep;
	_pmom += mchit.momentum().mag()*edep;
      }
      if(esum > 0.0){
	_pdist /= esum;
	_pperp /= esum;
	_pmom /= esum;
      }
      // summarize the MC truth for this strawhit
      if(_kfitmc.mcData()._mcsteps != 0){
	const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
	_mcnunique = mcsum.size();
	// compute energy sum
	_mcedep = 0.0;
	for(std::vector<MCHitSum>::const_iterator isum=mcsum.begin(); isum != mcsum.end(); ++isum){
	  _mcedep += isum->_esum;
	}
	// first entry
	_mcemax = mcsum[0]._esum;
	_mcnmax = mcsum[0]._count;
	_mcpdg = mcsum[0]._pdgid;
	_mcgen = mcsum[0]._gid;
	_mcproc = mcsum[0]._pid;
	_mctime = mcsum[0]._time;
	_mcshp = mcsum[0]._pos;
	_mcop = det->toDetector(mcsum[0]._spp->startPosition());
	_mcshlen = (mcsum[0]._pos-straw.getMidPoint()).dot(straw.getDirection());
	bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
	if(conversion){
	  ++_nchit;
	}
      }
      _esel     = _flags->at(istr).hasAllProperties(StrawHitFlag::energysel);
      _rsel     = _flags->at(istr).hasAllProperties(StrawHitFlag::radsel);
      _timesel  = _flags->at(istr).hasAllProperties(StrawHitFlag::timesel);
      _stereo   = _flags->at(istr).hasAllProperties(StrawHitFlag::stereo);
      _isolated = _flags->at(istr).hasAllProperties(StrawHitFlag::isolated);
      _delta    = _flags->at(istr).hasAllProperties(StrawHitFlag::delta);
      _shpres   = _shpcol->at(istr).posRes(StrawHitPosition::phi);
      _shrres   = _shpcol->at(istr).posRes(StrawHitPosition::rho);
      _shmct0   = _kfitmc.MCT0(KalFitMC::trackerMid);
      _shmcmom  = _kfitmc.MCMom(KalFitMC::trackerMid);
      _shmctd   = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;

      _shchisq  = -1.0;
      _shdt     = 0.0;
      _shdist   = -1.0;

      // compare to different time peaks
      _ntpeak = _tpeaks->size();
      //      _nshtpeak = 0;
      //      _shtpeak = -1.0;
      //      _ishpeak = -1;
      //      hitIndex myindex(istr);
//       if(_shmcmom >0){
// 	for(unsigned ipeak=0;ipeak<_tpeaks->size();++ipeak){
// 	  std::vector<hitIndex>::iterator ifind =
// 	    std::find(_tpeaks[ipeak]._trkptrs.begin(),_tpeaks->at(ipeak)._trkptrs.end(),myindex);
// 	  if(ifind != _tpeaks[ipeak]._trkptrs.end()){
// 	    _ishpeak = ipeak;
// 	    break;
// 	  }
// 	}
//       }
//      if(_ishpeak>=0){
//	_nshtpeak = _tpeaks[_ishpeak]._trkptrs.size();
	//	_shtpeak = _tpeaks[_ishpeak]._tpeak;
      //      }
      _shdiag->Fill();
    }
  }

  void CalPatRec::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp;

    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    char tdsname[100];
    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(csname,100,"convtspectrum%i",_iev);
    snprintf(tsname,100,"tighttspectrum%i",_iev);
    snprintf(lsname,100,"loosetspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltatspectrum%i",_iev);
    ttsp = tfs->make<TH1F>(tsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ttsp->SetLineColor(kCyan);
    ltsp = tfs->make<TH1F>(lsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ltsp->SetLineColor(kGreen);
    rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kBlue);
    ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kRed);
    tdtsp = tfs->make<TH1F>(tdsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    tdtsp->SetLineColor(kOrange);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _shcol->at(istr).time();
      bool conversion(false);
      // summarize the MC truth for this strawhit
      if(_kfitmc.mcData()._mcsteps != 0) {
	const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
	conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
      }
      // fill plots
      rtsp->Fill(time);
      if(_flags->at(istr).hasAllProperties(_tsel)){
	ttsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	tdtsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_ksel) && !_flags->at(istr).hasAnyProperty(_bkgsel)){
	ltsp->Fill(time);
      }
      if(conversion){
	ctsp->Fill(time);
      }
    }
// 					// find peaks, so they show up on diagnostic plot too
//     TSpectrum tspec(_maxnpeak);
//     Double_t mb = tdtsp->GetMaximum();
//     double thresh(0.1);
//     if(mb > _1dthresh) thresh = _1dthresh/mb;
//     tspec.Search(tdtsp,1,_dtspecpar.c_str(),thresh);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalPatRec::fillFitDiag (int                       ipeak   , 
			       HelixFitHackResult const& helixfit,
			       KalFitResult const&       seedfit , 
			       KalFitResult const&       kalfit  ) {
    // convenience numbers
    static const double pi     (M_PI);
    static const double twopi  (2*pi);
    static const double halfpi (0.5*pi);
    // initialize some variables
    _ipeak = ipeak;
    _nmc = 0;
    if(ipeak >= 0){
      CalTimePeak& tpeak = _tpeaks->at(ipeak);

					// time peak information
      _tpeak   = tpeak._tpeak;
      _npeak   = tpeak._trkptrs.size();
      for(std::vector<hitIndex>::const_iterator istr= tpeak._trkptrs.begin(); istr != tpeak._trkptrs.end(); ++istr){
	// summarize the MC truth for this strawhit
	if(_kfitmc.mcData()._mcsteps != 0) {
	  const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr->_index); 
	  if(mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2)
	    ++_nmc;
	}
      } 
    } else {
      _tpeak   = -1.0;
      _npeak   = -1;
    }
    // fit status 
    _helixfail = helixfit._fit.failure();
    _seedfail  = seedfit._fit.failure();
    _kalfail   = kalfit._fit.failure();
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
      _snhits = seedfit._tdef.strawHitIndices().size();
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
    if(_kfitmc.trkFromMC(itrk,mctrk)){
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
    for(std::vector<TrkStrawHit*>::const_iterator ish=kalfit._hits.begin();ish!=kalfit._hits.end();++ish){
      if((*ish)->usability()==3)++_nadd;
    }
    // fill kalman fit info.  This needs to be last, as it calls TTree::Fill().
    _kfitmc.kalDiag(kalfit._krep);
  }

//-----------------------------------------------------------------------------
// 2014-04-03 P.M.: as far as I understand, this routine is not being used
//-----------------------------------------------------------------------------
  void CalPatRec::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {

    const Tracker& tracker = getTrackerOrThrow();

    StrawHit const& sh          = _shcol->at(ish);
    StrawHitPosition const& shp = _shpcol->at(ish);

    shinfo._pos   = shp.pos();
    shinfo._time  = sh.time();
    shinfo._rho   = shp.pos().perp();
    shinfo._pres  = shp.posRes(StrawHitPosition::phi);
    shinfo._rres  = shp.posRes(StrawHitPosition::rho);

    shinfo._chisq = -1.0;
    shinfo._stdt  = 0.0;
    shinfo._dist  = -1.0;

    shinfo._edep  = sh.energyDep();
    const Straw& straw = tracker.getStraw( sh.strawIndex() );

    shinfo._device = straw.id().getDevice();
    shinfo._sector = straw.id().getSector();
    shinfo._layer  = straw.id().getLayer();
    shinfo._straw  = straw.id().getStraw();
    shinfo._esel   = shp.flag().hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel   = shp.flag().hasAllProperties(StrawHitFlag::radsel);
    shinfo._delta  = shp.flag().hasAllProperties(StrawHitFlag::delta);
    shinfo._stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);

    if(_kfitmc.mcData()._mcsteps != 0) {
      const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      shinfo._mcpdg  = mcsum[0]._pdgid;
      shinfo._mcgen  = mcsum[0]._gid;
      shinfo._mcproc = mcsum[0]._pid;
      shinfo._mcpos  = mcsum[0]._pos;
      shinfo._mctime = mcsum[0]._time;
      shinfo._mcedep = mcsum[0]._esum;
      shinfo._mct0   = _kfitmc.MCT0(KalFitMC::trackerMid);
      shinfo._mcmom  = _kfitmc.MCMom(KalFitMC::trackerMid);
      shinfo._mctd   = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
    }
  }
}
using mu2e::CalPatRec;
DEFINE_ART_MODULE(CalPatRec);
