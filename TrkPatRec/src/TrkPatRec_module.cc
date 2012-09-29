//
// $Id: TrkPatRec_module.cc,v 1.43 2012/09/29 18:30:07 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/09/29 18:30:07 $
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
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/DeltaHitInfo.hh"
#include "TrkPatRec/inc/HelixFit.hh"
#include "TrkBase/TrkPoca.hh"
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
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TMVA/Reader.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
using namespace std; 

namespace mu2e 
{
// structs for MVAs
  struct DeltaHitMVA {
    Float_t _ddphi; // delta-phi
    Float_t _ddrho; // delta-rho
    Float_t _ddt;  // delta time
  };

  struct DeltaPeakMVA {
    Float_t _rfold;
    Float_t _zmin, _zmax, _zgap;
    Float_t _sphi;
    Float_t _nsmiss, _nsh;
  };

  struct DeltaInfo {
// peak information
    double _tpeak, _ppeak;
    double _mindist, _mindt, _mindphi;
// median information
    double _tmed, _pmed, _rmed;
    double _tmean, _pmean, _rmean;
    double _trms, _prms, _rrms;
//summary information
    unsigned _ngoodhits;
    double _zmin, _zmax, _zgap;
    unsigned _ismin, _ismax, _ns, _nsmiss;
    bool _isdelta, _oldisdelta;
    double _pmvaout;
// information about the hits in the delta
    std::vector<size_t> _hindex;
    std::vector<double> _htime;
    std::vector<double> _hphi;
    std::vector<double> _hrho;
    std::vector<double> _hz;
    std::vector<double> _hgd;
    std::vector<int> _hflag;
  };

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
    // configuration parameters
    int _diag,_debug;
    int _printfreq;
    bool _addhits;
    // event object labels
    std::string _strawhitslabel;
    // cut variables
    double _edept, _edepl, _edepvl;
    double _rmint, _rminl;
    double _rmaxt, _rmaxl;
    double _maxdt, _maxdtmiss;
  // delta-ray removal parameters
    bool _filterdeltas;
    double _max2ddt,_maxdp;
    unsigned _maxndelta, _npbins, _ntbins;
    int _nselbins;
    double _2dthresh, _2dsigma;
    double _fbf,_mindp;
    double _maxzgap,_maxnsmiss;
    double _minrho, _maxrho;
    std::string _dhittype, _dhitweights;
    std::string _dpeaktype, _dpeakweights;
    double _dhitmvacut, _dpeakmvacut;
    // time spectrum parameters
    unsigned _maxnpeak;
    unsigned _minnhits;
    double _tmin;
    double _tmax;
    double _tbin;
    unsigned _nbins;
    double _ymin;
    double _1dthresh;
    double _tpeakerr;
    // outlier cuts
    double _maxseeddoca,_maxhelixdoca,_maxadddoca, _maxaddchi;
    TrkParticle _tpart; // particle type being searched for
    TrkFitDirection _fdir;  // fit direction in search
    // cache of event objects
    const StrawHitCollection* _strawhits;
   // Kalman fitters.  Seed fit has a special configuration
    KalFit _seedfit, _kfit;
  // robust helix fitter
    HelixFit _hfit;
  // cache of hit positions
    std::vector<CLHEP::Hep3Vector> _shpos;
  // cache of hit falgs
    std::vector<TrkHitFlag> _tflags;
  // cache of time peaks
    std::vector<TrkTimePeak> _tpeaks;
    std::string _iname; // data instance name
    //
    PayloadSaver _payloadSaver;
   // helper functions
    bool findData(const art::Event& e);
    bool tighthit(double edep, double rho);
    bool loosehit(double edep, double rho);
    bool veryloosehit(double edep, double rho);
    void findProximity( unsigned ish, double maxdist, unsigned& nprox, double& dmin);
    void findPositions();
    void preselectHits();
    void filterDeltas();
    void fillDeltaInfo(TSpectrum2 const& tspec2,std::vector<DeltaInfo>& dinfo);
    void fillDeltaSummary(DeltaInfo& delta);
    void fillDeltaDiag(std::vector<DeltaInfo> const& deltas);
    void findTimePeaks();
    void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,std::vector<TrkHitFilter>& thfvec);
    void findMissingHits(KalFitResult& kalfit, std::vector<hitIndex>& indices);
    void createDiagnostics();
    void fillStrawDiag();
    void fillTimeDiag();
    void fillFitDiag(int ipeak, HelixFitResult const& helixfit,
	KalFitResult const& seedfit,KalFitResult const& kalfit);
    void fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const;
    double correctedTime(size_t istr) const;
    void initializeReaders();
// MC tools
    KalFitMC _kfitmc;
// strawhit tuple variables
    TTree* _shdiag;
    Int_t _eventid;
    threevec _shp;
    Float_t _edep;
    Float_t _time, _corrtime;
    Float_t _dmin;
    Int_t _n50,_n100,_n150,_n200;
    Int_t _nmcsteps;
    Int_t _mcnunique,_mcnmax;
    Int_t _mcpdg,_mcgen,_mcproc;
    threevec _mcshp;
    Float_t _mcedep,_mcemax;
    Float_t _pdist,_pperp,_pmom;
    Float_t _mctime;
    Int_t _vloose, _loose, _tight, _delta;
    Int_t _device, _sector, _layer, _straw;
    Int_t _ishpeak, _ntpeak, _nshtpeak;
    Float_t _shtpeak;
    Float_t _shmct0, _shmcmom, _shmctd;
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
    std::vector<TrkHitFilter> _sfilt, _hfilt;
// delta removal diagnostics
    TTree* _ddiag;
    TMVA::Reader *_dhReader; // assign hits to delta peaks
    DeltaHitMVA _dhmva; // input variables to TMVA for delta hit selection
    DeltaPeakMVA _dpmva;
    TMVA::Reader* _dpReader; // classify delta peaks
    Int_t _ip, _iev;
    Bool_t _isdelta, _oldisdelta;
    Float_t _pphi, _pt, _prho;
    Float_t _zmin, _zmax, _zgap;
    Int_t _ns, _smin, _smax, _nsmiss;
    Int_t _nsh, _ndpeak, _ndmax; 
    Int_t _nconv, _ndelta, _ncompt, _ngconv, _nebkg, _nprot;
    std::vector<DeltaHitInfo> _phits;
    Float_t _dmct0, _dmcmom, _dmctd;
    Float_t _mindist, _mindt, _mindphi;
    Float_t _tmed, _pmed, _rmed;
    Float_t _tmean, _pmean, _rmean;
    Float_t _stime, _sphi, _srho;
    Float_t _pmvaout;
// flow diagnostic
    TH1F* _cutflow;
 };

  TrkPatRec::TrkPatRec(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel","makeSH")),
    _edept(pset.get<double>("EDep_tight",0.0045)),
    _edepl(pset.get<double>("EDep_loose",0.005)),
    _edepvl(pset.get<double>("EDep_veryloose",0.008)),
    _rmint(pset.get<double>("RMin_tight",420.0)),
    _rminl(pset.get<double>("RMin_loose",390.0)),
    _rmaxt(pset.get<double>("RMax_tight",630.0)),
    _rmaxl(pset.get<double>("RMax_loose",650.0)),
    _maxdt(pset.get<double>("DtMax",35.0)),
    _maxdtmiss(pset.get<double>("DtMaxMiss",55.0)),
    _filterdeltas(pset.get<bool>("FilterDeltas",true)),
    _max2ddt(pset.get<double>("Dt2DMax",75.0)),
    _maxdp(pset.get<double>("DPhiMax",0.21)),
    _maxndelta(pset.get<unsigned>("MaxNDeltas",200)),
    _npbins(pset.get<unsigned>("NPhiBins",100)),
    _ntbins(pset.get<unsigned>("NTimeBins",80)),
    _nselbins(pset.get<int>("NSelBins",3)),
    _2dthresh(pset.get<double>("TwoDPeakThreshold",3)),
    _2dsigma(pset.get<double>("TwoDPeakSigma",1.0)),
    _fbf(pset.get<double>("PhiEdgeBuffer",1.1)),
    _mindp(pset.get<double>("Min2dPeak",5)),
    _maxzgap(pset.get<double>("MaxZGap",0.0)),
    _maxnsmiss(pset.get<double>("MaxNMiss",4)),
    _minrho(pset.get<double>("MinRho",410.0)),
    _maxrho(pset.get<double>("MaxRho",660.0)),
    _dhittype(pset.get<std::string>("DeltaHitTMVAType","BDT method")),
    _dhitweights(pset.get<std::string>("DeltaHitTMVAWeights","TrkPatRec/test/deltahits_BDT.weights.xml")),
    _dpeaktype(pset.get<std::string>("DeltaPeakTMVAType","BDT method")),
    _dpeakweights(pset.get<std::string>("DeltaPeakTMVAWeights","TrkPatRec/test/deltapeak_BDT.weights.xml")),
    _dhitmvacut(pset.get<double>("DeltaHitMVACut",0.4)),
    _dpeakmvacut(pset.get<double>("DeltaPeakMVACut",0.1)),
    _maxnpeak(pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits(pset.get<unsigned>("MinNHits",0)),
    _tmin(pset.get<double>("tmin",0.0)),
    _tmax(pset.get<double>("tmax",2000.0)),
    _tbin(pset.get<double>("tbin",20.0)),
    _ymin(pset.get<double>("ymin",5)),
    _1dthresh(pset.get<double>("OneDPeakThreshold",5.0)),
    _tpeakerr(pset.get<double>("timepeakerr",-8.0)),
    _maxseeddoca(pset.get<double>("MaxSeedDoca",10.0)),
    _maxhelixdoca(pset.get<double>("MaxHelixDoca",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",5.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet())),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet())),
    _hfit(pset.get<fhicl::ParameterSet>("HelixFit",fhicl::ParameterSet())),
    _payloadSaver(pset),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC",fhicl::ParameterSet()))
  {
// tag the data product instance by the direction and particle type found by this fitter
    _iname = _fdir.name() + _tpart.name();
    produces<KalRepCollection>(_iname);
    produces<KalRepPayloadCollection>();
// set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
  }

  TrkPatRec::~TrkPatRec(){}

  void TrkPatRec::beginJob(){
// create diagnostics if requested
    if(_diag > 0)createDiagnostics();
// create a histogram of throughput: this is a basic diagnostic that should ALWAYS be on
    art::ServiceHandle<art::TFileService> tfs;
    _cutflow=tfs->make<TH1F>("cutflow","Cutflow",10,-0.5,9.5);
    _eventid = 0;
// initialize the TMVA readers
    initializeReaders();
  }

  void TrkPatRec::beginRun(art::Run& ){}

  void TrkPatRec::produce(art::Event& event ) {
    ++_eventid;
    _cutflow->Fill(0.0);
// create output
    auto_ptr<KalRepCollection> tracks(new KalRepCollection );
// event printout
    _iev=event.id().event();
    if((_iev%_printfreq)==0)cout<<"TrkPatRec: event="<<_iev<<endl;
// find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
// find mc truth if we're making diagnostics
    if(_diag > 0){
      if(!_kfitmc.findMCData(event)){
	cout<<"MC information missing "<< endl;
//	return;
      }
    }
//  find hit positions.  This uses conditions data, so it's not an attribute of the hits 
    findPositions();
// preselect hits based on internal properties (position, energy, ...)
    preselectHits();
// filter 'delta rays'
    if(_filterdeltas)filterDeltas(); 
// find the time peaks in the time spectrum of selected hits
    findTimePeaks();
// fill diagnostics if requested
    if(_diag > 2)fillTimeDiag();
    if(_diag > 1)fillStrawDiag();
// dummy objects
    static TrkDef dummydef;
    static HelixFitResult dummyhfit(dummydef);
    static KalFitResult dummykfit(dummydef);
// loop over the accepted time peaks
    if(_tpeaks.size()>0)_cutflow->Fill(1.0);
    bool findhelix(false), findseed(false), findkal(false);
    for(unsigned ipeak=0;ipeak<_tpeaks.size();++ipeak){
// create track definitions for the helix fit from this initial information 
      TrkDef helixdef(_strawhits,_tpeaks[ipeak]._trkptrs,_tpart,_fdir);
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
      if(_hfit.findHelix(helixfit)){
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
	      std::vector<hitIndex> misshits;
	      findMissingHits(kalfit,misshits);
	      if(misshits.size() > 0){
		_kfit.addHits(kalfit,_strawhits,misshits,_maxaddchi);
	      }
	    }
	  }
        }
      }
// fill fit diagnostics if requested
      if(_diag > 0)
	fillFitDiag(ipeak,helixfit,seedfit,kalfit);
      if(kalfit._fit.success()){
// save successful kalman fits in the event
	tracks->push_back( kalfit.stealTrack() );
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
    event.put(tracks,_iname);
  }

  void TrkPatRec::endJob()
  {
// does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

// find the input data objects 
  bool TrkPatRec::findData(const art::Event& evt){
    _strawhits = 0;
    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(evt.getByLabel(_strawhitslabel,strawhitsH))
      _strawhits = strawhitsH.product();
    return _strawhits != 0;
  }

  bool TrkPatRec::tighthit(double edep, double rho){
    return edep < _edept && rho > _rmint && rho < _rmaxt;
  }

  bool TrkPatRec::veryloosehit(double edep, double rho){
  // very loose cuts for adding hits to an existing track
    return edep < _edepvl;
  }

  bool TrkPatRec::loosehit(double edep, double rho){
  // looser cuts for pat. rec.
    return edep < _edepl && rho > _rminl && rho < _rmaxl; 
  }

  void
  TrkPatRec::findPositions(){
    _shpos.clear();
// tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      const CLHEP::Hep3Vector& mid = straw.getMidPoint();
      const CLHEP::Hep3Vector& wiredir = straw.getDirection();
    // get position from time division
      double tddist = tcal->TimeDiffToDistance(straw.index(),sh.dt());
      CLHEP::Hep3Vector pos = mid + tddist*wiredir;
      _shpos.push_back(pos);
    }
  }

  void
  TrkPatRec::findProximity(unsigned istr, double maxdist, unsigned& nprox, double& dmin) {
    unsigned nstrs = _strawhits->size();
    dmin=1e6;
    nprox = 0;
    for(unsigned jstr=0;jstr<nstrs;++jstr){
      if(jstr != istr){
        double dist = (_shpos[istr]-_shpos[jstr]).mag();
        dmin = std::min(dist,dmin);
        if(dist < maxdist)nprox++;
      }
    }
  }

  void
  TrkPatRec::preselectHits(){
    _tflags.clear();
    _tflags.reserve(_strawhits->size());
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      TrkHitFlag flag;
      if(veryloosehit(sh.energyDep(),_shpos[istr].rho()))flag.setVeryLoose(); 
      if(loosehit(sh.energyDep(),_shpos[istr].rho()))flag.setLoose(); 
      if(tighthit(sh.energyDep(),_shpos[istr].rho()))flag.setTight();
      _tflags.push_back(flag);
    }
  }

  void
  TrkPatRec::filterDeltas(){
// make a 2d plot of phi vs time to isolate the delta rays
    TSpectrum2 tspec2(_maxndelta);
    TH2F tpsp("tpsp","phi time spectrum",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
    double dbf = (_fbf-1.0)*M_PI;
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_tflags[istr].veryLoose()){
	double time = correctedTime(istr);
	double phi = _shpos[istr].phi();
// include buffer around phi to account for wrapping
	tpsp.Fill(time,phi);
	if(M_PI-phi<dbf)tpsp.Fill(time,phi-2*M_PI);
	if(phi+M_PI<dbf)tpsp.Fill(time,phi+2*M_PI);
      }
    }
// search for peaks.  Convert to an absolute threshold instead of a relative threshold
    Double_t bmax = tpsp.GetMaximum();
    double thresh(0.1);
    if(bmax > _2dthresh)thresh = _2dthresh/bmax;
    tspec2.Search(&tpsp,_2dsigma,"nobackgroundnomarkovgoff",thresh);
// fill hit information about peaks
    std::vector<DeltaInfo> dinfo;
    fillDeltaInfo(tspec2,dinfo);
// Loop over peaks
    for(size_t ip=0;ip<dinfo.size();++ip){
      DeltaInfo& delta = dinfo[ip];
// if the peak is a delta, flag its hits
      if(delta._isdelta){
	for(size_t ih =0;ih <  delta._hindex.size();++ih){
	  if(delta._hflag[ih]){
	    _tflags[delta._hindex[ih]].setDelta();
	  }
	}
      }
    }
// diagnostics
    if(_diag > 1)
      fillDeltaDiag(dinfo);
  }
  
  void 
  TrkPatRec::findTimePeaks() {
    _tpeaks.clear();
    TSpectrum tspec(_maxnpeak);
    TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
// loop over straws hits and fill time spectrum plot for tight hits
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_tflags[istr].tight()&&!_tflags[istr].delta()){
	double time = correctedTime(istr);
        timespec.Fill(time);
      }
    }
    Double_t mb = timespec.GetMaximum();
    double thresh(0.1);
    if(mb > _1dthresh) thresh = _1dthresh/mb;
    unsigned np = tspec.Search(&timespec,1,"nobackgroundnomarkovgoff",thresh);
    Float_t *xpeaks = tspec.GetPositionX();
    Float_t *ypeaks = tspec.GetPositionY();
// Loop over peaks, looking only at those with a minimum peak value
    for (unsigned ip=0; ip<np; ++ip) {
      Float_t xp = xpeaks[ip];
      Float_t yp = ypeaks[ip];
      TrkTimePeak tpeak(xp,yp);
      if(yp > _ymin){
// record hits in time with each peak, and accept them if they have a minimum # of hits
        for(unsigned istr=0; istr<nstrs;++istr){
	  if(_tflags[istr].tight()&&!_tflags[istr].delta()){
	    if(fabs(correctedTime(istr)-xp) < _maxdt)tpeak._trkptrs.push_back(istr);
	  }
	}
	if(tpeak._trkptrs.size() > _minnhits)_tpeaks.push_back(tpeak);
      }
    }
// sort the peaks so that the largest comes first
    std::sort(_tpeaks.begin(),_tpeaks.end(),greater<TrkTimePeak>());
  }
  
  double 
  TrkPatRec::correctedTime(size_t istr) const {
// correct the strawhit time for the z-dependence of the propagation time.  This depends on the pitch and
// particle assumption, FIXME!!!
    static const double slope(0.0048);
//    static const double slope(0.0);
    StrawHit const& sh = _strawhits->at(istr);
    double tcorr = sh.time()-slope*_shpos[istr].z();
    return tcorr;
  }

  void
  TrkPatRec::filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,std::vector<TrkHitFilter>& thfvec){
//  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
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

  void
  TrkPatRec::findMissingHits(KalFitResult& kalfit,std::vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    kalfit._krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_tflags[istr].veryLoose()){
	StrawHit const& sh = _strawhits->at(istr);
	if(fabs(correctedTime(istr)-kalfit._krep->t0()._t0) < _maxdtmiss) {
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


  void
  TrkPatRec::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
// straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("shpos",&_shp,"x/F:y/F:z/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("corrtime",&_corrtime,"corrtime/F");
    _shdiag->Branch("device",&_device,"device/I");
    _shdiag->Branch("sector",&_sector,"sector/I");
    _shdiag->Branch("layer",&_layer,"layer/I");
    _shdiag->Branch("straw",&_straw,"straw/I");
    _shdiag->Branch("ishpeak",&_ishpeak,"ishpeak/I");
    _shdiag->Branch("ntpeak",&_ntpeak,"ntpeak/I");
    _shdiag->Branch("tpeak",&_shtpeak,"tpeak/F");
    _shdiag->Branch("nshtpeak",&_nshtpeak,"nshtpeak/I");
    _shdiag->Branch("dmin",&_dmin,"dmin/F");
    _shdiag->Branch("n50",&_n50,"n50/I");
    _shdiag->Branch("n100",&_n100,"n100/I");
    _shdiag->Branch("n150",&_n150,"n150/I");
    _shdiag->Branch("n200",&_n200,"n200/I");
    _shdiag->Branch("mcshpos",&_mcshp,"x/F:y/F:z/F");
    _shdiag->Branch("mcedep",&_mcedep,"mcedep/F");
    _shdiag->Branch("mcemax",&_mcemax,"mcemax/F");
    _shdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/I");
    _shdiag->Branch("mcnunique",&_mcnunique,"mcnunique/I");
    _shdiag->Branch("mcnmax",&_mcnmax,"mcnmax/I");
    _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
    _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
    _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
    _shdiag->Branch("mctime",&_mctime,"mctime/F");
    _shdiag->Branch("vloose",&_vloose,"vloose/I");
    _shdiag->Branch("loose",&_loose,"loose/I");
    _shdiag->Branch("tight",&_tight,"tight/I");
    _shdiag->Branch("delta",&_delta,"delta/I");
    _shdiag->Branch("pdist",&_pdist,"pdist/F");
    _shdiag->Branch("pperp",&_pperp,"pperp/F");
    _shdiag->Branch("pmom",&_pmom,"pmom/F");
    _shdiag->Branch("mct0",&_shmct0,"mct0/F");
    _shdiag->Branch("mcmom",&_shmcmom,"mcmom/F");
    _shdiag->Branch("mctd",&_shmctd,"mctd/F");
// extend the KalFitMC track diagnostic tuple
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
// delta diagnostics
    _ddiag=tfs->make<TTree>("ddiag","delta diagnostics");
    _ddiag->Branch("iev",&_iev,"iev/I");
    _ddiag->Branch("ip",&_ip,"ip/I");
    _ddiag->Branch("isdelta",&_isdelta,"isdelta/B");
    _ddiag->Branch("oldisdelta",&_oldisdelta,"oldisdelta/B");
    _ddiag->Branch("nsh",&_nsh,"nsh/I");
    _ddiag->Branch("ndpeak",&_ndpeak,"ndpeak/I");
    _ddiag->Branch("ndmax",&_ndmax,"ndmax/I");
    _ddiag->Branch("nconv",&_nconv,"nconv/I");
    _ddiag->Branch("ndelta",&_ndelta,"ndelta/I");
    _ddiag->Branch("ncompt",&_ncompt,"ncompt/I");
    _ddiag->Branch("ngconv",&_ngconv,"ngconv/I");
    _ddiag->Branch("nebkg",&_nebkg,"nebkg/I");
    _ddiag->Branch("nprot",&_nprot,"nprot/I");
    _ddiag->Branch("ns",&_ns,"ns/I");
    _ddiag->Branch("smin",&_smin,"smin/I");
    _ddiag->Branch("smax",&_smax,"smax/I");
    _ddiag->Branch("nsmiss",&_nsmiss,"nsmiss/I");
    _ddiag->Branch("pphi",&_pphi,"pphi/F");
    _ddiag->Branch("pt",&_pt,"pt/F");
    _ddiag->Branch("prho",&_prho,"prho/F");
    _ddiag->Branch("zmin",&_zmin,"zmin/F");
    _ddiag->Branch("zmax",&_zmax,"zmax/F");
    _ddiag->Branch("zgap",&_zgap,"zgap/F");
    _ddiag->Branch("mindist",&_mindist,"mindist/F");
    _ddiag->Branch("mindt",&_mindt,"mindt/F");
    _ddiag->Branch("mindphi",&_mindphi,"mindphi/F");
    _ddiag->Branch("tmed",&_tmed,"tmed/F");
    _ddiag->Branch("pmed",&_pmed,"pmed/F");
    _ddiag->Branch("rmed",&_rmed,"rmed/F");
    _ddiag->Branch("tmean",&_tmean,"tmean/F");
    _ddiag->Branch("pmean",&_pmean,"pmean/F");
    _ddiag->Branch("rmean",&_rmean,"rmean/F");
    _ddiag->Branch("stime",&_stime,"stime/F");
    _ddiag->Branch("sphi",&_sphi,"sphi/F");
    _ddiag->Branch("srho",&_srho,"srho/F");
    _ddiag->Branch("phits",&_phits);
    _ddiag->Branch("mct0",&_dmct0,"mct0/F");
    _ddiag->Branch("mcmom",&_dmcmom,"mcmom/F");
    _ddiag->Branch("mctd",&_dmctd,"mctd/F");
    _ddiag->Branch("pmvaout",&_pmvaout,"pmvaout/F");
  }

  void
  TrkPatRec::fillStrawDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = getTrackerOrThrow();
    _nchit = 0;
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      _device = straw.id().getDevice();
      _sector = straw.id().getSector();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();

      _shp = _shpos[istr];
      _edep = sh.energyDep();
      _time = sh.time();
      _corrtime = correctedTime(istr);
     // find proximity for different radii
      double dmin(0.0);
//      findProximity(_shpos,istr,50.0,_n50,dmin);
//      findProximity(_shpos,istr,100.0,_n100,dmin);
//      findProximity(_shpos,istr,150.0,_n150,dmin);
//      findProximity(_shpos,istr,200.0,_n200,dmin);
      _dmin = dmin;
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
	bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
	if(conversion){
	  ++_nchit;
	}
      }
      _tight = _tflags[istr].tight();
      _delta = _tflags[istr].delta();
      _vloose = _tflags[istr].veryLoose();
      _loose = _tflags[istr].loose();
      _shmct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
      _shmcmom = _kfitmc.MCMom(KalFitMC::trackerMid);
      _shmctd = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
      // compare to different time peaks
      _ntpeak = _tpeaks.size();
      _nshtpeak = 0;
      _shtpeak = -1.0;
      _ishpeak = -1;
      hitIndex myindex(istr);
      if(_shmcmom >0){
	for(unsigned ipeak=0;ipeak<_tpeaks.size();++ipeak){
	  std::vector<hitIndex>::iterator ifind =
	    std::find(_tpeaks[ipeak]._trkptrs.begin(),_tpeaks[ipeak]._trkptrs.end(),myindex);
	  if(ifind != _tpeaks[ipeak]._trkptrs.end()){
	    _ishpeak = ipeak;
	    break;
	  }
	}
      }
      if(_ishpeak>=0){
	_nshtpeak = _tpeaks[_ishpeak]._trkptrs.size();
	_shtpeak = _tpeaks[_ishpeak]._tpeak;
      }
      _shdiag->Fill();
    }
  }

  void
  TrkPatRec::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp;
    TH2F *cptsp, *rptsp, *tptsp, *lptsp, *tdptsp;
    TH2F *crtsp, *rrtsp, *trtsp, *lrtsp;

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
    
    snprintf(rsname,100,"rawptspectrum%i",_iev);
    snprintf(csname,100,"convptspectrum%i",_iev);
    snprintf(tsname,100,"tightptspectrum%i",_iev);
    snprintf(lsname,100,"looseptspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltaptspectrum%i",_iev);
// buffer the range so that we don't loose any peaks
    tptsp = tfs->make<TH2F>(tsname,"time spectrum;nsec;#phi",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
    lptsp = tfs->make<TH2F>(lsname,"time spectrum;nsec;#phi",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
    rptsp = tfs->make<TH2F>(rsname,"time spectrum;nsec;#phi",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
    cptsp = tfs->make<TH2F>(csname,"time spectrum;nsec;#phi",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
    tdptsp = tfs->make<TH2F>(tdsname,"time spectrum;nsec;#phi",_ntbins,_tmin,_tmax,_npbins,-_fbf*M_PI,_fbf*M_PI);
 
    snprintf(rsname,100,"rawrtspectrum%i",_iev);
    snprintf(csname,100,"convrtspectrum%i",_iev);
    snprintf(tsname,100,"tightrtspectrum%i",_iev);
    snprintf(lsname,100,"loosertspectrum%i",_iev);
    trtsp = tfs->make<TH2F>(tsname,"time spectrum;nsec;r(cm)",_nbins,_tmin,_tmax,50,350,750);
    lrtsp = tfs->make<TH2F>(lsname,"time spectrum;nsec;r(cm)",_nbins,_tmin,_tmax,50,350,750);
    rrtsp = tfs->make<TH2F>(rsname,"time spectrum;nsec;r(cm)",_nbins,_tmin,_tmax,50,350,750);
    crtsp = tfs->make<TH2F>(csname,"time spectrum;nsec;r(cm)",_nbins,_tmin,_tmax,50,350,750);

    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = correctedTime(istr);
      double rad = _shpos[istr].perp();
      double phi = _shpos[istr].phi();
      bool conversion(false);
      // summarize the MC truth for this strawhit
      if(_kfitmc.mcData()._mcsteps != 0) {
	const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
	conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
      }
      // fill plots
      rtsp->Fill(time);
      rptsp->Fill(time,phi);
      rrtsp->Fill(time,rad);
      double dbf = (_fbf-1.0)*M_PI;
      if(M_PI-phi<dbf)rptsp->Fill(time,phi-2*M_PI);
      if(phi+M_PI<dbf)rptsp->Fill(time,phi+2*M_PI);
      if(_tflags[istr].tight()){
	ttsp->Fill(time);
	tptsp->Fill(time,phi);
	trtsp->Fill(time,rad);
	if(M_PI-phi<dbf)tptsp->Fill(time,phi-2*M_PI);
	if(phi+M_PI<dbf)tptsp->Fill(time,phi+2*M_PI);
      }
      if(_tflags[istr].tight()&&!_tflags[istr].delta()){
	tdtsp->Fill(time);
	tdptsp->Fill(time,phi);
	if(M_PI-phi<dbf)tdptsp->Fill(time,phi-2*M_PI);
	if(phi+M_PI<dbf)tdptsp->Fill(time,phi+2*M_PI);
      }
      if(_tflags[istr].veryLoose()){
	ltsp->Fill(time);
	lptsp->Fill(time,phi);
	lrtsp->Fill(time,rad);
      	if(M_PI-phi<dbf)lptsp->Fill(time,phi-2*M_PI);
	if(phi+M_PI<dbf)lptsp->Fill(time,phi+2*M_PI);
      }
      if(conversion){
	ctsp->Fill(time);
	cptsp->Fill(time,phi);
	crtsp->Fill(time,rad);
       	if(M_PI-phi<dbf)cptsp->Fill(time,phi-2*M_PI);
	if(phi+M_PI<dbf)cptsp->Fill(time,phi+2*M_PI);
      }
    }
    // find peaks, so they show up on diagnostic plot too
    TSpectrum tspec(_maxnpeak);
    Double_t mb = tdtsp->GetMaximum();
    double thresh(0.1);
    if(mb > _1dthresh) thresh = _1dthresh/mb;
    tspec.Search(tdtsp,1,"nobackgroundnomarkov",thresh);
  }

  void
  TrkPatRec::fillFitDiag(int ipeak,HelixFitResult const& helixfit,
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
      for(std::vector<hitIndex>::const_iterator istr= tpeak._trkptrs.begin(); istr != tpeak._trkptrs.end(); ++istr){
	// summarize the MC truth for this strawhit
	if(_kfitmc.mcData()._mcsteps != 0) {
	  const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr->_index); 
	  if(mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2)
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
    if(helixfit._fit.success()){
      HepVector hpar;
      HepVector hparerr;
      _hfit.helixParams(helixfit,hpar,hparerr);
      _hpar = helixpar(hpar);
      _hparerr = helixpar(hparerr);
      _hcx = helixfit._center.x(); _hcy = helixfit._center.y(); _hr = helixfit._radius;
      _hdfdz = helixfit._dfdz; _hfz0 = helixfit._fz0;
    }
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
    }
// use MC truth to define hits and seed helix
    TrkDef mctrk(_strawhits,_tpart,_fdir);
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
  
  void
  TrkPatRec::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {
    const Tracker& tracker = getTrackerOrThrow();
    const StrawHit& sh = _strawhits->at(ish);
    shinfo._pos = _shpos[ish];
    shinfo._time = sh.time();
    shinfo._corrtime = correctedTime(ish);
    shinfo._edep = sh.energyDep();
    const Straw& straw = tracker.getStraw( sh.strawIndex() );
    shinfo._device = straw.id().getDevice();
    shinfo._sector = straw.id().getSector();
    shinfo._layer = straw.id().getLayer();
    shinfo._straw = straw.id().getStraw();

    if(_kfitmc.mcData()._mcsteps != 0) {
      const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      shinfo._mcpdg = mcsum[0]._pdgid;
      shinfo._mcgen = mcsum[0]._gid;
      shinfo._mcproc = mcsum[0]._pid;
      shinfo._mcpos = mcsum[0]._pos;
      shinfo._mctime = mcsum[0]._time;
      shinfo._mcedep = mcsum[0]._esum;
      shinfo._tight = _tflags[ish].tight();
      shinfo._delta = _tflags[ish].delta();
      shinfo._vloose = _tflags[ish].veryLoose();
      shinfo._loose = _tflags[ish].loose();
      shinfo._mct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
      shinfo._mcmom = _kfitmc.MCMom(KalFitMC::trackerMid);
      shinfo._mctd = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
    }
  }
  
  void
  TrkPatRec::initializeReaders() {
// define the hit classifier
  _dhReader = new TMVA::Reader();
  _dhReader->AddVariable("ddphi",&_dhmva._ddphi);
  _dhReader->AddVariable("ddrho",&_dhmva._ddrho);
  _dhReader->AddVariable("ddt",&_dhmva._ddt);
  _dhReader->BookMVA(_dhittype,_dhitweights);

// define the peak classifier
  _dpReader = new TMVA::Reader();
  _dpReader->AddVariable("rfold",&_dpmva._rfold);
  _dpReader->AddVariable("zmin",&_dpmva._zmin);
  _dpReader->AddVariable("zmax",&_dpmva._zmax);
  _dpReader->AddVariable("zgap",&_dpmva._zgap);
  _dpReader->AddVariable("nsmiss",&_dpmva._nsmiss);
  _dpReader->AddVariable("nsh",&_dpmva._nsh);
  _dpReader->AddVariable("sphi",&_dpmva._sphi);
// load the constants
  _dpReader->BookMVA(_dpeaktype,_dpeakweights);

  }

  void
  TrkPatRec::fillDeltaInfo(TSpectrum2 const& tspec2,std::vector<DeltaInfo>& dinfo) {
    unsigned np = tspec2.GetNPeaks();
    Float_t* tpeaks(tspec2.GetPositionX());
    Float_t* ppeaks(tspec2.GetPositionY());
    dinfo.reserve(np);
// loop over peaks and fill delta information
    for (unsigned ip=0; ip<np; ++ip) {
      DeltaInfo delta;
      delta._tpeak = tpeaks[ip]; 
      delta._ppeak = ppeaks[ip];
      dinfo.push_back(delta);
    }
    // 2nd loop over pairs of peaks
    for (unsigned ip=0; ip<np; ++ip) {
    // look at the distance to the nearest peak.  Use a metric based on the plot range
      double mindist2(FLT_MAX);
      dinfo[ip]._mindist = -1.0;
      dinfo[ip]._mindt = -1.0;
      dinfo[ip]._mindphi = -1.0;
      for (unsigned jp=0; jp<np; jp++) {
	if(jp != ip){
	  Float_t tp = tpeaks[jp];
	  Float_t pp = ppeaks[jp];
	  double dt = dinfo[ip]._tpeak - tp;
	  double dphi = dinfo[ip]._ppeak - pp;
	  if(dphi > M_PI){
	    dphi -= 2*M_PI;
	  } else if(dphi < -M_PI){
	    dphi += 2*M_PI;
	  }
    // normalize distance to the plot dimensions
	  double dist2 = pow(dt/(_tmax-_tmin),2) + pow(dphi/(2*_fbf*M_PI),2);
	  if(dist2 < mindist2){
	    mindist2 = dist2;
	    dinfo[ip]._mindist = sqrt(mindist2);
	    dinfo[ip]._mindt = dt;
	    dinfo[ip]._mindphi = dphi;
	  }
	}
      }
    }
    // Loop over all the hits and assign them to peaks
    unsigned nstrs = _strawhits->size();
    for(size_t istr=0; istr<nstrs;++istr){
      // only look at selected hits
      if(_tflags[istr].veryLoose()){
	double phi = _shpos[istr].phi();
	double time = correctedTime(istr);
// loop over all the peaks
	for (unsigned ip=0; ip<np; ++ip) {
	  DeltaInfo& delta = dinfo[ip];
// select on time first
	  if(fabs(time-delta._tpeak) < _max2ddt){
// compute phi, being careful of wrapping
	    double dphi = phi - delta._ppeak;
	    if(dphi > M_PI){
	      dphi -= 2*M_PI;
	    } else if(dphi < -M_PI){
	      dphi += 2*M_PI;
	    }
// select on phi
	    if(fabs(dphi) < _maxdp){
	      delta._hindex.push_back(istr);
	      delta._htime.push_back(time);
	      delta._hphi.push_back(dphi+delta._ppeak);
	      // accumulate other hit information
	      delta._hrho.push_back(_shpos[istr].perp());
	      delta._hz.push_back(_shpos[istr].z());
	    }
	  }
	}
      }
    }
// compute summary information for each peak
    for (unsigned ip=0; ip<np; ip++) {
      DeltaInfo& delta = dinfo[ip];
      fillDeltaSummary(delta);
    }
  }

  void
  TrkPatRec::fillDeltaSummary(DeltaInfo& delta) {
// tracker, to get StrawID later
    const TTracker& tracker = dynamic_cast<const TTracker&>(getTrackerOrThrow());
    unsigned ndevices = tracker.nDevices();
// median information from all info
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc, pacc, racc;
    tacc = std::for_each( delta._htime.begin(), delta._htime.end(), tacc );
    pacc = std::for_each( delta._hphi.begin(), delta._hphi.end(), pacc );
    racc = std::for_each( delta._hrho.begin(), delta._hrho.end(), racc );
    delta._tmed = extract_result<tag::median>(tacc);
    delta._pmed = extract_result<tag::median>(pacc);
    delta._rmed = extract_result<tag::median>(racc);
    // fill the BDT information about each hit, compared to the medians
    // also accumulate information about selected hits
    accumulator_set<double, stats<tag::variance(lazy)> > tacc2,pacc2,racc2;
    std::vector<bool> devices(ndevices,false);
    std::vector<double> hz;
    hz.reserve(delta._hindex.size());
    for(size_t ih =0;ih < delta._hindex.size();++ih){
      _dhmva._ddphi = delta._hphi[ih]-delta._pmed;
      _dhmva._ddrho = delta._hrho[ih]-delta._rmed;
      _dhmva._ddt = delta._htime[ih]-delta._tmed;
      double gd = _dhReader->EvaluateMVA(_dhittype);
      delta._hgd.push_back(gd);
      if(gd > _dhitmvacut){
	delta._hflag.push_back(1);
	tacc2(delta._htime[ih]);
	pacc2(delta._hphi[ih]);
	racc2(delta._hrho[ih]);
	hz.push_back(delta._hz[ih]);
	const StrawHit& sh = _strawhits->at(delta._hindex[ih]);
	unsigned idevice = (unsigned)(tracker.getStraw(sh.strawIndex()).id().getDeviceId());
	devices[idevice] = true;
      } else {
	delta._hflag.push_back(0);
      }
    }
// compute the mean and RMS of the selected hits
    delta._ngoodhits = extract_result<tag::count>(tacc2);
    if(delta._ngoodhits>0){
      delta._tmean = extract_result<tag::mean>(tacc2);
      delta._pmean = extract_result<tag::mean>(pacc2);
      delta._rmean = extract_result<tag::mean>(racc2);
      delta._trms = sqrt(extract_result<tag::variance>(tacc2));
      delta._prms = sqrt(extract_result<tag::variance>(pacc2));
      delta._rrms = sqrt(extract_result<tag::variance>(racc2));
      std::sort(hz.begin(),hz.end());
      delta._zmin = hz.front();
      delta._zmax = hz.back();
      // look for gaps in z
      delta._zgap = 0.0;
      for(unsigned iz=1;iz<hz.size();++iz){
	if(hz[iz]-hz[iz-1] > delta._zgap)delta._zgap = hz[iz]-hz[iz-1]; 
      }
      // count 'missing' devices between first and last
      delta._ismin = 0;
      delta._ismax = ndevices-1;
      delta._nsmiss = 0;
      while(!devices[delta._ismin])++delta._ismin;
      while(!devices[delta._ismax])--delta._ismax;
      delta._ns = delta._ismax-delta._ismin+1;
      for(unsigned is =delta._ismin;is<delta._ismax;++is){
	if(!devices[is])++delta._nsmiss;
      }
      // compute final delta selection based on the delta properties
      delta._oldisdelta = false;
      delta._isdelta = false;
      if(delta._ngoodhits >= _mindp){
	// decide if these hits are deltas: if so, flag them.  This algorithm should be a neural net, FIXME!!!!
	delta._oldisdelta = delta._zgap < _maxzgap || delta._nsmiss <= _maxnsmiss || delta._rmean > _maxrho || delta._rmean < _minrho;
	
	_dpmva._rfold = fabs(delta._rmean-535);
	_dpmva._zmin = delta._zmin;
	_dpmva._zmax = delta._zmax;
	_dpmva._zgap = delta._zgap;
	_dpmva._nsmiss = delta._nsmiss;
	_dpmva._nsh = delta._ns;
	_dpmva._sphi = delta._prms;
	
	delta._pmvaout = _dpReader->EvaluateMVA(_dpeaktype);
	delta._isdelta = delta._pmvaout > _dpeakmvacut;
      }
    } else {
      delta._isdelta = false;
    }
  }

  void
  TrkPatRec::fillDeltaDiag(std::vector<DeltaInfo> const& deltas) {
    for(size_t ip=0;ip<deltas.size();++ip){
      DeltaInfo const& delta(deltas[ip]);
      _ip = ip;
      _isdelta = delta._isdelta;
      _oldisdelta = delta._oldisdelta;
      _nsh = delta._ngoodhits;
      _ndpeak = delta._hindex.size();
      _mindt = delta._mindt;
      _mindphi = delta._mindphi;
      _mindist = delta._mindist;
      _pt = delta._tpeak;
      _pphi = delta._ppeak;
      _prho =  delta._rmean;
      _tmed = delta._tmed;
      _pmed =  delta._pmed;
      _rmed =  delta._rmed;
      _tmean = delta._tmean;
      _pmean =  delta._pmean;
      _rmean =  delta._rmean;
      _stime = delta._trms;;
      _sphi = delta._prms;
      _srho = delta._rrms;
      _zmin = delta._zmin;
      _zmax = delta._zmax;
      _zgap = delta._zgap;
      _smin = delta._ismin;
      _smax = delta._ismax;
      _ns = delta._ns;
      _nsmiss = delta._nsmiss;
      _pmvaout = delta._pmvaout;
  // find MC truth information
      _nconv = 0;
      _nprot = 0;
      _ndelta= 0;
      _ncompt = 0;
      _ngconv = 0;
      _nebkg = 0;
      _phits.clear();
      for(size_t ih=0;ih< delta._hindex.size(); ++ih){
	DeltaHitInfo shinfo;
	size_t ish = delta._hindex[ih];
	fillStrawHitInfo(ish,shinfo);
	// add the BDT hit classifier
	shinfo._hflag = delta._hflag[ih];
	shinfo._hgd = delta._hgd[ih];
	shinfo._dphi = delta._hphi[ih]-delta._pmean;
	shinfo._dt = delta._htime[ih]-delta._tmean;
	shinfo._drho = delta._hrho[ih]-delta._rmean;
	_phits.push_back(shinfo);
	if(shinfo._mcgen == 2)++_nconv;
	if(abs(shinfo._mcpdg) == PDGCode::e_minus && shinfo._mcgen <0){
	  _nebkg++;
	  if(shinfo._mcproc == ProcessCode::eIoni ||shinfo._mcproc == ProcessCode::hIoni ){
	    ++_ndelta;
	} else if(shinfo._mcproc == ProcessCode::compt){
	  ++_ncompt;
	} else if(shinfo._mcproc == ProcessCode::conv){
	  ++_ngconv;
	}
      }
      if(shinfo._mcpdg == 2212)++_nprot;
      }
      if(_nconv >= 1){
	_dmct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
	_dmcmom = _kfitmc.MCMom(KalFitMC::trackerMid);
	_dmctd = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
      } else {
	_dmct0 = -1;
	_dmcmom = 0.0;
	_dmctd = -100.0;
      }
      _ddiag->Fill();
    }
  }
}
using mu2e::TrkPatRec;
DEFINE_ART_MODULE(TrkPatRec);
