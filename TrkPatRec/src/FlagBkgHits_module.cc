// $Id: FlagBkgHits_module.cc,v 1.14 2013/04/15 22:30:10 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/04/15 22:30:10 $
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
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
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// Mu2e
#include "KalmanTests/inc/KalFitMC.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "TrkPatRec/inc/ClusterStrawHits.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
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
#include <vector>
#include <set>
#include <map>
using namespace std; 
using namespace boost::accumulators;

namespace mu2e 
{
  // structs for MVAs
  // select hits close to each other and/or to the cluster
  struct DeltaHitMVA {
    Float_t _dphi; // phi diff to peak
    Float_t _drho; // rho diff to peak
    Float_t _dt;  // time diff to peak
    Float_t _stereo;  // stereo flag
    Float_t _rho;  // hit transverse radius
  };

  struct DeltaPeakMVA {
    Float_t _prho;
    Float_t _srho;
    Float_t _zmin, _zmax, _zgap;
    Float_t _nsmiss;
    Float_t _sphi;
    Float_t _ncore, _nhalo;
  };

  struct DeltaHitInfo {
    DeltaHitInfo(size_t hindex,double htime, double hphi, double hrho, double hz,double hgd,int hflag) :
      _hindex(hindex), _htime(htime), _hphi(hphi), _hrho(hrho), _hz(hz), _hgd(hgd),_hflag(hflag) {}
    size_t _hindex;
    Float_t _htime, _hphi, _hrho, _hz, _hgd, _cdist;
    Int_t _hflag;
  };

  struct DeltaInfo {
    DeltaInfo() : _tpeak(0.0), _ppeak(0.0), _rpeak(0.0),
      _tmed(0.0), _pmed(0.0),_rmed(0.0),
      _tmean(0.0), _pmean(0.0), _rmean(0.0),
      _trms(0.0), _prms(0.0),_rrms(0.0),
      _nbox(0),_nwide(0),_ncore(0),
      _zmin(0.0),_zmax(0.0),_zgap(0.0),
      _ismin(-1),_ismax(-1),_ns(-1),_nsmiss(-1),
      _idmin(-1),_idmax(-1),_nd(-1),_ndmiss(-1),
      _isdelta(false),_pmvaout(-100.0){}
        // peak information
    double _tpeak, _ppeak, _rpeak;
    // median information
    boost::accumulators::accumulator_set<double, stats<tag::median(with_p_square_quantile) > > _tacc, _pacc, _racc;
    double _tmed, _pmed, _rmed;
    double _tmean, _pmean, _rmean;
    double _trms, _prms, _rrms;
    //summary information
    unsigned _nbox, _nwide, _ncore;
    double _zmin, _zmax, _zgap;
    unsigned _ismin, _ismax, _ns, _nsmiss;
    unsigned _idmin, _idmax, _nd, _ndmiss;
    bool _isdelta;
    double _pmvaout;
    // information about the hits in the delta
    std::vector<DeltaHitInfo> _dhinfo;
  };

  class FlagBkgHits : public art::EDProducer
  {
    public:
      explicit FlagBkgHits(fhicl::ParameterSet const&);
      virtual ~FlagBkgHits();
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
      void endJob();
    private:
      // configuration parameters
      int _diag,_debug;
      int _printfreq;
      // event object labels
      std::string _shlabel, _shplabel, _shflabel;
      // straw hit selection masks
      StrawHitFlag _stmask, _deltamask, _ismask;
      // delta-ray removal parameters
      double _bz; // average Z beta for deltas
      std::string _dhittype, _dhitweights;
      std::string _dpeaktype, _dpeakweights;
      double _gdcut, _gdcore;
      unsigned _mindp;
      double _dpeakmvacut;
      int _minflag;
// input collections
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StrawHitFlagCollection* _shfcol;
// output collections
      StrawHitFlagCollection* _bkgfcol;
// internal helper functions
      void filterDeltas();
      void fillDeltaInfo(std::list<StrawHitCluster> const& clusters,std::vector<DeltaInfo>& dinfo);
      void fillDeltaSummary(DeltaInfo& delta);
      void fillDeltaDiag(std::vector<DeltaInfo> const& deltas);
      void findPrimary(DeltaInfo const& delta,art::Ptr<SimParticle>& pptr) const;
// correct the time for propagation through the tracker, given the Z component of beta
      double deltaPhi(double phi1,double phi2) const;
      void initializeReaders();
      bool findData(const art::Event& evt);
      void createDiagnostics();
      void fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const;
      // MC tools
      KalFitMC _kfitmc;
      // clusterer
      ClusterStrawHits _clusterer;
      // delta removal diagnostics
      TTree* _ddiag;
      TMVA::Reader *_dhReader; // assign hits to delta peaks
      DeltaHitMVA _dhmva; // input variables to TMVA for delta hit selection
      DeltaPeakMVA _dpmva;
      TMVA::Reader* _dpReader; // classify delta peaks
      Int_t _ip, _iev;
      Bool_t _isdelta;
      Float_t _pphi, _pt, _prho;
      Float_t _zmin, _zmax, _zgap;
      Int_t _ns, _nd, _smin, _smax, _nsmiss, _ndmiss;
      Int_t _nbox,_nwide,_ncore;
      Int_t _pid;
      Int_t _nconv, _ndelta, _ncompt, _ngconv, _nebkg, _nprot, _nprimary;
      Int_t _ppdgid,_pgen, _pproc;
      std::vector<StrawHitInfo> _phits;
      Float_t _dmct0, _dmcmom, _dmctd;
      Float_t _tmed, _pmed, _rmed;
      Float_t _tmean, _pmean, _rmean;
      Float_t _stime, _sphi, _srho;
      Float_t _pmvaout;
      // histograms
      TH1F* _nhits, *_niter, *_nclusters, *_nchanged;
  };

  FlagBkgHits::FlagBkgHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _shlabel(pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
    _shplabel(pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shflabel(pset.get<std::string>("StrawHitFlagCollectionLabel","FlagStrawHits")),
    _stmask(StrawHitFlag::stereo),
    _deltamask(StrawHitFlag::delta),
    _ismask(StrawHitFlag::isolated),
    _dhittype(pset.get<std::string>("DeltaHitTMVAType","MLP method")),
    _dpeaktype(pset.get<std::string>("DeltaPeakTMVAType","MLP method")),
    _gdcut(pset.get<double>("DeltaHitMVACut",0.7)),
    _gdcore(pset.get<double>("DeltaHitMVACoreCut",0.85)),
    _mindp(pset.get<unsigned>("MinDeltaHits",6)),
    _dpeakmvacut(pset.get<double>("DeltaPeakMVACut",0.95)),
    _minflag(pset.get<int>("DeltaHitFlag",1)),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC",fhicl::ParameterSet())),
    _clusterer(pset.get<fhicl::ParameterSet>("ClusterStrawHits",fhicl::ParameterSet()))
  {
    // location-independent files
    ConfigFileLookupPolicy configFile;
    std::string dhitweights = pset.get<std::string>("DeltaHitTMVAWeights","TrkPatRec/test/StereoGDist_MLP.weights.xml");
    std::string dpeakweights = pset.get<std::string>("DeltaPeakTMVAWeights","TrkPatRec/test/StereoPeak_MLP.weights.xml");
    _dhitweights = configFile(dhitweights);
    _dpeakweights = configFile(dpeakweights);
    produces<StrawHitFlagCollection>();
// eventually, should also produce a collection of delta-rays and their properties, FIXME!!
  }

  FlagBkgHits::~FlagBkgHits(){}

  void FlagBkgHits::beginJob(){
    // create diagnostics if requested
    createDiagnostics();
    // initialize the TMVA readers
    initializeReaders();
  }

  void FlagBkgHits::beginRun(art::Run& ){}

  void FlagBkgHits::produce(art::Event& event ) {
    // create output
    unique_ptr<StrawHitFlagCollection> bkgfcol(new StrawHitFlagCollection);
    _bkgfcol = bkgfcol.get();
    // event printout
    _iev=event.id().event();
    if((_iev%_printfreq)==0)cout<<"FlagBkgHits: event="<<_iev<<endl;
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
// merge the input flags
    size_t nsh = _shcol->size();
    for(size_t ish=0;ish<nsh;++ish){
      StrawHitFlag flag(_shfcol->at(ish));;
      flag.merge(_shpcol->at(ish).flag());
      bkgfcol->push_back(flag);
    }
// find clusters in time/phi/rho space
    std::vector<DeltaInfo> dinfo;
// test of new straw hit clustering
    StrawHitClusterList clusters;
    _clusterer.findClusters(*_shcol,*_shpcol,*bkgfcol,clusters);
    _nhits->Fill(clusters._nhits);
    _niter->Fill(clusters._niter);
    _nchanged->Fill(clusters._nchanged);
    _nclusters->Fill(clusters._clist.size());
    fillDeltaInfo(clusters._clist,dinfo);
// loop over clusters and classify them
    for(size_t ip=0;ip<dinfo.size();++ip){
      DeltaInfo& delta = dinfo[ip];
  // fill summary information
      fillDeltaSummary(delta);
      // if the peak is a delta, flag its hits
      if(delta._isdelta){
	for(size_t ih =0;ih <  delta._dhinfo.size();++ih){
	  if(delta._dhinfo[ih]._hflag >= _minflag){
	    bkgfcol->at(delta._dhinfo[ih]._hindex).merge(_deltamask);
	  }
	}
      }
    }
    // flag the unclustered hits as 'isolated'
    for(size_t ish=0;ish<nsh;++ish){
      if(clusters._cids[ish]<0)bkgfcol->at(ish).merge(_ismask);
    }
    // diagnostics
    if(_diag > 0)
      fillDeltaDiag(dinfo);
    // put the background flag into the event
    event.put(std::move(bkgfcol));
  }

  void FlagBkgHits::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

  // find the input data objects 
  bool FlagBkgHits::findData(const art::Event& evt){
    _shcol = 0; _shpcol = 0;  _shfcol = 0;
    art::Handle<mu2e::StrawHitCollection> shH;
    if(evt.getByLabel(_shlabel,shH))
      _shcol = shH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shplabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shflabel,shflagH))
      _shfcol = shflagH.product();
    return _shcol != 0 && _shpcol != 0 && _shfcol != 0;
  }

  void FlagBkgHits::fillDeltaInfo(std::list<StrawHitCluster> const& clusters,std::vector<DeltaInfo>& dinfo){
    for (std::list<StrawHitCluster>::const_iterator icl=clusters.begin();icl!=clusters.end();++icl){
      StrawHitCluster const& cluster = *icl;
      DeltaInfo dp;
      dp._tpeak = dp._tmed = cluster.time();
      dp._ppeak = dp._pmed = cluster.pos().phi();
      dp._rpeak = dp._rmed = cluster.pos().perp();
      // find precise peak position
      for(size_t ich=0;ich<cluster.hits().size();++ich){
	size_t ish =cluster.hits()[ich]._index; 
	StrawHit const& sh = _shcol->at(ish);
	StrawHitPosition const& shp = _shpcol->at(ish);
	double ct = sh.time();
	double phi = dp._ppeak+deltaPhi(shp.pos().phi(),dp._ppeak);
	double rho = shp.pos().perp();
	double dphi = deltaPhi(phi,dp._ppeak);
	_dhmva._dphi = dphi;
	_dhmva._drho = rho - dp._rmed;
	_dhmva._dt = ct - dp._tmed;
	_dhmva._stereo = shp.flag().hasAllProperties(_stmask);
	_dhmva._rho = rho;
	double gd = _dhReader->EvaluateMVA(_dhittype);
	int iflag(0);
	if(gd > _gdcore)
	  iflag = 2;
	else if(gd > _gdcut)
	  iflag = 1;
	DeltaHitInfo dhinfo(ish,ct,phi,rho,shp.pos().z(),gd,iflag);
	dhinfo._cdist = cluster.hits()[ich]._dist;
	dp._dhinfo.push_back(dhinfo);
      }
      fillDeltaSummary(dp);
      dinfo.push_back(dp);
    }
  }
  
  void FlagBkgHits::fillDeltaSummary(DeltaInfo& delta) {
    // tracker, to get StrawID later
    const TTracker& tracker = dynamic_cast<const TTracker&>(getTrackerOrThrow());
    unsigned ndevices = tracker.nDevices();
    unsigned nstations = ndevices/2;
    // fill the BDT information about each hit, compared to the medians
    // also accumulate information about selected hits
    accumulator_set<double, stats<tag::mean,tag::variance(lazy)> > tacc,pacc,racc;
    std::vector<bool> devices(ndevices,false);
    std::vector<bool> stations(nstations,false);
    std::vector<double> hz;
    hz.reserve(delta._dhinfo.size());
    delta._ncore= delta._nwide=0;
    delta._nbox = delta._dhinfo.size();
    for(size_t ih =0;ih < delta._dhinfo.size();++ih){
      DeltaHitInfo& dhinfo = delta._dhinfo[ih];
      if(dhinfo._hgd > _gdcut){
	tacc(dhinfo._htime);
	pacc(dhinfo._hphi);
	racc(dhinfo._hrho);
	const StrawHit& sh = _shcol->at(dhinfo._hindex);
	unsigned idevice = (unsigned)(tracker.getStraw(sh.strawIndex()).id().getDeviceId());
	unsigned istation = idevice/2;
	devices[idevice] = true;
	stations[istation] = true;
	hz.push_back(dhinfo._hz);
	++delta._nwide;
	if(dhinfo._hgd > _gdcore)++delta._ncore;
      }
    }
    if(extract_result<tag::count>(tacc) > 0){
      delta._tmean = extract_result<tag::mean>(tacc);
      delta._pmean = extract_result<tag::mean>(pacc);
      delta._rmean = extract_result<tag::mean>(racc);
      delta._trms = extract_result<tag::variance>(tacc);
      delta._prms = extract_result<tag::variance>(pacc);
      delta._rrms = extract_result<tag::variance>(racc);
      delta._trms = (delta._trms>0.0) ? sqrt(delta._trms) : 0.0;
      delta._prms = (delta._prms>0.0) ? sqrt(delta._prms) : 0.0;
      delta._rrms = (delta._rrms>0.0) ? sqrt(delta._rrms) : 0.0;
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
      delta._ismax = ndevices/2-1;
      delta._idmin = 0;
      delta._idmax = ndevices-1;
      delta._nsmiss = 0;
      delta._ndmiss = 0;
      while(!stations[delta._ismin])++delta._ismin;
      while(!stations[delta._ismax])--delta._ismax;
      delta._ns = delta._ismax-delta._ismin+1;
      for(unsigned is =delta._ismin;is<delta._ismax;++is){
	if(!stations[is])++delta._nsmiss;
      }
      while(!devices[delta._idmin])++delta._idmin;
      while(!devices[delta._idmax])--delta._idmax;
      delta._nd = delta._idmax-delta._idmin+1;
      for(unsigned id =delta._idmin;id<delta._idmax;++id){
	if(!devices[id])++delta._ndmiss;
      }
      // compute final delta selection based on the delta properties
      delta._isdelta = false;
      _dpmva._prho = delta._rmed;
      _dpmva._srho = delta._rrms;
      _dpmva._zmin = delta._zmin;
      _dpmva._zmax = delta._zmax;
      _dpmva._zgap = delta._zgap;
      _dpmva._nsmiss = delta._nsmiss;
      _dpmva._sphi = delta._prms;
      _dpmva._ncore = delta._ncore;
      _dpmva._nhalo = delta._nwide-delta._ncore;

      delta._pmvaout = _dpReader->EvaluateMVA(_dpeaktype);
      delta._isdelta = delta._nwide >= _mindp && delta._pmvaout > _dpeakmvacut;
    } 
  }

 double FlagBkgHits::deltaPhi(double phi1,double phi2) const {
    double dphi = phi2-phi1;
    if(dphi > M_PI){
      dphi -= 2*M_PI;
    } else if(dphi < -M_PI){
      dphi += 2*M_PI;
    }
    return dphi;
  }

  void FlagBkgHits::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
    _nhits=tfs->make<TH1F>("nhits","N hits used in Clustering",100,0,5000);
    _niter=tfs->make<TH1F>("niter","N Cluster Iterations",20,-0.5,19.5);
    _nchanged=tfs->make<TH1F>("nchanged","N Hits changed in last iteration",100,-0.5,99.5);
    _nclusters=tfs->make<TH1F>("nclusters","N Clusters",100,0,200);
    if(_diag > 0){
      // detailed delta diagnostics
      _ddiag=tfs->make<TTree>("ddiag","delta diagnostics");
      _ddiag->Branch("iev",&_iev,"iev/I");
      _ddiag->Branch("ip",&_ip,"ip/I");
      _ddiag->Branch("isdelta",&_isdelta,"isdelta/B");
      _ddiag->Branch("nbox",&_nbox,"nbox/I");
      _ddiag->Branch("nwide",&_nwide,"nwide/I");
      _ddiag->Branch("ncore",&_ncore,"nwide/I");
      _ddiag->Branch("pid",&_pid,"pid/I");
      _ddiag->Branch("ppdgid",&_ppdgid,"ppdgid/I");
      _ddiag->Branch("pgen",&_pgen,"pgen/I");
      _ddiag->Branch("pproc",&_pproc,"pproc/I");
      _ddiag->Branch("nprimary",&_nprimary,"nprimary/I");
      _ddiag->Branch("nconv",&_nconv,"nconv/I");
      _ddiag->Branch("ndelta",&_ndelta,"ndelta/I");
      _ddiag->Branch("ncompt",&_ncompt,"ncompt/I");
      _ddiag->Branch("ngconv",&_ngconv,"ngconv/I");
      _ddiag->Branch("nebkg",&_nebkg,"nebkg/I");
      _ddiag->Branch("nprot",&_nprot,"nprot/I");
      _ddiag->Branch("ns",&_ns,"ns/I");
      _ddiag->Branch("nd",&_nd,"nd/I");
      _ddiag->Branch("smin",&_smin,"smin/I");
      _ddiag->Branch("smax",&_smax,"smax/I");
      _ddiag->Branch("nsmiss",&_nsmiss,"nsmiss/I");
      _ddiag->Branch("ndmiss",&_ndmiss,"ndmiss/I");
      _ddiag->Branch("pphi",&_pphi,"pphi/F");
      _ddiag->Branch("pt",&_pt,"pt/F");
      _ddiag->Branch("prho",&_prho,"prho/F");
      _ddiag->Branch("zmin",&_zmin,"zmin/F");
      _ddiag->Branch("zmax",&_zmax,"zmax/F");
      _ddiag->Branch("zgap",&_zgap,"zgap/F");
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
  }

  void FlagBkgHits::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHit const& sh = _shcol->at(ish);
    StrawHitPosition const& shp = _shpcol->at(ish);
    StrawHitFlag const& shflag = _bkgfcol->at(ish);
    shinfo._pos = shp.pos();
    shinfo._time = sh.time();
    shinfo._rho = shp.pos().perp();
    shinfo._pres = shp.posRes(StrawHitPosition::phi);
    shinfo._rres = shp.posRes(StrawHitPosition::rho);
    shinfo._chisq = shp.chisq();
    shinfo._edep = sh.energyDep();
    const Straw& straw = tracker.getStraw( sh.strawIndex() );
    shinfo._device = straw.id().getDevice();
    shinfo._sector = straw.id().getSector();
    shinfo._layer = straw.id().getLayer();
    shinfo._straw = straw.id().getStraw();
    shinfo._esel = shflag.hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel = shflag.hasAllProperties(StrawHitFlag::radsel);
    shinfo._delta = shflag.hasAllProperties(StrawHitFlag::delta);
    shinfo._stereo = shflag.hasAllProperties(StrawHitFlag::stereo);

    if(_kfitmc.mcData()._mcsteps != 0) {
      const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      shinfo._mcpdg = mcsum[0]._pdgid;
      shinfo._mcgen = mcsum[0]._gid;
      shinfo._mcproc = mcsum[0]._pid;
      shinfo._mcpos = mcsum[0]._pos;
      shinfo._mctime = mcsum[0]._time;
      shinfo._mcedep = mcsum[0]._esum;
      shinfo._mct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
      shinfo._mcmom = mcsum[0]._mom.mag();
      shinfo._mctd = mcsum[0]._mom.cosTheta();
    }
  }

  void FlagBkgHits::initializeReaders() {
    // define the hit classifier

    _dhReader = new TMVA::Reader();
    _dhReader->AddVariable("_dphi",&_dhmva._dphi);
    _dhReader->AddVariable("_drho",&_dhmva._drho);
    _dhReader->AddVariable("_dt",&_dhmva._dt);
    _dhReader->AddVariable("_stereo",&_dhmva._stereo);
    _dhReader->AddVariable("_rho",&_dhmva._rho);
    _dhReader->BookMVA(_dhittype,_dhitweights);

    // define the peak classifier
    _dpReader = new TMVA::Reader();
    _dpReader->AddVariable("prho",&_dpmva._prho);
    _dpReader->AddVariable("srho",&_dpmva._srho);
    _dpReader->AddVariable("zmin",&_dpmva._zmin);
    _dpReader->AddVariable("zmax",&_dpmva._zmax);
    _dpReader->AddVariable("zgap",&_dpmva._zgap);
    _dpReader->AddVariable("nsmiss",&_dpmva._nsmiss);
    _dpReader->AddVariable("sphi",&_dpmva._sphi);
    _dpReader->AddVariable("ncore",&_dpmva._ncore);
    _dpReader->AddVariable("nhalo",&_dpmva._nhalo);
    // load the constants
    _dpReader->BookMVA(_dpeaktype,_dpeakweights);

  }

  void FlagBkgHits::fillDeltaDiag(std::vector<DeltaInfo> const& deltas)
  {
    for(size_t ip=0;ip<deltas.size();++ip){
      DeltaInfo const& delta(deltas[ip]);
      _ip = ip;
      _isdelta = delta._isdelta;
      _nbox = delta._nbox;
      _nwide = delta._nwide;
      _ncore = delta._ncore;
      _pt = delta._tpeak;
      _pphi = delta._ppeak;
      _prho =  delta._rpeak;
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
      _nd = delta._nd;
      _nsmiss = delta._nsmiss;
      _ndmiss = delta._ndmiss;
      _pmvaout = delta._pmvaout;
      // find MC truth information for this delta
      art::Ptr<SimParticle> pptr;
      findPrimary(delta,pptr);
      if(pptr.isNonnull()){
	_pid = pptr->id().asInt();
	_ppdgid = pptr->pdgId();
	_pproc = pptr->creationCode();
	if( pptr->genParticle().isNonnull())
	  _pgen = pptr->genParticle()->generatorId().id();
	else
	  _pgen = -1;
      } else {
	_pid =-1;
	_ppdgid = -1;
	_pproc = -1;
	_pgen = -1;
      }
// loop over hits in this delta and classify them
      _nconv = 0;
      _nprot = 0;
      _ndelta= 0;
      _ncompt = 0;
      _ngconv = 0;
      _nebkg = 0;
      _nprimary = 0;
      _phits.clear();
      _phits.reserve(delta._dhinfo.size());
      for(size_t ih=0;ih< delta._dhinfo.size(); ++ih){
	DeltaHitInfo const& dhinfo = delta._dhinfo[ih];
	size_t ish = dhinfo._hindex;
	StrawHitInfo shinfo;
	fillStrawHitInfo(ish,shinfo);
	const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
	art::Ptr<SimParticle> spp = mcsum[0]._spp;
	shinfo._relation=KalFitMC::none;
	if(spp.isNonnull()){
	  shinfo._relation = KalFitMC::relationship(spp,pptr);
	}
	shinfo._hflag = dhinfo._hflag;
	shinfo._hgd = dhinfo._hgd;
	shinfo._cdist = dhinfo._cdist;
	shinfo._dphi = dhinfo._hphi-delta._ppeak;
	shinfo._dt = dhinfo._htime-delta._tmed;
	shinfo._drho = dhinfo._hrho-delta._rmed;
	if(shinfo._relation==0)++_nprimary;
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
	_phits.push_back(shinfo);
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

  void FlagBkgHits::findPrimary(DeltaInfo const& delta,art::Ptr<SimParticle>& pptr) const { 
    // find the unique simparticles which produced these hits
    std::set<art::Ptr<SimParticle> > pp;
    for(size_t ih=0;ih< delta._dhinfo.size(); ++ih){
      size_t ish = delta._dhinfo[ih]._hindex;
      const std::vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      art::Ptr<SimParticle> spp = mcsum[0]._spp;
      if(spp.isNonnull()){
	// if the particle is a conversion and the momentum is really low, don't count it as part
	// of the parentage
	if(spp->genParticle().isNonnull() &&
	    spp->genParticle()->generatorId()== GenId::conversionGun){
	  if(spp->genParticle()->momentum().vect().mag()-mcsum[0]._mom.mag() < 25.0){
	    pp.insert(spp);
	  }
	} else 
	  pp.insert(spp);
      }
    }
    // map these particles back to each other, to compress out particles generated inside the peak 
    std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> > spmap;
    // look for particles produced at the same point, like conversions.  It's not enough to look for the same parent,
    // as that parent could produce multiple daughters at different times.  Regardless of mechanism or genealogy, call these 'the same'
    // as they will contribute equally to the spiral
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> sppi = *ipp;
      spmap[sppi] = sppi;
    }
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> sppi = *ipp;
      if(sppi->genParticle().isNull()){
	std::set<art::Ptr<SimParticle> >::iterator jpp=ipp;++jpp;
	for(;jpp!=pp.end();++jpp){
	  art::Ptr<SimParticle> sppj = *jpp;
	  if(sppj->genParticle().isNull()){
	    // call the particles 'the same' if they are related and were produced near each other
	    KalFitMC::relation rel = KalFitMC::relationship(sppi,sppj);
	    if(rel==KalFitMC::daughter || rel == KalFitMC::udaughter){
	      spmap[sppi] = sppj;
	      break;
	    } else if(rel == KalFitMC::mother || rel == KalFitMC::umother){
	      spmap[sppj] = sppi;
	    } else if(rel == KalFitMC::sibling || rel == KalFitMC::usibling){
	      double dist = (sppj->startPosition() - sppi->startPosition()).mag();
	      if(dist < 10.0){
		if(sppi->id().asInt() > sppj->id().asInt())
		  spmap[sppi] = sppj;
		else
		  spmap[sppj] = sppi;
	      }
	    }
	  }
	}
      }
    }
    // check for remapping
    bool changed(true);
    while(changed){
      changed = false;
      for(std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator im = spmap.begin();im!=spmap.end();++im){
	std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator ifnd = spmap.find(im->second);
	if(ifnd->second != ifnd->first){
	  changed = true;
	  spmap[im->first] = ifnd->second;
	}
      }
    }
    // find the most likely ultimate parent for this peak.  Also fill general info
    std::map<int,int> mode;
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> spp = *ipp;
      int mcid(-1);
      // map back to the ultimate parent
      spp = spmap[spp];
      mcid = spp->id().asInt();
      std::map<int,int>::iterator ifnd = mode.find(mcid);
      if(ifnd != mode.end())
	++(ifnd->second);
      else
	mode[mcid] = 1;
    }
    int max(0);
    std::map<int,int>::iterator imax = mode.end();
    for(std::map<int,int>::iterator im=mode.begin();im!=mode.end();++im){
      if(im->second>max){
	imax=im;
	max = im->second;
      }
    }
    unsigned pid(0);
    if(imax != mode.end())
      pid=imax->first;
    for(std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator im = spmap.begin();im!=spmap.end();++im){
      if(im->first->id().asInt() == pid){
	pptr = im->first;
	break;
      }
    }
  }
}
using mu2e::FlagBkgHits;


DEFINE_ART_MODULE(FlagBkgHits);
