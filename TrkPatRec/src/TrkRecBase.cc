//
// TTracker commons for Kalman Fit
//
// $Id: TrkRecBase.cc,v 1.1 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
//
// Original author D. Brown and G. Tassielli
//

// data
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPoca.hh"
// Mu2e
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
// root
#include "TMarker.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
// C++
#include <iostream>

using namespace std;
using namespace boost::accumulators;

// Mu2e
#include "TrkPatRec/inc/TrkRecBase.hh"

namespace mu2e 
{

  TrkRecBase::TrkRecBase(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
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
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet())),
    _hfit(pset.get<fhicl::ParameterSet>("HelixFit",fhicl::ParameterSet())),
    _payloadSaver(pset),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC",fhicl::ParameterSet()))
  {
    // tag the data product instance by the direction and particle type found by this fitter
    _iname = _fdir.name() + _tpart.name();
//    produces<KalRepCollection>(_iname);
//    produces<KalRepPtrCollection>(_iname);
//    produces<KalRepPayloadCollection>();
//    produces<StrawHitFlagCollection>(_iname);
    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    // location-independent files
    ConfigFileLookupPolicy configFile;
    std::string weights = pset.get<std::string>("PeakMVAWeights","TrkPatRec/test/TimePeak.weights.xml");
    _PMVAWeights = configFile(weights);
  }

  TrkRecBase::~TrkRecBase(){}

  void TrkRecBase::bgnJob(){
    initializeReaders();
//    initializeReaders(_peakMVA, _pmva, _PMVAType, _PMVAWeights);
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

  // find the input data objects 
  bool TrkRecBase::findData(const art::Event& evt){
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

//  void TrkRecBase::filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec){
//    //  Trajectory info
//    Hep3Vector tdir;
//    HepPoint tpos;
//    traj.getInfo(0.0,tpos,tdir);
//    // tracker and conditions
//    const Tracker& tracker = getTrackerOrThrow();
//    ConditionsHandle<TrackerCalibrations> tcal("ignored");
//    const StrawHitCollection* hits = mytrk.strawHitCollection();
//    const vector<hitIndex>& indices = mytrk.strawHitIndices();
//    vector<hitIndex> goodhits;
//    for(unsigned ihit=0;ihit<indices.size();++ihit){
//      StrawHit const& sh = hits->at(indices[ihit]._index);
//      Straw const& straw = tracker.getStraw(sh.strawIndex());
//      CLHEP::Hep3Vector hpos = straw.getMidPoint();
//      CLHEP::Hep3Vector hdir = straw.getDirection();
//      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
//      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
//      TrkLineTraj htraj(spt,hdir,-20,20);
//      // estimate flightlength along track.  This assumes a constant BField!!!
//      double fltlen = (hpos.z()-tpos.z())/tdir.z();
//      TrkPoca hitpoca(traj,fltlen,htraj,0.0);
//      // flag hits with small residuals
//      if(fabs(hitpoca.doca()) < maxdoca){
//	goodhits.push_back(indices[ihit]);
//      }
//      // optional diagnostics
//      if(_diag > 0){
//	// summarize the MC truth for this strawhit
//	TrkHitFilter thfilter;
//	HepPoint tpos =  traj.position(hitpoca.flt1());
//	thfilter._pos = CLHEP::Hep3Vector(tpos.x(),tpos.y(),tpos.z());
//	thfilter._doca = hitpoca.doca();
//	if(_kfitmc.mcData()._mcsteps != 0){
//	  const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ihit);
//	  thfilter._mcpdg = mcsum[0]._pdgid;
//	  thfilter._mcgen = mcsum[0]._gid;
//	  thfilter._mcproc = mcsum[0]._pid;
//	}
//	thfvec.push_back(thfilter);
//      }
//    }
//    // update track
//    mytrk.setIndices(goodhits);
//  }

  void TrkRecBase::findMissingHits(KalFitResult& kalfit,vector<hitIndex>& misshits) {
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

  void TrkRecBase::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
    // straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("shpos",&_shp,"x/F:y/F:z/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("deltat",&_deltat,"deltat/F");
    _shdiag->Branch("rho",&_rho,"rho/F");
    _shdiag->Branch("device",&_device,"device/I");
    _shdiag->Branch("sector",&_sector,"sector/I");
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
    _shdiag->Branch("mcemax",&_mcemax,"mcemax/F");
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
    _shdiag->Branch("mct0",&_shmct0,"mct0/F");
    _shdiag->Branch("mcmom",&_shmcmom,"mcmom/F");
    _shdiag->Branch("mctd",&_shmctd,"mctd/F");
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
  }

  void TrkRecBase::fillStrawDiag() {
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
      _device = straw.id().getDevice();
      _sector = straw.id().getSector();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();

      _shp = shp.pos();
      _stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);
      _edep = sh.energyDep();
      _time = sh.time();
      _deltat = sh.dt();
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
      _mcemax = -1;
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
      if(_kfitmc.mcData()._mcsteps != 0){
	const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
	_mcnunique = mcsum.size();
	// compute energy sum
	_mcedep = 0.0;
	for(vector<MCHitSum>::const_iterator isum=mcsum.begin(); isum != mcsum.end(); ++isum){
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
	_mcoe = mcsum[0]._spp->startMomentum().e();
	_mcom = mcsum[0]._spp->startMomentum().vect().mag();
	_mcshlen = (mcsum[0]._pos-straw.getMidPoint()).dot(straw.getDirection());
	_mcshd = (mcsum[0]._pos-straw.getMidPoint()).dot(straw.getDirection().cross(mcsum[0]._mom).unit());
	bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2 && mcsum[0]._mom.mag()>90.0);
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
	art::Ptr<SimParticle> sp = mcptr[0]->simParticle();
	if(sp.isNonnull() && sp->parent().isNonnull()){
	  const art::Ptr<SimParticle>& psp = sp->parent();
	  _mcppdg = psp->pdgId();
	  _mcpproc = psp->creationCode();
	  _mcptime = psp->startGlobalTime();
	  _mcpop = det->toDetector(psp->startPosition());
	  _mcpoe = psp->startMomentum().e();
	  _mcpom = psp->startMomentum().vect().mag();
	}
// generator information
	if(sp.isNonnull()){
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
	_xtalk = mcsum[0]._mcsid != sh.strawIndex();
      }
      _esel = _flags->at(istr).hasAllProperties(StrawHitFlag::energysel);
      _rsel = _flags->at(istr).hasAllProperties(StrawHitFlag::radsel);
      _timesel = _flags->at(istr).hasAllProperties(StrawHitFlag::timesel);
      _stereo = _flags->at(istr).hasAllProperties(StrawHitFlag::stereo);
      _isolated = _flags->at(istr).hasAllProperties(StrawHitFlag::isolated);
      _delta = _flags->at(istr).hasAllProperties(StrawHitFlag::delta);
      _shpres = _shpcol->at(istr).posRes(StrawHitPosition::phi);
      _shrres = _shpcol->at(istr).posRes(StrawHitPosition::rho);
      _shmct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
      _shmcmom = _kfitmc.MCMom(KalFitMC::trackerMid);
      _shmctd = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
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

  void TrkRecBase::fillTimeDiag() {
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
      if(_kfitmc.mcData()._mcsteps != 0) {
	const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr); 
	conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2 && mcsum[0]._mom.mag()>90.0);
	proton = mcsum[0]._pdgid==2212;
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

  void TrkRecBase::fillFitDiag(int ipeak,HelixFitResult const& helixfit,
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
	if(_kfitmc.mcData()._mcsteps != 0) {
	  const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(istr->_index); 
	  if(mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2 && mcsum[0]._mom.mag()>90.0)
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
    for(vector<TrkStrawHit*>::const_iterator ish=kalfit._hits.begin();ish!=kalfit._hits.end();++ish){
      if((*ish)->usability()==3)++_nadd;
    }
    // fill kalman fit info.  This needs to be last, as it calls TTree::Fill().
    _kfitmc.kalDiag(kalfit._krep);
  }


  void TrkRecBase::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {
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
    shinfo._device = straw.id().getDevice();
    shinfo._sector = straw.id().getSector();
    shinfo._layer = straw.id().getLayer();
    shinfo._straw = straw.id().getStraw();
    shinfo._esel = shp.flag().hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel = shp.flag().hasAllProperties(StrawHitFlag::radsel);
    shinfo._delta = shp.flag().hasAllProperties(StrawHitFlag::delta);
    shinfo._stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);

    if(_kfitmc.mcData()._mcsteps != 0) {
      const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      shinfo._mcpdg = mcsum[0]._pdgid;
      shinfo._mcgen = mcsum[0]._gid;
      shinfo._mcproc = mcsum[0]._pid;
      shinfo._mcpos = mcsum[0]._pos;
      shinfo._mctime = mcsum[0]._time;
      shinfo._mcedep = mcsum[0]._esum;
      shinfo._mct0 = _kfitmc.MCT0(KalFitMC::trackerMid);
      shinfo._mcmom = _kfitmc.MCMom(KalFitMC::trackerMid);
      shinfo._mctd = _kfitmc.MCHelix(KalFitMC::trackerMid)._td;
    }
  }

  void TrkRecBase::fillPeakDiag(size_t ip,TrkTimePeak const& tp) {
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


      const vector<MCHitSum>& mcsum = _kfitmc.mcHitSummary(ish);
      if(mcsum[0]._gid == 2 && mcsum[0]._mom.mag()>90.0){
	_ncphits++;
	if(fabs(dt) > _cdtimemax)_cdtimemax = fabs(dt);
      }
      if(mcsum[0]._gid == 2 && mcsum[0]._mom.mag()>90.0 && dphi > _cdphimax)_cdphimax = dphi;

      _pmva._dt = dt;
      _pmva._dphi = dphi;
      _pmva._rho = pos.perp();
      double mvaout = _peakMVA->EvaluateMVA(_PMVAType);

      TimePeakHitInfo tph;
      tph._dt = dt;
      tph._dphi = dphi;
      tph._rho = pos.perp();
      tph._mva = mvaout;
      tph._mcpdg = mcsum[0]._pdgid;
      tph._mcpdg = mcsum[0]._pdgid;
      tph._mcpdg = mcsum[0]._pdgid;
      tph._mcgen = mcsum[0]._gid;
      tph._mcproc = mcsum[0]._pid;
      _tphinfo.push_back(tph);
    }
    _tpdiag->Fill();
    if(_ncphits > 0.5*_nchits)_icepeak = ip;
  }

  void TrkRecBase::initializeReaders() {
    _peakMVA = new TMVA::Reader();
    _peakMVA->AddVariable("_dt",&_pmva._dt);
    _peakMVA->AddVariable("_dphi",&_pmva._dphi);
    _peakMVA->AddVariable("_rho",&_pmva._rho);
    _peakMVA->BookMVA(_PMVAType,_PMVAWeights);
  }

}
