//
// Module to perform BaBar Kalman fit
//
// $Id: TrkPatRec_module.cc,v 1.2 2011/09/08 19:06:05 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/09/08 19:06:05 $
//
// framework
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh" 
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh"
#include "TrkPatRec/inc/TrkHelixFit.hh"
#include "TrkBase/TrkPoca.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TSpectrum.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>


using namespace std; 

namespace mu2e 
{
// struct to keep track of hits in a time peak
  struct TrkTimePeak {
    std::vector<size_t> _trkptrs;
    double _tpeak;
    double _peakmax;
    TrkTimePeak(double tpeak,double ymax) : _tpeak(tpeak),_peakmax(ymax) {}
    bool operator < (TrkTimePeak const& other ) const { return _trkptrs.size() < other._trkptrs.size(); }
  };

// struct for flagging hits
  struct TrkHitFlag {
#define TIGHTBIT 0x1
#define LOOSEBIT 0x2
    TrkHitFlag() : _iflag(0) {}
    void setLoose() { _iflag |= LOOSEBIT; }
    void setTight() { _iflag |= TIGHTBIT; }

    bool loose() const { return (_iflag & LOOSEBIT) != 0; }
    bool tight() const { return (_iflag & TIGHTBIT) != 0; }

    unsigned _iflag;
  };

  class TrkPatRec : public art::EDProducer
  {
  public:
    enum fittype { helix=0,seed,kal};
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
    // event object labels
    std::string _strawhitslabel;
    // cut variables
    double _edept, _edepl;
    double _rmint, _rminl;
    double _rmaxt, _rmaxl;
    double _maxdt;
    unsigned _maxnpeak;
    unsigned _minnhits;
    // time spectrum parameters
    double _tmin;
    double _tmax;
    double _tbin;
    unsigned _nbins;
    double _ymin;
    double _peakfrac;
    double _tpeakerr;
    double _tdriftmean;
    // used seed fit t0?
    bool _seedt0;
    // outlier cuts
    double _maxseeddoca,_maxkaldoca;
    // cache of event objects
    const StrawHitCollection* _strawhits;
   // Kalman fitters.  Seed fit has a special configuration
    KalFit _seedfit, _kfit;
  // robust helix fitter
    TrkHelixFit _hfit;
   // helper functions
    bool findData(const art::Event& e);
    bool tighthit(double edep, double rho);
    bool loosehit(double edep, double rho);
    void findProximity(std::vector<CLHEP::Hep3Vector> const& shpos, unsigned ish, double maxdist, unsigned& nprox, double& dmin);
    void findPositions(std::vector<CLHEP::Hep3Vector>& shpos);
    void preselectHits(std::vector<CLHEP::Hep3Vector> const& shpos, std::vector<TrkHitFlag>& tflags);
    void findTimePeaks(std::vector<TrkHitFlag> const& tflags, std::vector<TrkTimePeak>& tpeaks);
    void filterOutliers(TrkDef& mytrk,fittype ftype);  
    void createDiagnostics();
    void fillStrawDiag(std::vector<CLHEP::Hep3Vector> const& shpos,std::vector<TrkTimePeak>& tpeaks);
    void fillTimeDiag(unsigned iev,std::vector<TrkHitFlag> const& tflags);
    void fillFitDiag(unsigned ipeak,TrkTimePeak const& tpeak, TrkDef const& helixdef,TrkHelix const& helixfit,
	TrkDef const& seeddef, TrkKalFit const& seedfit, TrkDef const& kaldef, TrkKalFit const& kalfit);
// MC tools
    KalFitMC _kfitmc;
// strawhit tuple variables
    TTree* _shdiag;
    threevec _shpos;
    Float_t _edep;
    Float_t _time;
    Float_t _dmin;
    UInt_t _n50,_n100,_n150,_n200;
    UInt_t _nmcsteps;
    UInt_t _mcnunique,_mcnmax;
    Int_t _mcpdg,_mcgen,_mcproc;
    threevec _mcshpos;
    Float_t _mcedep,_mcemax;
    Float_t _pdist,_pperp,_pmom;
    Float_t _mctime;
    Int_t _loose, _tight;
    UInt_t _ntpeak;
    std::vector<Float_t> _tpeaks;
    std::vector<Int_t> _ntpeaks;
// fit tuple variables
    TTree* _fitdiag;
    Int_t _ipeak;
    Float_t _hcx, _hcy, _hr, _hdfdz, _hfz0;
    Float_t _mccx, _mccy, _mcr, _mcdfdz, _mcfz0;
    Int_t _helixfail,_seedfail,_kalfail;
    helixpar _hpar,_spar;
    helixpar _hparerr,_sparerr;
    Int_t _snhits, _snactive, _sniter, _sndof, _snweediter;
    Float_t _schisq, _st00, _st0;
    UInt_t _nchit;
    UInt_t _npeak, _nmc;
    Float_t _peakmax, _tpeak;
// hit filtering tuple variables
    TTree *_sdiag, *_hdiag;
    Float_t _sresid, _hresid;
    threevec _spos, _hpos;
    Int_t _smcpdg,_smcgen,_smcproc;
    Int_t _hmcpdg,_hmcgen,_hmcproc;
  };

  TrkPatRec::TrkPatRec(fhicl::ParameterSet const& pset) : 
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",10)),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel","makeSH")),
    _edept(pset.get<double>("EDep_tight",0.0045)),
    _edepl(pset.get<double>("EDep_loose",0.005)),
    _rmint(pset.get<double>("RMin_tight",420.0)),
    _rminl(pset.get<double>("RMin_loose",390.0)),
    _rmaxt(pset.get<double>("RMax_tight",630.0)),
    _rmaxl(pset.get<double>("RMax_loose",650.0)),
    _maxdt(pset.get<double>("DtMax",40.0)),
    _maxnpeak(pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits(pset.get<unsigned>("MinNHits",0)),
    _tmin(pset.get<double>("tmin",0.0)),
    _tmax(pset.get<double>("tmax",2000.0)),
    _tbin(pset.get<double>("tbin",20.0)),
    _ymin(pset.get<double>("ymin",5)),
    _peakfrac(pset.get<double>("peakfrac",0.1)),
    _tpeakerr(pset.get<double>("timepeakerr",-8.0)),
    _seedt0(pset.get<bool>("SeedT0",false)),
    _maxseeddoca(pset.get<double>("MaxSeedDoca",180.0)),
    _maxkaldoca(pset.get<double>("MaxKalDoca",180.0)),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    _hfit(pset.get<fhicl::ParameterSet>("HelixFit")),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC"))
  {
    produces<TrkRecoTrkCollection>();
// set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
  }

  TrkPatRec::~TrkPatRec(){}

  void TrkPatRec::beginJob(){
// compute the mean drift time.  This is half the straw radius divided by the drift velocity.
// These numbers should come from conditions objects, FIXME!!!
    _tdriftmean = 0.5*2.5/0.05;
// create diagnostics if requested
    if(_diag > 0)createDiagnostics();
  }

  void TrkPatRec::beginRun(art::Run& ){}

  void TrkPatRec::produce(art::Event& event ) {
// create output
    auto_ptr<TrkRecoTrkCollection> tracks(new TrkRecoTrkCollection );
// event printout
    int iev=event.id().event();
    if((iev%_printfreq)==0)cout<<"TrkPatRec: event="<<iev<<endl;
// find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
// find mc truth if we're making diagnostics
    if(_diag > 0){
      if(!_kfitmc.findMCData(event)){
	cout<<"MC information missing "<< endl;
	return;
      }
    }
//  find hit positions.  This uses conditions data, so it's not an attribute of the hits 
    std::vector<CLHEP::Hep3Vector> shpos;
    findPositions(shpos);
// preselect hits based on internal properties (position, energy, ...)
    std::vector<TrkHitFlag> hitflags;
    hitflags.reserve(_strawhits->size());
    preselectHits(shpos,hitflags);
// find the time peaks in the time spectrum of selected hits
    std::vector<TrkTimePeak> tpeaks;
    findTimePeaks(hitflags,tpeaks);
// fill diagnostics if requested
    if(_diag > 1)fillTimeDiag(iev,hitflags);
    if(_diag > 0)fillStrawDiag(shpos,tpeaks);
// loop over the accepted time peaks
    for(unsigned ipeak=0;ipeak<tpeaks.size();++ipeak){
// track fitting objects for this peak
      TrkHelix helixfit;
      TrkKalFit seedfit, kalfit;
// create track definitions for the different fits from this initial information 
      TrkDef helixdef(_strawhits,tpeaks[ipeak]._trkptrs);
      helixdef.setTrkT0(tpeaks[ipeak]._tpeak-_tdriftmean,_tpeakerr);
      TrkDef seeddef(helixdef);
      TrkDef kaldef(helixdef);
// robust helix fit
      if(_hfit.findHelix(helixdef,helixfit)){
// convert the result to standard helix parameters, and initialize the seed definition helix
	HepVector hpar;
	HepVector hparerr;
	helixfit.helixParams(hpar,hparerr);
	HepSymMatrix hcov = vT_times_v(hparerr);
	seeddef.setHelix(HelixTraj(hpar,hcov));
// Filter outliers using this helix
	filterOutliers(seeddef,seed);
// now, fit the seed helix from the filtered hits
	_seedfit.makeTrack(seeddef,seedfit);
	if(seedfit._fit.success()){
// find the helix parameters from the helix fit, and initialize the full Kalman fit with this
	  double midflt = 0.5*(seedfit._krep->lowFitRange() + seedfit._krep->hiFitRange());
	  double locflt;
	  const HelixTraj* shelix = dynamic_cast<const HelixTraj*>(seedfit._krep->localTrajectory(midflt,locflt));
	  kaldef.setHelix(*shelix);
// filter the outliers
	  filterOutliers(kaldef,kal);
// if requested, use the t0 values from the seed fit.  Otherwise, this is re-computed from the hits
	  if(_seedt0)kaldef.setTrkT0(seedfit._t0);
	  _kfit.makeTrack(kaldef,kalfit);
        }
      }
// fill fit diagnostics if requested
      if(_diag > 0)
	fillFitDiag(ipeak,tpeaks[ipeak],helixdef,helixfit,seeddef,seedfit,kaldef,kalfit);
      if(kalfit._fit.success()){
// save successful kalman fits in the event
	tracks->push_back( kalfit.stealTrack() );
      } else
	kalfit.deleteTrack();
// cleanup the seed fit
      seedfit.deleteTrack();
    }
// put the tracks into the event
    event.put(tracks);
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

  bool TrkPatRec::loosehit(double edep, double rho){
  // looser cuts for pat. rec.
    return edep < _edepl && rho > _rminl && rho < _rmaxl; 
  }

  void
  TrkPatRec::findPositions(std::vector<CLHEP::Hep3Vector>& shpos){
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
      shpos.push_back(pos);
    }
  }

  void
  TrkPatRec::findProximity(std::vector<CLHEP::Hep3Vector> const& shpos, unsigned istr, double maxdist, unsigned& nprox, double& dmin) {
    unsigned nstrs = _strawhits->size();
    dmin=1e6;
    nprox = 0;
    for(unsigned jstr=0;jstr<nstrs;++jstr){
      if(jstr != istr){
        double dist = (shpos[istr]-shpos[jstr]).mag();
        dmin = std::min(dist,dmin);
        if(dist < maxdist)nprox++;
      }
    }
  }

  void
  TrkPatRec::preselectHits(std::vector<CLHEP::Hep3Vector> const& shpos, std::vector<TrkHitFlag>& hitflags){
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      TrkHitFlag flag;
      if(loosehit(sh.energyDep(),shpos[istr].rho()))flag.setLoose();
      if(tighthit(sh.energyDep(),shpos[istr].rho()))flag.setTight();
      hitflags.push_back(flag);
    }
  }

  void 
  TrkPatRec::findTimePeaks(std::vector<TrkHitFlag> const& hitflags, std::vector<TrkTimePeak>& tpeaks) {
    TSpectrum tspec(_maxnpeak);
    TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
// loop over straws hits and fill time spectrum plot for tight hits
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(hitflags[istr].tight()){
	StrawHit const& sh = _strawhits->at(istr);
	double time = sh.time();
        timespec.Fill(time);
      }
    }
    unsigned np = tspec.Search(&timespec,1,"new",_peakfrac);
    Float_t *xpeaks = tspec.GetPositionX();
    Float_t *ypeaks = tspec.GetPositionY();
// Loop over peaks, looking only at those with a minimum peak value
    for (unsigned ip=0; ip<np; ip++) {
      Float_t xp = xpeaks[ip];
      Float_t yp = ypeaks[ip];
      TrkTimePeak tpeak(xp,yp);
      if(yp > _ymin){
// record hits in time with each peak, and accept them if they have a minimum # of hits
        for(unsigned istr=0; istr<nstrs;++istr){
	  if(hitflags[istr].loose()){
	    StrawHit const& sh = _strawhits->at(istr);
	    if(fabs(sh.time()-xp) < _maxdt)tpeak._trkptrs.push_back(istr);
	  }
	}
	if(tpeak._trkptrs.size() > _minnhits)tpeaks.push_back(tpeak);
      }
    }
// sort the peaks so that the largest comes first
    std::sort(tpeaks.begin(),tpeaks.end());
  }

  void
  TrkPatRec::filterOutliers(TrkDef& mytrk,fittype ftype){
// tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// set cuts
    double maxdoca(0.0);
    if(ftype == seed)
      maxdoca = _maxseeddoca;
    else
      maxdoca = _maxkaldoca;
    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const std::vector<size_t>& indices = mytrk.strawHitIndices();
    const HelixTraj& helix = mytrk.helix();
    std::vector<size_t> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = hits->at(indices[ihit]);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      double tddist = tcal->TimeDiffToDistance(straw.index(),sh.dt());
      CLHEP::Hep3Vector pos = straw.getMidPoint() + tddist*straw.getDirection();
  // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(pos.x(),pos.y(),pos.z());
  // estimate fltlen
      double fltlen = (pos.z()-helix.z0())/helix.sinDip();
      TrkPoca ptpoca(helix,fltlen,spt);
  // flag hits with small residuals
      if(ptpoca.doca() < maxdoca){
	goodhits.push_back(indices[ihit]);
      }
  // optional diagnostics
      if(_diag > 0){
  // summarize the MC truth for this strawhit
	PtrStepPointMCVector const& mcptr(_kfitmc.mcData()._mchitptr->at(indices[ihit]));
	std::vector<trksum> mcsum;
	KalFitMC::fillMCSummary(mcptr,mcsum); 
	if(ftype == seed){
	  _spos =pos;
	  _sresid = ptpoca.doca();
	  _smcpdg = mcsum[0]._pdgid;
	  _smcgen = mcsum[0]._gid;
	  _smcproc = mcsum[0]._pid;
	  _sdiag->Fill();
	} else {
	  _hpos =pos;
	  _hresid = ptpoca.doca();
	  _hmcpdg = mcsum[0]._pdgid;
	  _hmcgen = mcsum[0]._gid;
	  _hmcproc = mcsum[0]._pid;
	  _hdiag->Fill();
	}
      }
    }
// update track
    mytrk.setIndices(goodhits);
  }

  void
  TrkPatRec::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
// straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("shpos",&_shpos,"x/F:y/F:z/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("ntpeak",&_ntpeak,"ntpeak/i");
    _shdiag->Branch("tpeaks",&_tpeaks);
    _shdiag->Branch("ntpeaks",&_ntpeaks);
    _shdiag->Branch("dmin",&_dmin,"dmin/F");
    _shdiag->Branch("n50",&_n50,"n50/i");
    _shdiag->Branch("n100",&_n100,"n100/i");
    _shdiag->Branch("n150",&_n150,"n150/i");
    _shdiag->Branch("n200",&_n200,"n200/i");
    _shdiag->Branch("mcshpos",&_mcshpos,"x/F:y/F:z/F");
    _shdiag->Branch("mcedep",&_mcedep,"mcedep/F");
    _shdiag->Branch("mcemax",&_mcemax,"mcemax/F");
    _shdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/i");
    _shdiag->Branch("mcnunique",&_mcnunique,"mcnunique/i");
    _shdiag->Branch("mcnmax",&_mcnmax,"mcnmax/i");
    _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
    _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
    _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
    _shdiag->Branch("mctime",&_mctime,"mctime/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("loose",&_loose,"loose/I");
    _shdiag->Branch("tight",&_tight,"tight/I");
    _shdiag->Branch("pdist",&_pdist,"pdist/F");
    _shdiag->Branch("pperp",&_pperp,"pperp/F");
    _shdiag->Branch("pmom",&_pmom,"pmom/F");
// extend the KalFitMC track diagnostic tuple
    TTree* trkdiag = _kfitmc.createTrkDiag();
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
    trkdiag->Branch("st00",&_st00,"st00/F");
    trkdiag->Branch("st0",&_st0,"st0/F");
    trkdiag->Branch("snhits",&_snhits,"snhits/I");
    trkdiag->Branch("sndof",&_sndof,"sndof/I");
    trkdiag->Branch("sniter",&_sniter,"sniter/i");
    trkdiag->Branch("snweediter",&_snweediter,"snweediter/i");
    trkdiag->Branch("snactive",&_snactive,"snactive/I");
    trkdiag->Branch("schisq",&_schisq,"schisq/F");
    trkdiag->Branch("nchit",&_nchit,"nchit/i");
    trkdiag->Branch("npeak",&_npeak,"npeak/i");
    trkdiag->Branch("tpeak",&_tpeak,"tpeak/F");
    trkdiag->Branch("nmc",&_nmc,"nmc/i");
// seed filtering tuple
    _sdiag = tfs->make<TTree>("sdiag","seed diagnostics");
    _sdiag->Branch("resid",&_sresid,"resid/F");
    _sdiag->Branch("spos",&_spos,"x/F:y/F:z/F");
    _sdiag->Branch("mcpdg",&_smcpdg,"mcpdg/I");
    _sdiag->Branch("mcgen",&_smcgen,"mcgen/I");
    _sdiag->Branch("mcproc",&_smcproc,"mcproc/I");
// helix filtering tuple
    _hdiag = tfs->make<TTree>("hdiag","helix diagnostics");
    _hdiag->Branch("resid",&_hresid,"resid/F");
    _hdiag->Branch("spos",&_hpos,"x/F:y/F:z/F");
    _hdiag->Branch("mcpdg",&_hmcpdg,"mcpdg/I");
    _hdiag->Branch("mcgen",&_hmcgen,"mcgen/I");
    _hdiag->Branch("mcproc",&_hmcproc,"mcproc/I");
  }

  void
  TrkPatRec::fillStrawDiag(std::vector<CLHEP::Hep3Vector> const& shpos,std::vector<TrkTimePeak>& tpeaks) {
    GeomHandle<DetectorSystem> det;
    _nchit = 0;
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      _shpos = shpos[istr];
      _edep = sh.energyDep();
      _time = sh.time();
      // compare to different time peaks
      _ntpeak = tpeaks.size();
      _tpeaks.clear();
      _ntpeaks.clear();
      for(unsigned ipeak=0;ipeak<tpeaks.size();++ipeak){
	_tpeaks.push_back(tpeaks[ipeak]._tpeak);
	_ntpeaks.push_back(tpeaks[ipeak]._trkptrs.size());
      }
      // find proximity for different radii
      double dmin;
      findProximity(shpos,istr,50.0,_n50,dmin);
      findProximity(shpos,istr,100.0,_n100,dmin);
      findProximity(shpos,istr,150.0,_n150,dmin);
      findProximity(shpos,istr,200.0,_n200,dmin);
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
      for( size_t imc=0; imc< _nmcsteps; ++imc ) {
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
      std::vector<trksum> mcsum;
      KalFitMC::fillMCSummary(mcptr,mcsum); 
      _mcnunique = mcsum.size();
      // compute energy sum
      _mcedep = 0.0;
      for(std::vector<trksum>::iterator isum=mcsum.begin(); isum != mcsum.end(); isum++){
	_mcedep += isum->_esum;
      }
      // first entry
      _mcemax = mcsum[0]._esum;
      _mcnmax = mcsum[0]._count;
      _mcpdg = mcsum[0]._pdgid;
      _mcgen = mcsum[0]._gid;
      _mcproc = mcsum[0]._pid;
      _mctime = mcsum[0]._time;
      _mcshpos = mcsum[0]._pos;
      _tight = tighthit(sh.energyDep(),shpos[istr].rho());
      _loose = loosehit(sh.energyDep(),shpos[istr].rho());
      _shdiag->Fill();
      bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
      if(conversion){
	++_nchit;
      }
    }
  }

  void
  TrkPatRec::fillTimeDiag(unsigned iev,std::vector<TrkHitFlag> const& hitflags) {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp;
    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    snprintf(rsname,100,"rawtspectrum%i",iev);
    snprintf(csname,100,"seltspectrum%i",iev);
    snprintf(tsname,100,"convtspectrum%i",iev);
    snprintf(lsname,100,"convtspectrum%i",iev);
    ttsp = tfs->make<TH1F>(tsname,"time spectrum",_nbins,_tmin,_tmax);
    ttsp->SetLineColor(kCyan);
    ltsp = tfs->make<TH1F>(lsname,"time spectrum",_nbins,_tmin,_tmax);
    ltsp->SetLineColor(kGreen);
    rtsp = tfs->make<TH1F>(rsname,"time spectrum",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kBlue);
    ctsp = tfs->make<TH1F>(csname,"time spectrum",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kRed);
    unsigned nstrs = _strawhits->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _strawhits->at(istr);
      double time = sh.time();
      // summarize the MC truth for this strawhit
      PtrStepPointMCVector const& mcptr(_kfitmc.mcData()._mchitptr->at(istr));
      std::vector<trksum> mcsum;
      KalFitMC::fillMCSummary(mcptr,mcsum); 
      bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
      // fill plots
      rtsp->Fill(time);
      if(hitflags[istr].tight())ttsp->Fill(time);
      if(hitflags[istr].loose())ltsp->Fill(time);
      if(conversion)ctsp->Fill(time);
    }
    // find peaks, so they show up on diagnostic plot too
    static TSpectrum tspec(_maxnpeak);
    tspec.Search(ttsp,1,"new",_peakfrac);
  }

  void
  TrkPatRec::fillFitDiag(unsigned ipeak,TrkTimePeak const& tpeak, TrkDef const& helixdef,TrkHelix const& helixfit,
  TrkDef const& seeddef, TrkKalFit const& seedfit, TrkDef const& kaldef, TrkKalFit const& kalfit) {
// convenience numbers
    static const double pi(M_PI);
    static const double twopi(2*pi);
    static const double halfpi(0.5*pi);
// initialize some variables
    _ipeak = ipeak;
    _nmc = 0;
// time peak information
    _peakmax = tpeak._peakmax;
    _tpeak = tpeak._tpeak;
    _npeak = tpeak._trkptrs.size();
// fit status 
    _helixfail = helixfit._fit.failure();
    _seedfail = seedfit._fit.failure();
    _kalfail = kalfit._fit.failure();
// helix information
    if(helixfit._fit.success()){
      HepVector hpar;
      HepVector hparerr;
      helixfit.helixParams(hpar,hparerr);
      _hpar = helixpar(hpar);
      _hparerr = helixpar(hparerr);
      _hcx = helixfit._center.x(); _hcy = helixfit._center.y(); _hr = helixfit._radius;
      _hdfdz = helixfit._dfdz; _hfz0 = helixfit._fz0;
    }
// seed fit information
    if(seedfit._fit.success()){
      _snhits = seeddef.strawHitIndices().size();
      _snactive = seedfit._krep->nActive();
      _sniter = seedfit._krep->iterations();
      _sndof = seedfit._krep->nDof();
      _schisq = seedfit._krep->chisq();
      _st00 = seedfit._t00.t0();
      _st0 = seedfit._t0.t0();
      _snweediter = seedfit._nweediter;
      double loclen;
      const TrkSimpTraj* ltraj = seedfit._krep->localTrajectory(0.0,loclen);
      _spar = helixpar(ltraj->parameters()->parameter());
      _sparerr = helixpar(ltraj->parameters()->covariance());
    }
// use MC truth to define hits and seed helix
    TrkDef mctrk(_strawhits);
    // should be chosing the track ID for conversion a better way, FIXME!!!
    cet::map_vector_key itrk(1);
    if(_kfitmc.trkFromMC(itrk,mctrk)){
      _nmc = 1;
// find true center, radius
      double rtrue = fabs(1.0/mctrk.helix().omega());
      double rad = 1.0/mctrk.helix().omega() + mctrk.helix().d0();
      double cx = -rad*sin(mctrk.helix().phi0());
      double cy = rad*cos(mctrk.helix().phi0());
      _mccx = cx; _mccy = cy; _mcr = rtrue;
      _mcdfdz = fabs(mctrk.helix().omega())/mctrk.helix().tanDip();
      // fix loop for MC values
      _mcfz0 = -mctrk.helix().z0()*mctrk.helix().omega()/mctrk.helix().tanDip() + mctrk.helix().phi0() - halfpi;
      int nloop = (int)rint((helixfit._fz0 - _mcfz0)/twopi);
      _mcfz0 += nloop*twopi;
    }
// fill kalman fit info
    _kfitmc.trkDiag(kalfit);
  }
}

using mu2e::TrkPatRec;
DEFINE_ART_MODULE(TrkPatRec);
