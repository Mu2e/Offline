//
// Straw hit diagnostics.  Split out of TrkPatRec
//
// Original author D. Brown
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// root 
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// diagnostics
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "DataProducts/inc/threevec.hh"
using namespace std; 

namespace mu2e 
{
  class StrawHitDiag : public art::EDAnalyzer {
    public:
      explicit StrawHitDiag(fhicl::ParameterSet const&);
      virtual ~StrawHitDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);

    private:
      // helper functions
      void fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const;
      void fillStrawHitDiag();
      bool findData(const art::Event& e);
      // event object labels
      string _shLabel;
      string _shpLabel;
      string _shfLabel;
      string _stLabel;
      std::string _mcdigislabel;
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StereoHitCollection* _stcol;
      const StrawHitFlagCollection* _shfcol;
      const StrawDigiMCCollection *_mcdigis;

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
  };


  StrawHitDiag::StrawHitDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel(pset.get<string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _stLabel(pset.get<string>("StereoHitCollectionLabel","MakeStereoHits")),
    _mcdigislabel(pset.get<string>("StrawHitMCLabel","makeSH")) 
  {}

  StrawHitDiag::~StrawHitDiag(){}

  void StrawHitDiag::analyze(const art::Event& event ) {
    _eventid = event.event();
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: data missing or incomplete"<< endl;
    }
    fillStrawHitDiag();
  }

  void StrawHitDiag::fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const {
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

    if(_mcdigis != 0) {

      StrawDigiMC const& mcdigi = _mcdigis->at(ish);
      // use TDC channel 0 to define the MC match
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
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

  bool StrawHitDiag::findData(const art::Event& evt){
    _shcol = 0;
    _shpcol = 0;
    _shfcol = 0;
    _stcol = 0;
    _mcdigis = 0;

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(evt.getByLabel(_shLabel,strawhitsH))
      _shcol = strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
    art::Handle<mu2e::StereoHitCollection> stH;
    if(evt.getByLabel(_stLabel,stH))
      _stcol = stH.product();
    art::Handle<StrawDigiMCCollection> mcdigisHandle;
    if(evt.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle))
      _mcdigis = mcdigisHandle.product();
// don't require stereo hits or MC
    return _shcol != 0 && _shpcol != 0&& _shfcol != 0;
  }

  void StrawHitDiag::beginJob(){
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
  }


  void StrawHitDiag::fillStrawHitDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = getTrackerOrThrow();
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
      _mcpproc=-1;
      _mcptime=0.0;
      _mcpop = threevec();
      _mcpoe = _mcpom = -1.0;
      _xtalk = false;
      if(_mcdigis != 0){
        StrawDigiMC const& mcdigi = _mcdigis->at(istr);
        // use TDC channel 0 to define the MC match
        StrawDigi::TDCChannel itdc = StrawDigi::zero;
        if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
        art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
        art::Ptr<SimParticle> const& spp = spmcp->simParticle();
        CLHEP::Hep3Vector dprod = spmcp->position()-det->toDetector(spp->startPosition());
        static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
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
      _esel = _shfcol->at(istr).hasAllProperties(StrawHitFlag::energysel);
      _rsel = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      _timesel = _shfcol->at(istr).hasAllProperties(StrawHitFlag::timesel);
      _stereo = _shfcol->at(istr).hasAllProperties(StrawHitFlag::stereo);
      _isolated = _shfcol->at(istr).hasAllProperties(StrawHitFlag::isolated);
      _delta = _shfcol->at(istr).hasAllProperties(StrawHitFlag::delta);
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
  }
}  // end namespace mu2e

// Part of the magic that makes this class a module.
using mu2e::StrawHitDiag;
DEFINE_ART_MODULE(StrawHitDiag);


