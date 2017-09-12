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
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// diagnostics
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
      void fillStrawHitDiag();
      bool findData(const art::Event& e);
      // control flags
      bool _mcdiag;
  // data tags
      art::InputTag _shTag;
      art::InputTag _shpTag;
      art::InputTag _shfTag;
      art::InputTag _stTag;
      art::InputTag _mcdigisTag;
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StereoHitCollection* _stcol;
      const StrawHitFlagCollection* _shfcol;
      const StrawDigiMCCollection *_mcdigis;
      // time offset
      SimParticleTimeOffset _toff;
      // strawhit tuple variables
      TTree *_shdiag;
      Int_t _eventid;
      threevec _shp;
      Float_t _shlen, _slen; 
      Float_t _edep;
      Float_t _time, _bkgt, _rho;
      Int_t _mcnunique,_mcnmax;
      Int_t _mcpdg,_mcgen,_mcproc;
      Int_t _mcppdg,_mcpproc;
      Int_t _mcgid, _mcgpdg;
      Float_t _mcge, _mcgt;
      threevec _mcshp, _mcop, _mcpop, _mcgpos;
      Float_t _mcoe, _mcpoe, _mcom, _mcpom;
      Float_t _mcshlen,_mcshd;
      Float_t _mcedep, _mcetrig;
      Float_t _pdist,_pperp,_pmom;
      Double_t _mctime, _mcptime;
      Int_t _esel,_rsel, _tsel,  _bkgclust, _bkg, _stereo, _tdiv, _isolated, _strawxtalk, _elecxtalk, _calosel;
      Int_t _plane, _panel, _layer, _straw;
      Float_t _shwres, _shtres, _shchisq, _shdt, _shdist;
      Bool_t _mcxtalk;
  };

  StrawHitDiag::StrawHitDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _shTag(pset.get<string>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollectionTag","MakeStereoHits")),
    _shfTag(pset.get<string>("StrawHitFlagCollectionTag","FlagBkgHits")),
    _stTag(pset.get<string>("StereoHitCollectionTag","MakeStereoHits")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {}

  StrawHitDiag::~StrawHitDiag(){}

  void StrawHitDiag::analyze(const art::Event& event ) {
    _eventid = event.event();
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: data missing or incomplete"<< endl;
    }
    fillStrawHitDiag();
  }
  bool StrawHitDiag::findData(const art::Event& evt){
    _shcol = 0;
    _shpcol = 0;
    _shfcol = 0;
    _stcol = 0;
    _mcdigis = 0;
    // nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto hsH = evt.getValidHandle<StereoHitCollection>(_stTag);
    _stcol = hsH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      // update time offsets
      _toff.updateMap(evt);
    }
    return _shcol != 0 && _shpcol != 0 && _shfcol != 0 
      && (_mcdigis != 0  || !_mcdiag);
  }

  void StrawHitDiag::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    // straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("shpos",&_shp,"x/F:y/F:z/F");
    _shdiag->Branch("shlen",&_shlen,"shlen/F");
    _shdiag->Branch("slen",&_slen,"slen/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"time/F");
    _shdiag->Branch("deltat",&_bkgt,"deltat/F");
    _shdiag->Branch("rho",&_rho,"rho/F");
    _shdiag->Branch("plane",&_plane,"plane/I");
    _shdiag->Branch("panel",&_panel,"panel/I");
    _shdiag->Branch("layer",&_layer,"layer/I");
    _shdiag->Branch("straw",&_straw,"straw/I");
    _shdiag->Branch("esel",&_esel,"esel/I");
    _shdiag->Branch("rsel",&_rsel,"rsel/I");
    _shdiag->Branch("tsel",&_tsel,"tsel/I");
    _shdiag->Branch("bkgclust",&_bkgclust,"bkgclust/I");
    _shdiag->Branch("bkg",&_bkg,"bkg/I");
    _shdiag->Branch("stereo",&_stereo,"stereo/I");
    _shdiag->Branch("tdiv",&_tdiv,"tdiv/I");
    _shdiag->Branch("strawxtalk",&_strawxtalk,"strawxtalk/I");
    _shdiag->Branch("elecxtalk",&_elecxtalk,"elecxtalk/I");
    _shdiag->Branch("isolated",&_isolated,"isolated/I");
    _shdiag->Branch("calosel",&_calosel,"calosel/I");
    _shdiag->Branch("pdist",&_pdist,"pdist/F");
    _shdiag->Branch("pperp",&_pperp,"pperp/F");
    _shdiag->Branch("pmom",&_pmom,"pmom/F");
    _shdiag->Branch("wres",&_shwres,"wres/F");
    _shdiag->Branch("tres",&_shtres,"tres/F");
    _shdiag->Branch("shchisq",&_shchisq,"shchisq/F");
    _shdiag->Branch("shdt",&_shdt,"shdt/F");
    _shdiag->Branch("shdist",&_shdist,"shdist/F");
    if(_mcdiag){
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
      _shdiag->Branch("mcetrig",&_mcetrig,"mcetrig/F");
      _shdiag->Branch("mcnunique",&_mcnunique,"mcnunique/I");
      _shdiag->Branch("mcnmax",&_mcnmax,"mcnmax/I");
      _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
      _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
      _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
      _shdiag->Branch("mctime",&_mctime,"mctime/D");
      _shdiag->Branch("mcppdg",&_mcppdg,"mcpdg/I");
      _shdiag->Branch("mcpproc",&_mcpproc,"mcpproc/I");
      _shdiag->Branch("mcptime",&_mcptime,"mcptime/D");
      _shdiag->Branch("mcgid",&_mcgid,"mcgid/I");
      _shdiag->Branch("mcgpdg",&_mcgpdg,"mcgpdg/I");
      _shdiag->Branch("mcge",&_mcge,"mcge/F");
      _shdiag->Branch("mcgt",&_mcgt,"mcgt/F");
      _shdiag->Branch("mcgpos",&_mcgpos,"x/F:y/F:z/F");
      _shdiag->Branch("mcxtalk",&_mcxtalk,"mcxtalk/B");
    }
  }

  void StrawHitDiag::fillStrawHitDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = getTrackerOrThrow();
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _shcol->at(istr);
      StrawHitPosition const& shp = _shpcol->at(istr);
      StrawHitFlag const& shf = _shfcol->at(istr);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      _plane = straw.id().getPlane();
      _panel = straw.id().getPanel();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();
      _edep = sh.energyDep();
      _time = sh.time();
      _bkgt = sh.dt();
      _shp = shp.pos();
      _shlen =(shp.pos()-straw.getMidPoint()).dot(straw.getDirection());
      _slen = straw.getHalfLength(); 
      _stereo = shp.flag().hasAllProperties(StrawHitFlag::stereo);
      _tdiv = shp.flag().hasAllProperties(StrawHitFlag::tdiv);
      _esel = shf.hasAllProperties(StrawHitFlag::energysel);
      _rsel = shf.hasAllProperties(StrawHitFlag::radsel);
      _tsel = shf.hasAllProperties(StrawHitFlag::timesel);
      _calosel = shf.hasAllProperties(StrawHitFlag::calosel);
      _strawxtalk = shf.hasAllProperties(StrawHitFlag::strawxtalk);
      _elecxtalk = shf.hasAllProperties(StrawHitFlag::elecxtalk);
      _isolated = shf.hasAllProperties(StrawHitFlag::isolated);
      _bkg = shf.hasAllProperties(StrawHitFlag::bkg);
      _bkgclust = shf.hasAllProperties(StrawHitFlag::bkgclust);
      _rho = shp.pos().perp();
      // summarize the MC truth for this strawhit.  Preset the values in case MC is missing/incomplete
      _mcgid = -1;
      _mcgpdg = -1;
      _mcge = -1.0;
      _mcgt = -1.0;
      _mcgpos = threevec();
      _mcppdg=0;
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
      _mcxtalk = false;
      if(_mcdigis != 0){
        StrawDigiMC const& mcdigi = _mcdigis->at(istr);
        // use TDC channel 0 to define the MC match
        StrawDigi::TDCChannel itdc = StrawDigi::zero;
        if(!mcdigi.hasTDC(StrawDigi::zero)) itdc = StrawDigi::one;
        art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
        art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	SimParticle const& osp = spp->originParticle();
	CLHEP::Hep3Vector dprod = spmcp->position()-det->toDetector(osp.startPosition());
	static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
        _pdist = dprod.mag();
        _pperp = dprod.perp(zdir);
        _pmom = spmcp->momentum().mag();
        _mcnunique = mcdigi.stepPointMCs().size();
        // compute energy sum
        _mcedep = mcdigi.energySum();
        _mcetrig = mcdigi.triggerEnergySum();
        _mcnmax = mcdigi.stepPointMCs().size();
        _mcpdg = osp.pdgId();
        _mcproc = osp.creationCode();
        _mcgen = -1;
        if(osp.genParticle().isNonnull())
          _mcgen = osp.genParticle()->generatorId().id();
        _mctime = _toff.timeWithOffsetsApplied(*spmcp);
        _mcshp = spmcp->position();
        _mcop = det->toDetector(osp.startPosition());
        _mcoe = osp.startMomentum().e();
        _mcom = osp.startMomentum().vect().mag();
        _mcshlen = (spmcp->position()-straw.getMidPoint()).dot(straw.getDirection());
        _mcshd = (spmcp->position()-straw.getMidPoint()).dot(straw.getDirection().cross(spmcp->momentum().unit()));
  // immediate parent information
	art::Ptr<SimParticle> psp = osp.parent();
	if(psp.isNonnull()){
	  SimParticle const& posp =psp->originParticle();
          _mcppdg = posp.pdgId();
          _mcpproc = posp.creationCode();
          _mcptime = _toff.totalTimeOffset(psp) + psp->startGlobalTime();
          _mcpop = det->toDetector(posp.startPosition());
          _mcpoe = posp.startMomentum().e();
          _mcpom = posp.startMomentum().vect().mag();
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
        _mcxtalk = spmcp->strawIndex() != sh.strawIndex();

      }
      _shwres = _shpcol->at(istr).posRes(StrawHitPosition::wire);
      _shtres = _shpcol->at(istr).posRes(StrawHitPosition::trans);
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


