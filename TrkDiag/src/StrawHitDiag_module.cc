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
#include "art_root_io/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
// data
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "DataProducts/inc/EventWindowMarker.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
using namespace std;
using CLHEP::Hep3Vector;

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
      bool _mcdiag, _useshfcol;
  // data tags
      art::InputTag _shTag;
      art::InputTag _chTag;
      art::InputTag _shfTag;
      art::InputTag _ewMarkerTag;
      art::InputTag _stTag;
      art::InputTag _mcdigisTag;
      // cache of event objects
      const StrawHitCollection* _shcol;
      const ComboHitCollection* _chcol;
      const StrawHitFlagCollection* _shfcol;
      const StrawDigiMCCollection *_mcdigis;
      // time offset
      SimParticleTimeOffset _toff;
      // strawhit tuple variables
      TTree *_shdiag;
      Int_t _eventid, _subrunid, _runid;
      Hep3Vector _shp;
      Float_t _shlen, _slen;
      Float_t _edep;
      Float_t _time[2], _tot[2];
      Float_t _correcttime, _ctime, _dtime, _ptime;
      Float_t _rho;
      Int_t _mcnsteps;
      Int_t _mcpdg,_mcgen,_mcproc, _mcid;
      Int_t _mcppdg,_mcpproc;
      Int_t _mcgid, _mcgpdg;
      Float_t _mcge, _mcgt;
      Hep3Vector _mcshp, _mcop, _mcpop, _mcgpos;
      Float_t _mcoe, _mcpoe, _mcom, _mcpom;
      Float_t _mcshlen,_mcshd, _mcsphi, _mcplen;
      Float_t _mcedep, _mcetrig;
      Float_t _mcct[2], _mccphi[2], _mccd[2];
      Float_t _pdist,_pperp,_pmom;
      Float_t _mcwt[2];
      Double_t _mcsptime;
      Double_t _mcptime;
      Int_t _esel,_rsel, _tsel,  _bkgclust, _bkg, _stereo, _tdiv, _isolated, _strawxtalk, _elecxtalk, _calosel;
      Int_t _sid, _plane, _panel, _layer, _straw;
      Float_t _shwres, _shtres;
      Bool_t _mcxtalk;
      Float_t _ewMarkerOffset;
      // helper array
      StrawEnd _end[2];
  };

  StrawHitDiag::StrawHitDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _useshfcol(pset.get<bool>("UseStrawHitFlagCollection",true)),
    _shTag(pset.get<string>("StrawHitCollection","makeSH")),
    _chTag(pset.get<string>("ComboHitCollection","makeSH")),
    _shfTag(pset.get<string>("StrawHitFlagCollection")),
    _ewMarkerTag(pset.get<art::InputTag>("EventWindowMarkerLabel","EWMProducer")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets")),
    _end{StrawEnd::cal,StrawEnd::hv}
  {
    if(pset.get<bool>("TestStrawId",false)) {
      for(uint16_t plane = 0; plane < StrawId::_nplanes; ++plane){
	StrawId sid(plane,0,0);
	std::cout << "Plane StrawId " << sid.asUint16() << " plane " << sid.plane() << std::endl;
	for(uint16_t panel = 0; panel < StrawId::_npanels; ++panel){
	  StrawId sid(plane,panel,0);
	  std::cout << "Panel StrawId " << sid.asUint16() << " panel " << sid.uniquePanel() << std::endl;
	  for(uint16_t straw = 0; straw < StrawId::_nstraws; ++straw){
	    StrawId sid(plane,panel,straw);
	      std::cout << "Straw StrawId " << sid.asUint16() << " unique straw " << sid.uniqueStraw() << std::endl;
	  }
	}
      }
    }
  }

  StrawHitDiag::~StrawHitDiag(){}

  void StrawHitDiag::analyze(const art::Event& event ) {
    _eventid = event.event();
    _subrunid = event.subRun();
    _runid = event.run();
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkPatRec: data missing or incomplete"<< endl;
    }
    fillStrawHitDiag();
  }

  bool StrawHitDiag::findData(const art::Event& evt){
    _shcol = 0;
    _chcol = 0;
    _shfcol = 0;
    _mcdigis = 0;
    _ewMarkerOffset = 0;
    // nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
    _chcol = chH.product();
    if(_useshfcol){
      auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
      _shfcol = shfH.product();
    }
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      // update time offsets
      _toff.updateMap(evt);
    }
    auto ewMarkerHandle = evt.getValidHandle<EventWindowMarker>(_ewMarkerTag);
    auto ewMarker = ewMarkerHandle.product();
    _ewMarkerOffset = ewMarker->timeOffset();
    return _shcol != 0 && _chcol != 0 && (_shfcol != 0 || !_useshfcol) && (_mcdigis != 0  || !_mcdiag);
  }

  void StrawHitDiag::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    // straw hit tuple
    _shdiag=tfs->make<TTree>("shdiag","strawhit diagnostics");
    _shdiag->Branch("eventid",&_eventid,"eventid/I");
    _shdiag->Branch("subrunid",&_subrunid,"subrunid/I");
    _shdiag->Branch("runid",&_runid,"runid/I");
    _shdiag->Branch("shpos.",&_shp);
    _shdiag->Branch("shlen",&_shlen,"shlen/F");
    _shdiag->Branch("slen",&_slen,"slen/F");
    _shdiag->Branch("edep",&_edep,"edep/F");
    _shdiag->Branch("time",&_time,"tcal/F:thv/F");
    _shdiag->Branch("ctime",&_ctime,"ctime/F");
    _shdiag->Branch("dtime",&_dtime,"dtime/F");
    _shdiag->Branch("ptime",&_ptime,"ptime/F");
    _shdiag->Branch("correcttime",&_correcttime,"correcttime/F");
    _shdiag->Branch("tot",&_tot,"totcal/F:tothv/F");
    _shdiag->Branch("rho",&_rho,"rho/F");
    _shdiag->Branch("sid",&_sid,"sid/I");
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
    if(_mcdiag){
      _shdiag->Branch("mcshpos.",&_mcshp);
      _shdiag->Branch("mcopos.",&_mcop);
      _shdiag->Branch("mcpopos.",&_mcpop);
      _shdiag->Branch("mcct",&_mcct,"mcctcal/F:mccthv/F");
      _shdiag->Branch("mcoe",&_mcoe,"mcoe/F");
      _shdiag->Branch("mcom",&_mcom,"mcom/F");
      _shdiag->Branch("mcpoe",&_mcpoe,"mcpoe/F");
      _shdiag->Branch("mcpom",&_mcpom,"mcpom/F");
      _shdiag->Branch("mcshlen",&_mcshlen,"mcshlen/F");
      _shdiag->Branch("mcshd",&_mcshd,"mcshd/F");
      _shdiag->Branch("mcsphi",&_mcsphi,"mcsphi/F");
      _shdiag->Branch("mccphi",&_mccphi,"mccphical/F:,mccphihv/F");
      _shdiag->Branch("mccd",&_mccd,"mccdcal/F:,mccdhv/F");
      _shdiag->Branch("mcplen",&_mcplen,"mcplen/F");
      _shdiag->Branch("mcedep",&_mcedep,"mcedep/F");
      _shdiag->Branch("mcetrig",&_mcetrig,"mcetrig/F");
      _shdiag->Branch("mcnsteps",&_mcnsteps,"mcnsteps/I");
      _shdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
      _shdiag->Branch("mcgen",&_mcgen,"mcgen/I");
      _shdiag->Branch("mcproc",&_mcproc,"mcproc/I");
      _shdiag->Branch("mcsptime",&_mcsptime,"mcsptime/D");
      _shdiag->Branch("mcwt",&_mcwt,"mcwtcal/F:mcwthv/F");
      _shdiag->Branch("mcppdg",&_mcppdg,"mcppdg/I");
      _shdiag->Branch("mcpproc",&_mcpproc,"mcpproc/I");
      _shdiag->Branch("mcptime",&_mcptime,"mcptime/D");
      _shdiag->Branch("mcgid",&_mcgid,"mcgid/I");
      _shdiag->Branch("mcid",&_mcid,"mcid/I");
      _shdiag->Branch("mcgpdg",&_mcgpdg,"mcgpdg/I");
      _shdiag->Branch("mcge",&_mcge,"mcge/F");
      _shdiag->Branch("mcgt",&_mcgt,"mcgt/F");
      _shdiag->Branch("mcgpos.",&_mcgpos);
      _shdiag->Branch("mcxtalk",&_mcxtalk,"mcxtalk/B");
    }
  }

  void StrawHitDiag::fillStrawHitDiag() {
    GeomHandle<DetectorSystem> det;
    const Tracker& tracker = *GeomHandle<Tracker>();
    static const double rstraw = tracker.getStraw(StrawId(0,0,0)).getRadius();
    unsigned nstrs = _chcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawHit const& sh = _shcol->at(istr);
      ComboHit const& ch = _chcol->at(istr);
      StrawHitFlag shf = ch.flag();
      if(_useshfcol) shf.merge(_shfcol->at(istr));
      const Straw& straw = tracker.getStraw( ch.strawId() );
      _sid = straw.id().asUint16();
      _plane = straw.id().getPlane();
      _panel = straw.id().getPanel();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();
      _edep = ch.energyDep();
      for(size_t iend=0;iend<2;++iend){
	_time[iend] = sh.time(_end[iend]);
	_tot[iend] = sh.TOT(_end[iend]);
      }
      _correcttime = ch.correctedTime();
      _ctime = ch.time();
      _dtime = ch.driftTime();
      _ptime = ch.propTime();
      _shp = ch.posCLHEP();
      _shlen =(ch.posCLHEP()-straw.getMidPoint()).dot(straw.getDirection());
      _slen = straw.halfLength();
      _stereo = ch.flag().hasAllProperties(StrawHitFlag::stereo);
      _tdiv = ch.flag().hasAllProperties(StrawHitFlag::tdiv);
      _esel = shf.hasAllProperties(StrawHitFlag::energysel);
      _rsel = shf.hasAllProperties(StrawHitFlag::radsel);
      _tsel = shf.hasAllProperties(StrawHitFlag::timesel);
      _calosel = shf.hasAllProperties(StrawHitFlag::calosel);
      _strawxtalk = shf.hasAllProperties(StrawHitFlag::strawxtalk);
      _elecxtalk = shf.hasAllProperties(StrawHitFlag::elecxtalk);
      _isolated = shf.hasAllProperties(StrawHitFlag::isolated);
      _bkg = shf.hasAllProperties(StrawHitFlag::bkg);
      _bkgclust = shf.hasAllProperties(StrawHitFlag::bkgclust);
      _rho = ch.posCLHEP().perp();
      // summarize the MC truth for this strawhit.  Preset the values in case MC is missing/incomplete
      _mcgid = -1;
      _mcid = -1;
      _mcgpdg = -1;
      _mcge = -1.0;
      _mcgt = -1.0;
      _mcgpos = Hep3Vector();
      _mcppdg=0;
      _mcpproc=-1;
      _mcptime=0.0;
      _mcpop = Hep3Vector();
      _mcpoe = _mcpom = -1.0;
      _mcnsteps = -1;
      _mcpdg = -1;
      _mcgen = -1;
      _mcproc = -1;
      _mcsptime = -1.0;
      _mcwt[0] = _mcwt[1] = -1.0;
      _mcshp = Hep3Vector();
      _mcop = Hep3Vector();
      _mcoe = -1;
      _mcom = -1;
      _mcshlen = -1;
      _mcshd = -1;
      _mcsphi = 0.0;
      _mccphi[0] = _mccphi[1] = 0.0;
      _mccd[0] = _mccd[1] = -1;
      _mcplen = -1;
      _mcpproc=-1;
      _mcptime=0.0;
      _mcpoe = _mcpom = -1.0;
      _mcxtalk = false;
      if(_mcdigis != 0){
        StrawDigiMC const& mcdigi = _mcdigis->at(istr);
        // use TDC channel 0 to define the MC match
        StrawEnd itdc;
        auto const& spmcp = mcdigi.strawGasStep(itdc);
        art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	SimParticle const& osp = spp->originParticle();
	Hep3Vector dprod = spmcp->position()-det->toDetector(osp.startPosition());
	static Hep3Vector zdir(0.0,0.0,1.0);
        _pdist = dprod.mag();
        _pperp = dprod.perp(zdir);
        _pmom = sqrt(spmcp->momentum().mag2());
        _mcnsteps = 2; // FIXME!
        // compute energy sum
        _mcedep = mcdigi.energySum();
        _mcetrig = mcdigi.triggerEnergySum(StrawEnd::cal);
        _mcpdg = osp.pdgId();
        _mcproc = osp.creationCode();
        _mcgen = -1;
        if(osp.genParticle().isNonnull())
          _mcgen = osp.genParticle()->generatorId().id();
        _mcsptime = _toff.timeWithOffsetsApplied(*spmcp) - _ewMarkerOffset;
	for(size_t iend=0;iend<2; ++iend){
	  _mcwt[iend] = mcdigi.wireEndTime(_end[iend]);
	  _mcct[iend] = mcdigi.clusterPosition(_end[iend]).t();
          Hep3Vector cpos = mcdigi.clusterPosition(_end[iend]).vect();
          Hep3Vector cdir = (cpos-straw.getMidPoint());
          cdir -= straw.getDirection()*(cdir.dot(straw.getDirection()));
          _mccphi[iend] = cdir.theta();
          _mccd[iend] = min(cdir.perp(straw.getDirection()),straw.innerRadius());
	}
        _mcshp = spmcp->position();
        _mcop = det->toDetector(osp.startPosition());
        _mcoe = osp.startMomentum().e();
        _mcom = osp.startMomentum().vect().mag();
        _mcshlen = (spmcp->position()-straw.getMidPoint()).dot(straw.getDirection());
	Hep3Vector mdir = Geom::Hep3Vec(spmcp->momentum()).unit();
	Hep3Vector tdir = (straw.getDirection().cross(mdir)).unit();
        _mcshd = (spmcp->position()-straw.getMidPoint()).dot(tdir);
	double scos = mdir.dot(straw.getDirection());
        _mcplen = 2.0*sqrt( (rstraw*rstraw -_mcshd*_mcshd)/(1.0-scos*scos) );
	_mcsphi = atan2(tdir.perp(),tdir.z()); // 'azimuth' around the straw of the POCA
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
	  _mcid = spp->id().asInt();
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
        _mcxtalk = spmcp->strawId() != sh.strawId();
      }
      _shwres = _chcol->at(istr).posRes(ComboHit::wire);
      _shtres = _chcol->at(istr).posRes(ComboHit::trans);
//  Info depending on stereo hits
      _shdiag->Fill();
    }
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::StrawHitDiag);
