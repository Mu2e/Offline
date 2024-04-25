//
// Combo hit diagnostics.  Split out of MakeComboHits
//
// Original author D. Brown
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"
#include "TTree.h"
// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
// Utilities
#include "Offline/TrkDiag/inc/TrkMCTools.hh"
// diagnostics
#include "Offline/TrkDiag/inc/ComboHitInfo.hh"
#include <map>
#include <limits>

namespace mu2e
{
  class ComboHitDiag : public art::EDAnalyzer {
    public:
      using SPM = std::map<art::Ptr<SimParticle>,unsigned>;

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag"),0 };
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug"),0 };
        fhicl::Atom<bool> mcdiag{ Name("MCDiag"), Comment("MonteCarlo Diag"), true };
        fhicl::Atom<bool> digidiag{ Name("digiDiag"), Comment("Digi Diag"), true };
        fhicl::Atom<art::InputTag> ComboHitCollection{   Name("ComboHitCollection"),   Comment("ComboHit collection name") };
        fhicl::Atom<art::InputTag> StrawDigiCollection{   Name("StrawDigiCollection"),   Comment("StrawDigi collection name") };
        fhicl::Atom<art::InputTag> StrawDigiMCCollection{   Name("StrawDigiMCCollection"),   Comment("StrawDigiMC collection name") };
        fhicl::Atom<art::InputTag> MCPrimary{ Name("MCPrimary"),Comment("MC Primary Particle") };
     };

      explicit ComboHitDiag(const art::EDAnalyzer::Table<Config>& config);
      virtual ~ComboHitDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      // configuration
      int _diag;
      bool _mcdiag, _digidiag;
      // event object Tags
      art::ProductToken<ComboHitCollection> _chToken;
      art::ProductToken<StrawDigiMCCollection> _mcdigisToken;
      art::ProductToken<StrawDigiCollection> _digisToken;
      art::ProductToken<StrawDigiADCWaveformCollection> _digiadcsToken;
      art::ProductToken<PrimaryParticle> _mcprimaryToken;
     // event data cache
      const ComboHitCollection* _chcol;
      const StrawDigiMCCollection *_mcdigis;
      const StrawDigiCollection *_digis;
      const StrawDigiADCWaveformCollection *_digiadcs;
      const PrimaryParticle *_mcprimary;
     // diagnostics
      TTree *_chdiag;
      int _evt; // add event id
      XYZVectorF _pos; // average position
      XYZVectorF _udir; // direction at this position (typically the wire direction)
      XYZVectorF _hdir; // hit direction (generally Z)
      float _wdist; // distance from wire center along this direction
      float _ures, _vres, _wres, _tres, _hcostres, _hphires; // estimated resolution
      float _etime[2]; // end times
      float _ctime; // corrected time
      float _dtime, _ptime; // drift and propagation times: these should be end-specific, TODO
      float _edep; // average energy deposit for these
      float _qual; // quality of combination
      float _dz; // z extent
      float _tot[2]; // tot values
      int _nsh, _nch; // number of associated straw hits
      int _strawid, _straw, _panel, _plane, _level; // strawid info
      int _eend;
      int _esel,_rsel, _tsel, _nsel,  _bkgclust, _bkg, _sth, _dead, _ph, _tdiv, _isolated, _strawxtalk, _elecxtalk, _calosel, _sline;
      // mc diag
      XYZVectorF _mcpos, _mcmom;
      float _mctime, _mcwdist;
      int _mcpdg, _mcproc, _mcgen, _mcndigi;
      int _prel;

      float _mcfrac;

      float _threshold[2], _adcgain;
      std::vector<short unsigned> _digiadc;
      int _digifwpmp;
      int _digipeak;
      float _digipedestal;
      int _digitdc[2], _digitot[2];
      ProditionsHandle<StrawElectronics> _strawele_h;

      // per-hit diagnostics
      std::vector<ComboHitInfo> _chinfo;
      std::vector<ComboHitInfoMC> _chinfomc;
      // helper functions
      bool findData(const art::Event& evt);
  };

  ComboHitDiag::ComboHitDiag(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config),
    _diag( config().diag() ),
    _mcdiag( config().mcdiag() ),
    _digidiag( config().digidiag() ),
    _chToken{ consumes<ComboHitCollection>(config().ComboHitCollection() ) },
    _mcdigisToken{ consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection() ) },
    _digisToken{ consumes<StrawDigiCollection>(config().StrawDigiCollection() ) },
    _digiadcsToken{ consumes<StrawDigiADCWaveformCollection>(config().StrawDigiCollection() ) },
    _mcprimaryToken{ consumes<PrimaryParticle>(config().MCPrimary() ) }
 {}

  ComboHitDiag::~ComboHitDiag(){}

  void ComboHitDiag::beginJob() {
    // create diagnostics if requested
    art::ServiceHandle<art::TFileService> tfs;
    // detailed diagnostics
    _chdiag=tfs->make<TTree>("chdiag","combo hit diagnostics");
    _chdiag->Branch("evt",&_evt,"evt/I");  // add event id
    _chdiag->Branch("wdist",&_wdist,"wdist/F");
    _chdiag->Branch("pos.",&_pos);
    _chdiag->Branch("udir.",&_udir);
    _chdiag->Branch("hdir.",&_hdir);
    _chdiag->Branch("ures",&_ures,"ures/F");
    _chdiag->Branch("vres",&_vres,"vres/F");
    _chdiag->Branch("wres",&_wres,"wres/F");
    _chdiag->Branch("hcostres",&_hcostres,"hcostres/F");
    _chdiag->Branch("hphires",&_hphires,"hphires/F");
    _chdiag->Branch("tres",&_tres,"tres/F");
    _chdiag->Branch("etime",&_etime,"etimecal/F:etimehv");
    _chdiag->Branch("ctime",&_ctime,"ctime/F");
    _chdiag->Branch("dtime",&_dtime,"dtime/F");
    _chdiag->Branch("ptime",&_ptime,"ptime/F");
    _chdiag->Branch("edep",&_edep,"edep/F");
    _chdiag->Branch("qual",&_qual,"qual/F");
    _chdiag->Branch("dz",&_dz,"dz/F");
    _chdiag->Branch("tot",&_tot,"totcal/F:tothv/F");
    _chdiag->Branch("nsh",&_nsh,"nsh/I");
    _chdiag->Branch("nch",&_nch,"nch/I");
    _chdiag->Branch("esel",&_esel,"esel/I");
    _chdiag->Branch("rsel",&_rsel,"rsel/I");
    _chdiag->Branch("tsel",&_tsel,"tsel/I");
    _chdiag->Branch("nsel",&_nsel,"nsel/I");
    _chdiag->Branch("bkgclust",&_bkgclust,"bkgclust/I");
    _chdiag->Branch("bkg",&_bkg,"bkg/I");
    _chdiag->Branch("sth",&_sth,"sth/I");
    _chdiag->Branch("dead",&_dead,"dead/I");
    _chdiag->Branch("ph",&_ph,"ph/I");
    _chdiag->Branch("tdiv",&_tdiv,"tdiv/I");
    _chdiag->Branch("strawxtalk",&_strawxtalk,"strawxtalk/I");
    _chdiag->Branch("elecxtalk",&_elecxtalk,"elecxtalk/I");
    _chdiag->Branch("isolated",&_isolated,"isolated/I");
    _chdiag->Branch("calosel",&_calosel,"calosel/I");
    _chdiag->Branch("sline",&_sline,"sline/I");
    _chdiag->Branch("strawid",&_strawid,"strawid/I");
    _chdiag->Branch("straw",&_straw,"straw/I");
    _chdiag->Branch("panel",&_panel,"panel/I");
    _chdiag->Branch("plane",&_plane,"plane/I");
    _chdiag->Branch("level",&_level,"level/I");
    _chdiag->Branch("earlyend",&_eend,"eend/I");
    if(_diag > 0)
      _chdiag->Branch("chinfo",&_chinfo);
    if(_mcdiag){
      _chdiag->Branch("mcpos.",&_mcpos);
      _chdiag->Branch("mcmom",&_mcmom);
      _chdiag->Branch("mctime",&_mctime,"mctime/F");
      _chdiag->Branch("mcudist",&_mcwdist,"mcudist/F");
      _chdiag->Branch("mcfrac",&_mcfrac,"mcfrac/F");
      _chdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
      _chdiag->Branch("mcproc",&_mcproc,"mcproc/I");
      _chdiag->Branch("mcgen",&_mcgen,"mcgen/I");
      _chdiag->Branch("prel",&_prel,"prel/I");
      _chdiag->Branch("mcndigi",&_mcndigi,"mcndigi/I");
      if(_diag > 0)
        _chdiag->Branch("chinfomc",&_chinfomc);
    }
    if (_digidiag){
      _chdiag->Branch("digitdc",&_digitdc,"digitdccal/I:digitdchv/I");
      _chdiag->Branch("digitot",&_digitot,"digitotcal/I:digitothv/I");
      _chdiag->Branch("digifwpmp",&_digifwpmp,"digifwpmp/I");
      _chdiag->Branch("digipeak",&_digipeak,"digipeak/I");
      _chdiag->Branch("digipedestal",&_digipedestal,"digipedestal/F");
      _chdiag->Branch("digiadc",&_digiadc);
      _chdiag->Branch("threshold",&_threshold,"thresholdcal/F:thresholdhv/F");
      _chdiag->Branch("adcgain",&_adcgain,"adcgain/F");
    }
  }

  void ComboHitDiag::analyze(const art::Event& evt ) {
    // conditions
    StrawElectronics const& strawele = _strawele_h.get(evt.id());
    // find data in event
    findData(evt);
    _evt = evt.id().event();  // add event id
    // loop over combo hits
    for(size_t ich = 0;ich < _chcol->size(); ++ich){
      ComboHit const& ch =(*_chcol)[ich];
      // general information
      _nsh = ch.nStrawHits();
      _nch = ch.nCombo();
      _strawid = ch.strawId().asUint16();
      _straw = ch.strawId().straw();
      _panel = ch.strawId().panel();
      _plane = ch.strawId().plane();
      _level = ch.mask().level();
      _pos = ch.pos();
      _udir = ch.uDir();
      _hdir = ch.hDir();
      _wdist = ch.uPos();
      _ures = sqrt(ch.uVar());
      _vres = sqrt(ch.vVar());
      _hcostres = sqrt(ch.hcostVar());
      _hphires = sqrt(ch.hphiVar());
      _tres = sqrt(ch.timeVar());
      _eend = ch.earlyEnd().end();
      _etime[StrawEnd::cal] = ch.endTime(StrawEnd::cal);
      _etime[StrawEnd::hv] = ch.endTime(StrawEnd::hv);
      _tot[StrawEnd::cal] = ch.TOT(StrawEnd::cal);
      _tot[StrawEnd::hv] = ch.TOT(StrawEnd::hv);
      _ctime = ch.correctedTime();
      _dtime = ch.driftTime();
      _ptime = ch.propTime();
      _edep = ch.energyDep();
      _qual = ch.qual();
      auto const& flag = ch.flag();
      _sth = flag.hasAllProperties(StrawHitFlag::stereo);
      _dead = flag.hasAllProperties(StrawHitFlag::dead);
      _ph = flag.hasAllProperties(StrawHitFlag::panelcombo);
      _tdiv = flag.hasAllProperties(StrawHitFlag::tdiv);
      _esel = flag.hasAllProperties(StrawHitFlag::energysel);
      _rsel = flag.hasAllProperties(StrawHitFlag::radsel);
      _tsel = flag.hasAllProperties(StrawHitFlag::timesel);
      _nsel = flag.hasAllProperties(StrawHitFlag::nhitsel);
      _calosel = flag.hasAllProperties(StrawHitFlag::calosel);
      _sline = flag.hasAllProperties(StrawHitFlag::sline);
      _strawxtalk = flag.hasAllProperties(StrawHitFlag::strawxtalk);
      _elecxtalk = flag.hasAllProperties(StrawHitFlag::elecxtalk);
      _isolated = flag.hasAllProperties(StrawHitFlag::isolated);
      _bkg = flag.hasAllProperties(StrawHitFlag::bkg);
      _bkgclust = flag.hasAllProperties(StrawHitFlag::bkgclust);
      _dz = 0.0;
      // center of this wire
      XYZVectorF cpos = _pos - _wdist*_udir;
      // now hit-by-hit info
      _chinfo.clear();
      if(_diag > 1){
        // loop over components (also ComboHits)
        ComboHitCollection::CHCIter compis;
        _chcol->fillComboHits(ich,compis);
        float minz(1.0e6),maxz(-1.0);
        for(auto compi : compis) {
          ComboHit const& comp = *compi;
          if(comp.pos().z() > maxz) maxz = comp.pos().z();
          if(comp.pos().z() < minz) minz = comp.pos().z();
          ComboHitInfo chi;

          chi._pos= comp.pos();
          chi._udir = comp.uDir();
          auto dp = comp.pos()-ch.pos();
          chi._du = dp.Dot(comp.uDir());
          chi._dv = dp.Dot(comp.vDir());
          chi._dw = dp.Dot(comp.wDir());
          chi._uvar = comp.uVar();
          chi._vvar = comp.vVar();
          chi._tvar = comp.timeVar();
          chi._thit = comp.time();
          chi._straw = comp.strawId().straw();
          chi._upanel = comp.strawId().uniquePanel();
          chi._nch = comp.nCombo();
          chi._nsh = comp.nStrawHits();
          _chinfo.push_back(chi);
        }
        _dz = maxz-minz;
      }
      if(_mcdiag){
        _chinfomc.clear();
        _mcpdg = _mcgen = _mcproc = 0;
        _prel=-1;
        // get the StrawDigi indices associated with this ComboHit
        std::vector<StrawDigiIndex> shids;
        _chcol->fillStrawDigiIndices(ich,shids);
        if(shids.size() != ch.nStrawHits())
          throw cet::exception("DIAG")<<"mu2e::ComboHitDiag: invalid ComboHit Nesting, nshids = "
            << shids.size() << " , n strawhits = " << ch.nStrawHits() << std::endl;
        // find the SimParticle responsable for most of the hits
        SPM spmap;
        for(auto shi : shids) {
          ComboHitInfoMC chimc;
          StrawDigiMC const& mcd = _mcdigis->at(shi);
          auto const& sgsp = mcd.earlyStrawGasStep();
          auto const&  spp = sgsp->simParticle();
          auto fnd = spmap.find(spp);
          if(fnd == spmap.end())
            spmap[spp] = 1;
          else
            ++fnd->second;
        }
        auto ispmax = spmap.begin();
        for(auto isp=spmap.begin();isp != spmap.end(); isp++){
          if(isp->second > ispmax->second) ispmax = isp;
        }
        auto spmax = ispmax->first;
        // use that to define the relationships
        _mcpdg = spmax->pdgId();
        _mcproc = spmax->creationCode();
        _mcfrac = ispmax->second/static_cast<float>(shids.size());
        auto upar = spmax;
        while (upar->genParticle().isNull() && upar->parent().isNonnull()) {
          upar = upar->parent();
        }
        if(upar->genParticle().isNonnull())_mcgen = upar->genParticle()->generatorId().id();
// relationship to MCPrimary
        for(auto const& mcmptr : _mcprimary->primarySimParticles()){
          MCRelationship rel(mcmptr,spmax);
          if(rel.relationship() > MCRelationship::none){
            if(_prel > MCRelationship::none)
              _prel = std::min(_prel,(int)rel.relationship());
            else
              _prel = rel.relationship();
          }
        }
        // find the relation with each individual straw hit
        // average the primary MC particle properties
        _mcpos = _mcmom = XYZVectorF();
        _mctime = 0.0;
        unsigned npri=0;
        for(auto shi : shids) {
          ComboHitInfoMC chimc;
          StrawDigiMC const& mcd = _mcdigis->at(shi);
          MCRelationship rel(mcd,spmax);
          chimc._rel = rel.relationship();
          _chinfomc.push_back(chimc);
          if(chimc._rel==0){
            npri++;
            auto const& sgsp = mcd.earlyStrawGasStep();
            _mcmom += sgsp->momentum();
            _mcpos += 0.5*(sgsp->startPosition() + sgsp->endPosition());
            _mctime += sgsp->time();
          }
        }
        if(npri > 0){
          _mcmom /= npri;
          _mcpos /= npri;
          _mctime /= npri;
        }
        _mcwdist = (_mcpos - cpos).Dot(_udir);
        // count mcdigis with the same StrawId as this hit
        _mcndigi = 0;
        for(auto const& mcdigi : *_mcdigis) {
          if(ch._mask.equal(mcdigi.strawId(),ch.strawId()) && spmax == mcdigi.earlyStrawGasStep()->simParticle())++_mcndigi;
        }
      }
      if (_digis != 0){
        std::vector<StrawDigiIndex> shids;
        _chcol->fillStrawDigiIndices(ich,shids);
        // use the 1st hit to define the MC match; this is arbitrary should be an average FIXME!
        auto digi = _digis->at(shids[0]);
        for(size_t iend=0;iend<2;++iend){
          _digitdc[iend] = digi.TDC(static_cast<StrawEnd::End>(iend));
          _digitot[iend] = digi.TOT(static_cast<StrawEnd::End>(iend));
        }
        _digifwpmp = digi.PMP();

        for (size_t iend=0;iend<2;++iend){
          _threshold[iend] = strawele.threshold(ch.strawId(), static_cast<StrawEnd::End>(iend));
        }
        _adcgain = strawele.currentToVoltage(ch.strawId(),StrawElectronics::adc);
      }
      if (_digiadcs != 0){
        std::vector<StrawDigiIndex> shids;
        _chcol->fillStrawDigiIndices(ich,shids);
        // use the 1st hit to define the MC match; this is arbitrary should be an average FIXME!
        auto digiadc = _digiadcs->at(shids[0]);
        _digiadc = digiadc.samples();
        _digipedestal = 0;
        _digipeak = 0;
        for (size_t i=0;i<_digiadc.size();++i){
          if (i<strawele.nADCPreSamples()){
            _digipedestal += _digiadc[i];
          }
          if (_digiadc[i] > _digipeak)
            _digipeak = _digiadc[i];
        }
        _digipedestal /= (float) strawele.nADCPreSamples();
      }
      _chdiag->Fill();
    }
  }

  bool ComboHitDiag::findData(const art::Event& evt){
    _chcol = 0; _mcdigis = 0; _digis = 0; _digiadcs = 0;
    auto chH = evt.getValidHandle(_chToken);
    _chcol = chH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle(_mcdigisToken);
      _mcdigis = mcdH.product();
       auto mcpH = evt.getValidHandle(_mcprimaryToken);
      _mcprimary = mcpH.product();
   }
    if (_digidiag){
      auto dH = evt.getValidHandle(_digisToken);
      _digis = dH.product();
      auto daH = evt.getValidHandle(_digiadcsToken);
      _digiadcs = daH.product();
    }
    return _chcol != 0 && ( (_mcdigis != 0 && _mcprimary != 0)  || !_mcdiag)
      && ((_digis != 0 && _digiadcs != 0) || !_digidiag);
  }
}
// Part of the magic that makes this class a module.
using mu2e::ComboHitDiag;
DEFINE_ART_MODULE(ComboHitDiag)
