//
// Stereo hit diagnostics.  Split out of MakeComboHits
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
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
// diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
using namespace std; 

namespace mu2e 
{
 class ComboHitDiag : public art::EDAnalyzer {
    public:
      explicit ComboHitDiag(fhicl::ParameterSet const&);
      virtual ~ComboHitDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      // configuration
      int _diag;
      bool _mcdiag;
      // event object Tags
      art::InputTag   _shtag;
      art::InputTag   _shptag;
      art::InputTag   _chtag;
      art::InputTag   _mcdigistag;
      // event data cache
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const ComboHitCollection* _chcol;
      const StrawDigiMCCollection *_mcdigis;
        // time offset
      SimParticleTimeOffset _toff;
      // diagnostics
      TTree *_chdiag;
      XYZVec _pos; // average position
      XYZVec _wdir; // direction at this position (typically the wire direction)
      Float_t _wdist; // distance from wire center along this direction
      Float_t _wres; // estimated error along this direction
      Float_t _tres; // estimated error perpendicular to this direction
      Float_t _time; // Average time for these
      Float_t _edep; // average energy deposit for these
      Float_t _qual; // quality of combination
      Int_t _nsh; // number of associated straw hits
      // mc diag
      XYZVec _mcpos; // average MC hit position
      Float_t _mctime, _mcdist;
      Int_t _mcpdg, _mcproc, _mcgen; 
// per-hit diagnostics
      vector<ComboHitInfo> _chinfo;
      vector<ComboHitInfoMC> _chinfomc;
      // helper functions
      bool findData(const art::Event& evt);
  };

  ComboHitDiag::ComboHitDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag		(pset.get<int>("diagLevel",1)),
    _mcdiag		(pset.get<bool>("MCdiag",true)),
    _shtag		(pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shptag		(pset.get<art::InputTag>("StrawHitPositionCollection","MakeStrawHitPositions")),
    _chtag		(pset.get<art::InputTag>("ComboHitCollection","CombineHits")),
    _mcdigistag		(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))

  {}

  ComboHitDiag::~ComboHitDiag(){}

  void ComboHitDiag::beginJob() {
    // create diagnostics if requested
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      // detailed diagnostics
      _chdiag=tfs->make<TTree>("chdiag","combo hit diagnostics");
      _chdiag->Branch("pos",&_pos);
      _chdiag->Branch("wdir",&_wdir);
      _chdiag->Branch("wdist",&_wdist,"wdist/F");
      _chdiag->Branch("wres",&_wres,"wres/F");
      _chdiag->Branch("tres",&_tres,"tres/F");
      _chdiag->Branch("time",&_time,"time/F");
      _chdiag->Branch("edep",&_edep,"edep/F");
      _chdiag->Branch("qual",&_qual,"qual/F");
      _chdiag->Branch("nsh",&_nsh,"nsh/I");
      if(_diag > 1)
	_chdiag->Branch("chinfo",&_chinfo);
      if(_mcdiag){
	_chdiag->Branch("mcpos",&_mcpos);
	_chdiag->Branch("mctime",&_mctime,"mctime/F");
	_chdiag->Branch("mcdist",&_mcdist,"mcdist/F");
	_chdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
	_chdiag->Branch("mcproc",&_mcproc,"mcproc/I");
	_chdiag->Branch("mcgen",&_mcgen,"mcgen/I");
	if(_diag > 1)
	  _chdiag->Branch("chinfomc",&_chinfomc);
      }
    }
  }

  void ComboHitDiag::analyze(const art::Event& evt ) {
    // find data in event
    findData(evt);
    // loop over combo hits
    for(auto ch : *_chcol) {
      // general information
      _nsh = ch._nsh;
      _pos = ch._pos;
      _wdir = ch._wdir;
      _wdist = ch._wdist;
      _wres = ch._wres;
      _tres = ch._tres;
      _time = ch._time;
      _edep = ch._edep;
      _qual = ch._qual;
      // center of this wire
      XYZVec cpos = _pos - _wdist*_wdir;
      // now hit-by-hit info
      _chinfo.clear();
      if(_diag > 1){
	for(size_t ish=0;ish < ch._nsh; ++ish) {
	  ComboHitInfo chi;
	  StrawHit const& sh = _shcol->at(ch._sh[ish]);
	  StrawHitPosition const& shp = _shpcol->at(ch._sh[ish]);
	  XYZVec chpos(shp.pos().x(), shp.pos().y(), shp.pos().z());
	  chi._dwire = (chpos-ch._pos).Dot(ch._wdir);
	  chi._dwerr = shp.posRes(StrawHitPosition::wire);
	  chi._dperp = sqrt((chpos-ch._pos).mag2() - chi._dwire*chi._dwire);
	  chi._dtime = std::min(sh.time(TrkTypes::cal),sh.time(TrkTypes::hv)) - ch._time; 
	  chi._dedep = sh.energyDep() - ch._edep;
	  _chinfo.push_back(chi);
	}
      }
      if(_mcdiag){
	_chinfomc.clear();
	// use the 1st hit to define the MC match; this is arbitrary but adequate
	StrawDigiMC const& mcd1 = _mcdigis->at(ch._sh[0]);
	art::Ptr<StepPointMC> const& spmcp = mcd1.stepPointMC(TrkTypes::cal);
	art::Ptr<SimParticle> spp = spmcp->simParticle();
	_mctime = _toff.timeWithOffsetsApplied(*spmcp);
	_mcpdg = spp->pdgId();
	_mcproc = spp->creationCode();
	if(spp->genParticle().isNonnull())
	  _mcgen = spp->genParticle()->generatorId().id();
	// find the relation with each hit
	_mcpos = XYZVec(0.0,0.0,0.0);
	for(size_t ish=0;ish < ch._nsh; ++ish) {
	  ComboHitInfoMC chimc;
	  StrawDigiMC const& mcd = _mcdigis->at(ch._sh[ish]);
	  chimc._rel = MCRelationship::relationship(mcd,mcd1);
	  _chinfomc.push_back(chimc);
	  // find average MC properties
	  _mcpos += XYZVec(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
	}
	_mcpos /= ch._nsh;
	_mcdist = (_mcpos - cpos).Dot(_wdir);

      }
      _chdiag->Fill();
    }
  }

  bool ComboHitDiag::findData(const art::Event& evt){
    _shcol = 0; _shpcol = 0; _chcol = 0;  _mcdigis = 0;
    auto shH = evt.getValidHandle<StrawHitCollection>(_shtag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shptag);
    _shpcol = shpH.product();
    auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
    _chcol = chH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
      _mcdigis = mcdH.product();
      // update time offsets
      _toff.updateMap(evt);
    }
    return _shcol != 0 && _chcol != 0 && (_mcdigis != 0 || !_mcdiag);
  }
}
// Part of the magic that makes this class a module.
using mu2e::ComboHitDiag;
DEFINE_ART_MODULE(ComboHitDiag);

