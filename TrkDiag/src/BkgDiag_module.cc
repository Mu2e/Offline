//
// Low-energy electron background diagnostics.  Split out of FlagBkgHits
//
// Original author D. Brown
//
//
/// framework
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
#include "TTree.h"
#include "Math/VectorUtil.h"
// data
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgQual.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "TrkDiag/inc/BkgHitInfo.hh"
/// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// diagnostics
#include "TrkDiag/inc/StrawHitInfo.hh"
using namespace std;
using CLHEP::Hep3Vector;
using namespace ROOT::Math::VectorUtil;
namespace mu2e 
{

  class BkgDiag : public art::EDAnalyzer {
    public:
      explicit BkgDiag(fhicl::ParameterSet const&);
      virtual ~BkgDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);

    private:
      // helper functions
      void fillStrawHitInfo(size_t ich, StrawHitInfo& bkghinfo) const;
      void fillStrawHitInfoMC(StrawDigiMC const& mcdigi, art::Ptr<SimParticle>const& pptr, StrawHitInfo& shinfo) const;
      bool findData(const art::Event& e);
      void findPrimary(std::vector<StrawDigiIndex>const& dids, art::Ptr<SimParticle>& pptr,double& pmom) const;

      // control flags
      int _diag,_debug;
      bool _mcdiag, _useflagcol;
      float _maxdt, _maxdrho;
  // data tags
      art::InputTag _chTag;
      art::InputTag _shfTag;
      art::InputTag _bkgcTag;
      art::InputTag _bkgqTag;
      art::InputTag _mcdigisTag;
      // time offset
      SimParticleTimeOffset _toff;
      // cache of event objects
      const ComboHitCollection* _chcol;
      const StrawHitFlagCollection* _shfcol;
      const StrawDigiMCCollection *_mcdigis;
      const BkgClusterCollection *_bkgccol;
      const BkgQualCollection *_bkgqcol;

      // background diagnostics
      TTree* _bdiag;
      Int_t _iev;
      CLHEP::Hep3Vector _cpos;
      Float_t _ctime;
      Float_t _mindt, _mindrho;
      Bool_t _isbkg, _isref, _isolated, _stereo;
      Int_t _nactive, _nchits, _nshits, _nstereo, _nsactive, _nbkg;
// BkgQual vars
      float _bkgqualvars[BkgQual::n_vars];
      Int_t _mvastat;
      Float_t _mvaout;

// MC truth variables
      Int_t _ppid, _ppdg, _pgen, _pproc;
      Float_t _pmom;
      Int_t _nconv, _ndelta, _ncompt, _ngconv, _nebkg, _nprot, _nprimary;
      std::vector<BkgHitInfo> _bkghinfo;
  };
 

  BkgDiag::BkgDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("diagLevel",1)),
    _debug(pset.get<int>("debugLevel",0)),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _useflagcol(pset.get<bool>("UseFlagCollection")),
    _maxdt(pset.get<double>("MaxTimeDifference",50.0)), // Maximum time difference (nsec)
    _maxdrho(pset.get<double>("MaxRhoDifference",50.0)), // Maximum transverse distance difference (mm)
    _chTag(pset.get<string>("ComboHitCollection")),
    _shfTag(pset.get<string>("StrawHitFlagCollection","FlagBkgHits")),
    _bkgcTag(pset.get<string>("BackgroundClusterCollection","FlagBkgHits")),
    _bkgqTag(pset.get<string>("BackgroundQualCollection","FlagBkgHits")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {}

  BkgDiag::~BkgDiag(){}

  void BkgDiag::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    if(_diag > 0){
      // detailed delta diagnostics
      _bdiag=tfs->make<TTree>("bkgdiag","background diagnostics");
      // general branches
      _bdiag->Branch("iev",&_iev,"iev/I");
      // cluster info branches
      _bdiag->Branch("cpos",&_cpos,"dx/D:dy/D:dz/D");
      _bdiag->Branch("ctime",&_ctime,"ctime/F");
      _bdiag->Branch("isbkg",&_isbkg,"isbkg/B");
      _bdiag->Branch("isref",&_isref,"isref/B");
      _bdiag->Branch("isolated",&_isolated,"isolated/B");
      _bdiag->Branch("stereo",&_stereo,"stereo/B");
      _bdiag->Branch("mindt",&_mindt,"mindt/F");
      _bdiag->Branch("mindrho",&_mindrho,"mindrho/F");
      _bdiag->Branch("nchits",&_nchits,"nchits/I");
      _bdiag->Branch("nshits",&_nshits,"nshits/I");
      _bdiag->Branch("nactive",&_nactive,"nactive/I");
      _bdiag->Branch("nstereo",&_nstereo,"nstereo/I");
      _bdiag->Branch("nsactive",&_nsactive,"nsactive/I");
      _bdiag->Branch("nbkg",&_nbkg,"nbkg/I");
      // cluster hit info branch
      if(_diag > 1)
	_bdiag->Branch("bkghinfo",&_bkghinfo);
      // Bkg qual info
      for(int ivar=0;ivar < BkgQual::n_vars; ++ivar){
	string vname = BkgQual::varName(static_cast<BkgQual::MVA_varindex>(ivar));
	string bname = vname+string("/F");
	_bdiag->Branch(vname.c_str(),&_bkgqualvars[ivar],bname.c_str());
      }
      _bdiag->Branch("mvaout", &_mvaout,"mvaout/F");
      _bdiag->Branch("mvastat", &_mvastat,"mvastat/I");
      // mc truth branches
      if(_mcdiag){
	_bdiag->Branch("pmom",&_pmom,"pmom/F");
	_bdiag->Branch("ppid",&_ppid,"ppid/I");
	_bdiag->Branch("ppdg",&_ppdg,"ppdg/I");
	_bdiag->Branch("pgen",&_pgen,"pgen/I");
	_bdiag->Branch("pproc",&_pproc,"pproc/I");
	_bdiag->Branch("nprimary",&_nprimary,"nprimary/I");
	_bdiag->Branch("nconv",&_nconv,"nconv/I");
	_bdiag->Branch("ndelta",&_ndelta,"ndelta/I");
	_bdiag->Branch("ncompt",&_ncompt,"ncompt/I");
	_bdiag->Branch("ngconv",&_ngconv,"ngconv/I");
	_bdiag->Branch("nebkg",&_nebkg,"nebkg/I");
	_bdiag->Branch("nprot",&_nprot,"nprot/I");
      }
    }
  }

  void BkgDiag::analyze(const art::Event& event ) {
    _iev = event.event();
    if(!findData(event))
      throw cet::exception("RECO")<<"mu2e::BkgDiag: data missing or incomplete"<< endl;
    // check consistency
    if(_bkgccol->size() != _bkgqcol->size())
      throw cet::exception("RECO")<<"mu2e::BkgDiag: data inconsistent"<< endl;
    // loop over background clusters
    for(size_t ibkg=0;ibkg<_bkgccol->size();++ibkg){
      BkgCluster const& cluster = _bkgccol->at(ibkg);
      BkgQual const& qual = _bkgqcol->at(ibkg);
      // fill cluster info
      _cpos = Geom::Hep3Vec(cluster.pos()); 
      _ctime = cluster.time();
      _isbkg = cluster.flag().hasAllProperties(BkgClusterFlag::bkg);
      _isref = cluster.flag().hasAllProperties(BkgClusterFlag::refined);
      _isolated = cluster.flag().hasAllProperties(BkgClusterFlag::iso);
      _stereo = cluster.flag().hasAllProperties(BkgClusterFlag::stereo);
      // fill Bkg qual info
      for(int ivar=0;ivar < BkgQual::n_vars; ++ivar){
	_bkgqualvars[ivar] = qual[static_cast<BkgQual::MVA_varindex>(ivar)];
      }
      _mvaout = qual.MVAOutput();
      _mvastat = qual.status();
      // info on nearest cluster
      _mindt = _mindrho = 1.0e3;
      for(size_t jbkg = 0; jbkg < _bkgccol->size(); ++jbkg){
	if(ibkg != jbkg){
	  BkgCluster const& ocluster = _bkgccol->at(jbkg);
	  double dt = fabs(ocluster.time() - cluster.time());
	  double drho = sqrt((ocluster.pos()-cluster.pos()).Perp2());
	  // only look at differences whtn the other dimension difference is small
	  if(drho < _maxdrho && dt < _mindt) _mindt = dt;
	  if(dt < _maxdt && drho < _mindrho) _mindrho = drho;
	}
      }
      // fill mc info
      art::Ptr<SimParticle> pptr;
      // loop over hits in this cluster and classify them
      _nconv = 0;
      _nprot = 0;
      _ndelta= 0;
      _ncompt = 0;
      _ngconv = 0;
      _nebkg = 0;
      _nprimary = 0;
      _pmom = 0.0;
      _ppid = _ppdg = _pgen = _pproc = 0;
      if(_mcdiag){
      // fill vector of indices to all digis used in this cluster's hits 
      // this goes recursively through the ComboHit chain
	std::vector<StrawDigiIndex> cdids;
	for(auto const& chit : cluster.hits()){
	  size_t ich = chit.index();
	  // get the list of StrawHit indices associated with this ComboHit
	  _chcol->fillStrawDigiIndices(event,ich,cdids);
	}	
	double pmom(0.0);
	findPrimary(cdids,pptr,pmom);
	_pmom = pmom;
	if(pptr.isNonnull()){
	  _ppid = pptr->id().asInt();
	  _ppdg = pptr->pdgId();
	  _pproc = pptr->creationCode();
	  if( pptr->genParticle().isNonnull())
	    _pgen = pptr->genParticle()->generatorId().id();
	}
      }
      // fill cluster hit info
      _bkghinfo.clear();
      _bkghinfo.reserve(cluster.hits().size());
      _nchits = cluster.hits().size();
      _nshits = 0;
      _nactive = _nstereo = _nsactive = _nbkg = 0;
      bool pce = _pgen==2; // primary from a CE
      for(auto const& chit : cluster.hits()){
	size_t ich = chit.index();
	ComboHit const& ch = _chcol->at(ich);
	_nshits += ch.nStrawHits();
	StrawHitFlag const& shf = chit.flag();
	if(shf.hasAllProperties(StrawHitFlag::active)){
	  _nactive += ch.nStrawHits();
	  if(shf.hasAllProperties(StrawHitFlag::stereo))_nsactive+= ch.nStrawHits();
	}
	if(shf.hasAllProperties(StrawHitFlag::stereo))_nstereo+= ch.nStrawHits();
	if(shf.hasAllProperties(StrawHitFlag::bkg))_nbkg+= ch.nStrawHits();
	// fill hit-specific information
	BkgHitInfo bkghinfo;
	// fill basic straw hit info
	fillStrawHitInfo(chit.index(),bkghinfo);
	if(_mcdiag){
	  std::vector<StrawDigiIndex> dids;
	  _chcol->fillStrawDigiIndices(event,ich,dids);
	  StrawDigiMC const& mcdigi = _mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
	  fillStrawHitInfoMC(mcdigi,pptr,bkghinfo);
	}
	// background hit specific information
	bkghinfo._active = shf.hasAllProperties(StrawHitFlag::active);
	bkghinfo._cbkg = shf.hasAllProperties(StrawHitFlag::bkg);
	bkghinfo._gdist = chit.distance();
	bkghinfo._index = ich;
	// calculate separation to cluster
	Hep3Vector psep = Geom::Hep3Vec(PerpVector(ch.pos()-cluster.pos(),Geom::ZDir()));
	double rho = psep.mag();
	Hep3Vector pdir = psep.unit();
	bkghinfo._rpos = psep;
	bkghinfo._rrho = rho;
	bkghinfo._rerr = std::max(2.5,ch.posRes(ComboHit::wire)*fabs(pdir.dot(ch.wdirCLHEP())));
	//global counting for the cluster: count signal hits only, but background from background is OK
	if(pce){
	  if(bkghinfo._relation==0) _nprimary += ch.nStrawHits(); // couunt only true primary
	} else {
	  if(bkghinfo._relation>=0 && bkghinfo._relation <=3) _nprimary += ch.nStrawHits(); // count primar + mother/daughter/sibling
	}
	if(bkghinfo._mcgen == 2)_nconv += ch.nStrawHits();
	if(abs(bkghinfo._mcpdg) == PDGCode::e_minus && bkghinfo._mcgen <0){
	  _nebkg += ch.nStrawHits();
	  if(bkghinfo._mcproc == ProcessCode::eIoni ||bkghinfo._mcproc == ProcessCode::hIoni ){
	    _ndelta += ch.nStrawHits();
	  } else if(bkghinfo._mcproc == ProcessCode::compt){
	    _ncompt += ch.nStrawHits();
	  } else if(bkghinfo._mcproc == ProcessCode::conv){
	    _ngconv += ch.nStrawHits();
	  }
	}
	if(bkghinfo._mcpdg == 2212)_nprot+= ch.nStrawHits();
	_bkghinfo.push_back(bkghinfo);
      }
      _bdiag->Fill();
    }
  }

  bool BkgDiag::findData(const art::Event& evt){
    _chcol = 0; _shfcol = 0; _bkgccol = 0; _bkgqcol = 0; _mcdigis = 0;
    // nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
    _chcol = chH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto bkgcH = evt.getValidHandle<BkgClusterCollection>(_bkgcTag);
    _bkgccol = bkgcH.product();
    auto bkgqH = evt.getValidHandle<BkgQualCollection>(_bkgqTag);
    _bkgqcol = bkgqH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      // update time offsets
      _toff.updateMap(evt);
    }
    return _chcol != 0 && _shfcol != 0 && _bkgccol != 0 && _bkgqcol != 0
      && (_mcdigis != 0  || !_mcdiag);
  }


  void BkgDiag::findPrimary(std::vector<uint16_t>const& dids, art::Ptr<SimParticle>& pptr,double& pmom) const { 
    static StrawEnd itdc(TrkTypes::cal);
    // find the unique simparticles which produced these hits
    std::set<art::Ptr<SimParticle> > pp;
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      if(spp.isNonnull()){
	pp.insert(spp);
      }
    }
    // map these particles back to each other, to compress out particles generated inside the cluster 
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
	    MCRelationship::relation rel = MCRelationship::relationship(sppi,sppj);
	    if(rel==MCRelationship::daughter || rel == MCRelationship::udaughter){
	      spmap[sppi] = sppj;
	      break;
	    } else if(rel == MCRelationship::mother || rel == MCRelationship::umother){
	      spmap[sppj] = sppi;
	    } else if(rel == MCRelationship::sibling || rel == MCRelationship::usibling){
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
	if( !(ifnd->second == ifnd->first)){
	  changed = true;
	  spmap[im->first] = ifnd->second;
	}
      }
    }
    // find the most likely ultimate parent for this cluster.  Also fill general info
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
    // find the momentum for the first step point from the primary particle in this delta
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      // use TDC channel 0 to define the MC match
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();

      if(spp == pptr){
	pmom = spmcp->momentum().mag();
	break;
      }
    }
  }

  void BkgDiag::fillStrawHitInfoMC(StrawDigiMC const& mcdigi, art::Ptr<SimParticle>const& pptr, StrawHitInfo& shinfo) const {
    // use TDC channel 0 to define the MC match
    StrawEnd itdc(TrkTypes::cal);
    art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
    art::Ptr<SimParticle> const& spp = spmcp->simParticle();
    shinfo._mct0 = _toff.timeWithOffsetsApplied(*spmcp);
    shinfo._mcht = mcdigi.wireEndTime(itdc);
    shinfo._mcpdg = spp->pdgId();
    shinfo._mcproc = spp->creationCode();
    shinfo._mcedep = mcdigi.energySum();
    shinfo._mcgen = -1;
    if(spp->genParticle().isNonnull())
      shinfo._mcgen = spp->genParticle()->generatorId().id();

    shinfo._mcpos = spmcp->position();
    shinfo._mctime = shinfo._mct0;
    shinfo._mcedep = mcdigi.energySum();;
    shinfo._mcmom = spmcp->momentum().mag();
    double cosd = spmcp->momentum().cosTheta();
    shinfo._mctd = cosd/sqrt(1.0-cosd*cosd);
    // relationship to parent
    shinfo._relation=MCRelationship::none;
    if(spmcp.isNonnull() && pptr.isNonnull()){
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      if(spp.isNonnull()){
	shinfo._relation = MCRelationship::relationship(spp,pptr);
      }
    }
  }

  void BkgDiag::fillStrawHitInfo(size_t ich, StrawHitInfo& shinfo) const {
    ComboHit const& ch = _chcol->at(ich);
    StrawHitFlag shf;
    if(_useflagcol)
      shf = _shfcol->at(ich);
    else
      shf = ch.flag();

    shinfo._stereo = shf.hasAllProperties(StrawHitFlag::stereo);
    shinfo._tdiv = shf.hasAllProperties(StrawHitFlag::tdiv);
    shinfo._esel = shf.hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel = shf.hasAllProperties(StrawHitFlag::radsel);
    shinfo._tsel = shf.hasAllProperties(StrawHitFlag::timesel);
    shinfo._strawxtalk = shf.hasAllProperties(StrawHitFlag::strawxtalk);
    shinfo._elecxtalk = shf.hasAllProperties(StrawHitFlag::elecxtalk);
    shinfo._isolated = shf.hasAllProperties(StrawHitFlag::isolated);
    shinfo._bkg = shf.hasAllProperties(StrawHitFlag::bkg);

    shinfo._pos = ch.posCLHEP();
    shinfo._time = ch.time();
    shinfo._rho = ch.posCLHEP().perp();
    shinfo._wres = ch.posRes(ComboHit::wire);
    shinfo._tres = ch.posRes(ComboHit::trans);
    // info depending on stereo hits
      shinfo._chisq = ch.qual();
    shinfo._edep = ch.energyDep();
    StrawId const& sid = ch.sid();
    shinfo._plane = sid.plane();
    shinfo._panel = sid.panel();
    shinfo._layer = sid.layer();
    shinfo._straw = sid.straw();
    shinfo._stereo = ch.flag().hasAllProperties(StrawHitFlag::stereo);
    shinfo._tdiv = ch.flag().hasAllProperties(StrawHitFlag::tdiv);

  }
} // mu2e namespace

// Part of the magic that makes this class a module.
using mu2e::BkgDiag;
DEFINE_ART_MODULE(BkgDiag);

