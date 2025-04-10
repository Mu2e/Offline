//
// Low-energy electron background diagnostics.  Split out of FlagBkgHits
//
// Original author D. Brown
//
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/TrkDiag/inc/BkgHitInfo.hh"

#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
#include "Math/VectorUtil.h"

using namespace ROOT::Math::VectorUtil;



namespace mu2e
{

  class BkgDiag2 : public art::EDAnalyzer {
    public:

      struct Config {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>    ComboHitCollection       {Name("ComboHitCollection"),        Comment("ComboHit collection name")             };
        fhicl::Atom<art::InputTag>    BkgClusterCollection     {Name("BkgClusterCollection"),      Comment("BackgroundCluster collection name")    };
        fhicl::Atom<art::InputTag>    BkgClusterHitCollection  {Name("BkgClusterHitCollection"),   Comment("BackgroundClusterHit collection name") };
        fhicl::Atom<art::InputTag>    StrawDigiMCCollection    {Name("StrawDigiMCCollection"),     Comment("StrawDigiMC collection name")          };
        fhicl::Atom<art::InputTag>    MCPrimary                {Name("MCPrimary"),                 Comment("MC Primary Particle")                  };
        fhicl::Sequence<std::string>  bkgmsk                   {Name("BackgroundMask"),            Comment("Bkg hit selection mask")               };
        fhicl::Sequence<std::string>  sigmsk                   {Name("SignalMask"),                Comment("Signal hit selection mask")            };
        fhicl::Atom<float>            coreR                    {Name("CoreRadius"),                Comment("Core radius in cluster")               };
        fhicl::Atom<bool>             mapSimAncestor           {Name("MapSimAncestor"),            Comment("Remap sims to ancestors")              };
        fhicl::Atom<bool>             mcdiag                   {Name("MCDiag"),                    Comment("Monte Carlo Diag")                     };
        fhicl::Atom<bool>             hdiag                    {Name("HitDiag"),                   Comment("Hit-level Diag tuple")                 };
        fhicl::Atom<bool>             bkdiag                   {Name("BkgHitDiag"),                Comment("Bgk hit Diag ntuple")                  };
        fhicl::Atom<int>              debug                    {Name("debugLevel"),                Comment("Debug Level"),0                        };
     };

      explicit     BkgDiag2(const art::EDAnalyzer::Table<Config>& config);
      virtual     ~BkgDiag2();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);


    private:
      void fillStrawHitInfo(size_t ich, StrawHitInfo& bkghinfo) const;
      void fillStrawHitInfoMC(const StrawDigiMC& mcdigi, art::Ptr<SimParticle>const& mptr, StrawHitInfo& shinfo) const;
      bool findData(const art::Event& e);
      void findSimAncestor(const BkgCluster& cluster, art::Ptr<SimParticle>& mptr, StrawHitIndex& shid, std::vector<int>& icontrib) const;


      art::ProductToken<ComboHitCollection>      _chToken;
      art::ProductToken<BkgClusterCollection>    _bkgcToken;
      art::ProductToken<BkgClusterHitCollection> _bkghToken;
      art::ProductToken<StrawDigiMCCollection>   _mcdigisToken;
      art::ProductToken<PrimaryParticle>         _mcprimaryToken;
      StrawHitFlag                               _bkgmask;
      StrawHitFlag                               _sigmask;
      float                                      _coreR2;
      bool                                       _mapSimAncestor;
      bool                                       _bkdiag;
      bool                                       _hdiag;
      bool                                       _mcdiag;
      int                                        _debugLevel;

      // cache of event objects
      const ComboHitCollection      *_chcol;
      const StrawDigiMCCollection   *_mcdigis;
      const PrimaryParticle         *_mcprimary;
      const BkgClusterCollection    *_bkgccol;
      const BkgClusterHitCollection *_bkghitcol;

      // Ntuple - background diagnostics
      TTree *_bcdiag, *_bhdiag;
      float _rmscposx,_rmscposy,_rmscrho,_rmsctime,_cchi2,_ccons,_ctime,_cpitch,_cqual,_csthqual,_surface;
      float _ccomqual,_cecc,_cdeltaR,_cdeltaP,_avecedep,_mindt,_mindrho,_crho,_zmin,_zmax,_zgap,_pfrac,_kQ;
      bool  _isinit,_isbkg,_isref,_isolated,_stereo;
      int   _iev,_cluIdx, _nactive, _nch, _nsh, _nsth, _nsha, _nbkg,_cndof,_np, _fp, _lp, _pgap, _nCoreMax;
      XYZVectorF _cpos;

      // Ntuple - MC truth variables
      int _ndelta, _ncompt, _ngconv, _nprot, _nmain, _nsibling, _nrel;
      int _mpdg, _mproc, _ncontrib, _icontrib[512], _prel;
      float _mmomMag;
      XYZVectorF _mmom, _mpos;
      std::vector<BkgHitInfo> _bkghinfo;
      std::vector<int>   _chindices;

      //Ntuple -  All hits variables
      Int_t              _nhits;
      std::vector<int>   _hitPdg,_hitproc,_hitnch,_hitnsh,_hitid,_hitflag;
      std::vector<float> _hitPosX,_hitPosY,_hitPosZ,_hitTime,_hiteDep,_hitDirX;
      std::vector<float> _hitDirY,_hitDirZ,_hitDirQ,_hitMCmomX,_hitMCmomY,_hitMCmomZ;
  };


  BkgDiag2::BkgDiag2(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _chToken        { consumes<ComboHitCollection>(config().ComboHitCollection() ) },
    _bkgcToken      { consumes<BkgClusterCollection>(config().BkgClusterCollection() ) },
    _bkghToken      { consumes<BkgClusterHitCollection>(config().BkgClusterHitCollection() ) },
    _mcdigisToken   { consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection() ) },
    _mcprimaryToken { consumes<PrimaryParticle>(config().MCPrimary() ) },
    _bkgmask        { config().bkgmsk() },
    _sigmask        { config().sigmsk() },
    _coreR2         { config().coreR()*config().coreR()},
    _mapSimAncestor { config().mapSimAncestor() },
    _bkdiag         { config().bkdiag() },
    _hdiag          { config().hdiag()  },
    _mcdiag         { config().mcdiag() },
    _debugLevel     { config().debug()  }
  {}

  BkgDiag2::~BkgDiag2(){}

  void BkgDiag2::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    _iev=0;
    _bcdiag=tfs->make<TTree>("bkgcdiag","background cluster diagnostics");
    _bcdiag->Branch("iev",      &_iev,"iev/I");
    _bcdiag->Branch("cpos",     &_cpos);
    _bcdiag->Branch("rmscposx", &_rmscposx, "rmscposx/F");
    _bcdiag->Branch("rmscposy", &_rmscposy, "rmscposy/F");
    _bcdiag->Branch("rmscrho",  &_rmscrho,  "rmscrho/F");
    _bcdiag->Branch("cpitch",   &_cpitch,   "cpitch/F");
    _bcdiag->Branch("cqual",    &_cqual,    "cqual/F");
    _bcdiag->Branch("csthqual", &_csthqual, "csthqual/F");
    _bcdiag->Branch("ccomqual", &_ccomqual, "ccomqual/F");
    _bcdiag->Branch("cecc",     &_cecc,     "cecc/F");
    _bcdiag->Branch("cchi2",    &_cchi2,    "cchi2/F");
    _bcdiag->Branch("ccons",    &_ccons,    "ccons/F");
    _bcdiag->Branch("ctime",    &_ctime,    "ctime/F");
    _bcdiag->Branch("cdeltaR",  &_cdeltaR,  "cdeltaR/F");
    _bcdiag->Branch("cdeltaP",  &_cdeltaP,  "cdeltaP/F");
    _bcdiag->Branch("rmsctime", &_rmsctime, "rmsctime/F");
    _bcdiag->Branch("avecedep", &_avecedep, "avecedep/F");
    _bcdiag->Branch("isinit",   &_isinit,   "isinit/B");
    _bcdiag->Branch("isbkg",    &_isbkg,    "isbkg/B");
    _bcdiag->Branch("isref",    &_isref,    "isref/B");
    _bcdiag->Branch("isolated", &_isolated, "isolated/B");
    _bcdiag->Branch("stereo",   &_stereo,   "stereo/B");
    _bcdiag->Branch("mindt",    &_mindt,    "mindt/F");
    _bcdiag->Branch("mindrho",  &_mindrho,  "mindrho/F");
    _bcdiag->Branch("csurface", &_surface,  "csurface/F");
    _bcdiag->Branch("ncoreMax", &_nCoreMax, "ncoreMax/I");
    _bcdiag->Branch("nch",      &_nch,      "nch/I");
    _bcdiag->Branch("nsh",      &_nsh,      "nsh/I");
    _bcdiag->Branch("nsth",     &_nsth,     "nsth/I");
    _bcdiag->Branch("nactive",  &_nactive,  "nactive/I");
    _bcdiag->Branch("nsha",     &_nsha,     "nsha/I");
    _bcdiag->Branch("nbkg",     &_nbkg,     "nbkg/I");
    _bcdiag->Branch("cndof",    &_cndof,    "cndof/I");
    _bcdiag->Branch("cluIdx",   &_cluIdx,   "cluIdx/I");
    _bcdiag->Branch("crho",     &_crho,     "crho/F");
    _bcdiag->Branch("zmin",     &_zmin,     "zmin/F");
    _bcdiag->Branch("zmax",     &_zmax,     "zmax/F");
    _bcdiag->Branch("zgap",     &_zgap,     "zgap/F");
    _bcdiag->Branch("np",       &_np,       "np/I");
    _bcdiag->Branch("pfrac",    &_pfrac,    "pfrac/F");
    _bcdiag->Branch("kQ",       &_kQ,       "kQ/F");
    _bcdiag->Branch("fp",       &_fp,       "fp/I");
    _bcdiag->Branch("lp",       &_lp,       "lp/I");
    _bcdiag->Branch("pgap",     &_pgap,     "pgap/I");
    _bcdiag->Branch("nhits",    &_nhits,    "nhits/I");
    _bcdiag->Branch("chindices",&_chindices);

    // cluster hit info branch
    if( _bkdiag) _bcdiag->Branch("bkghinfo",&_bkghinfo);

    // mc truth branches
    if (_mcdiag) {
      _bcdiag->Branch("mmom",    &_mmom);
      _bcdiag->Branch("mopos",   &_mpos);
      _bcdiag->Branch("mmomMag", &_mmomMag, "mmomMag/F");
      _bcdiag->Branch("mpdg",    &_mpdg,    "mpdg/I");
      _bcdiag->Branch("mproc",   &_mproc,   "mproc/I");
      _bcdiag->Branch("prel",    &_prel,    "prel/I");
      _bcdiag->Branch("nmain",   &_nmain,   "nmain/I");
      _bcdiag->Branch("nsibling",&_nsibling,"nsibling/I");
      _bcdiag->Branch("nrel",    &_nrel,    "nrel/I");
      _bcdiag->Branch("ndelta",  &_ndelta,  "ndelta/I");
      _bcdiag->Branch("ncompt",  &_ncompt,  "ncompt/I");
      _bcdiag->Branch("ngconv",  &_ngconv,  "ngconv/I");
      _bcdiag->Branch("nprot",   &_nprot,   "nprot/I");
      _bcdiag->Branch("ncontrib",&_ncontrib,"ncontrib/I");
      _bcdiag->Branch("icontrib",&_icontrib,"icontrib[ncontrib]/I");
    }

    if (_hdiag){
      _bhdiag = tfs->make<TTree>("bkghdiag","background hit diagnostics");
      _bhdiag->Branch("iev",     &_iev,       "iev/I");
      _bhdiag->Branch("nhits",   &_nhits,     "nhits/I");
      _bhdiag->Branch("posX",    &_hitPosX    );
      _bhdiag->Branch("posY",    &_hitPosY    );
      _bhdiag->Branch("posZ",    &_hitPosZ    );
      _bhdiag->Branch("dirX",    &_hitDirX    );
      _bhdiag->Branch("dirY",    &_hitDirY    );
      _bhdiag->Branch("dirZ",    &_hitDirZ    );
      _bhdiag->Branch("dirQ",    &_hitDirQ    );
      _bhdiag->Branch("time",    &_hitTime    );
      _bhdiag->Branch("eDep",    &_hiteDep    );
      _bhdiag->Branch("nch",     &_hitnch     );
      _bhdiag->Branch("nsh",     &_hitnsh     );
      _bhdiag->Branch("flag",    &_hitflag    );
      if (_mcdiag){
        _bhdiag->Branch("mcpdg",  &_hitPdg     );
        _bhdiag->Branch("mcproc", &_hitproc    );
        _bhdiag->Branch("mcid",   &_hitid      );
        _bhdiag->Branch("mcmomX", &_hitMCmomX  );
        _bhdiag->Branch("mcmomY", &_hitMCmomY  );
        _bhdiag->Branch("mcmomZ", &_hitMCmomZ  );
      }
    }
  }




  void BkgDiag2::analyze(const art::Event& event ) {

    if (!findData(event))
      throw cet::exception("RECO")<<"mu2e::BkgDiag2: data missing or incomplete"<< std::endl;


    // Ntuple - all hits
    //----------------------------------
    if (_hdiag){
      _nhits = 0;
      _hitPosX.clear(); _hitPosY.clear();_hitPosZ.clear();
      _hitDirX.clear(); _hitDirY.clear();_hitDirZ.clear();  _hitDirQ.clear();
      _hitTime.clear(); _hiteDep.clear();_hitMCmomX.clear();_hitMCmomY.clear();_hitMCmomZ.clear();
      _hitnch.clear();  _hitnsh.clear(); _hitPdg.clear();   _hitproc.clear();  _hitid.clear(); _hitflag.clear();

      for (size_t ich=0;ich<_chcol->size();++ich){

        bool passFlag = _chcol->at(ich).flag().hasAllProperties(_sigmask) && !_chcol->at(ich).flag().hasAnyProperty(_bkgmask);
        int flagVal   = passFlag ? 1 : 0;

        _hitPosX.push_back(_chcol->at(ich).pos().x());
        _hitPosY.push_back(_chcol->at(ich).pos().y());
        _hitPosZ.push_back(_chcol->at(ich).pos().z());
        _hitDirX.push_back(_chcol->at(ich).hDir().x());
        _hitDirY.push_back(_chcol->at(ich).hDir().y());
        _hitDirZ.push_back(_chcol->at(ich).hDir().z());
        _hitDirQ.push_back(_chcol->at(ich).qual());
        _hitTime.push_back(_chcol->at(ich).correctedTime());
        _hiteDep.push_back(_chcol->at(ich).energyDep());
        _hitnch.push_back(_chcol->at(ich).nCombo());
        _hitnsh.push_back(_chcol->at(ich).nStrawHits());
        _hitflag.push_back(flagVal);

        art::Ptr<SimParticle> spp; // main SimParticle for this
        if (_mcdiag) {
          std::vector<StrawDigiIndex> dids;
          _chcol->fillStrawDigiIndices(ich,dids);

          auto const& mcdigi = _mcdigis->at(dids[0]);
          auto const& spp    = mcdigi.earlyStrawGasStep()->simParticle();
          _hitPdg.push_back(spp->pdgId());
          _hitproc.push_back(spp->creationCode());
          _hitid.push_back(spp->id().asInt());
          _hitMCmomX.push_back(mcdigi.earlyStrawGasStep()->momentum().x());
          _hitMCmomY.push_back(mcdigi.earlyStrawGasStep()->momentum().y());
          _hitMCmomZ.push_back(mcdigi.earlyStrawGasStep()->momentum().z());
          }
        ++_nhits;
      }

      _bhdiag->Fill();
    }




    // Ntuple - clusters
    //----------------------------------

    _cluIdx=0;
    for (size_t ibkg=0;ibkg < _bkgccol->size();++ibkg) {
      const auto& cluster = _bkgccol->at(ibkg);
      const auto& hits    = cluster.hits();

      //-- MC information
      art::Ptr<SimParticle> mptr;
      _mpos     = XYZVectorF();
      _mmom     = XYZVectorF();
      _mmomMag  = 0;
      _mpdg     = 0;
      _mproc    = 0;
      _ncontrib = 0;
      _prel     = -1;

      if (_mcdiag){

        std::vector<int> icontrib;
        StrawHitIndex shid(0);
        findSimAncestor(cluster, mptr, shid, icontrib);

        for (int ic : icontrib) {_icontrib[_ncontrib]=ic; ++_ncontrib;}

        if (mptr.isNonnull()) {
          _mpdg    = mptr->pdgId();
          _mproc   = mptr->creationCode();
          _mpos    = mptr->startPosXYZ();

          const auto& mcdigi = _mcdigis->at(shid);
          const auto& sgsp   = mcdigi.earlyStrawGasStep();
          _mmom    = sgsp->momentum();
          _mmomMag = sqrt(_mmom.Mag2());

          for (auto const& mcmptr : _mcprimary->primarySimParticles()) {
            MCRelationship rel(mcmptr,mptr);
            if (rel.relationship() > MCRelationship::none){
              if (_prel > MCRelationship::none) _prel = std::min(_prel,(int)rel.relationship());
              else                              _prel = rel.relationship();
            }
          }
        }
      }


      //-- Global information
      _nch                = cluster.hits().size();
      _crho               = sqrtf(cluster.pos().perp2());
      _cpos               = cluster.pos();
      _ctime              = cluster.time();
      _kQ                 = cluster.getKerasQ();
      _isinit             = cluster.flag().hasAllProperties(BkgClusterFlag::init);
      _isbkg              = cluster.flag().hasAllProperties(BkgClusterFlag::bkg);
      _isref              = cluster.flag().hasAllProperties(BkgClusterFlag::refined);
      _isolated           = cluster.flag().hasAllProperties(BkgClusterFlag::iso);
      _stereo             = cluster.flag().hasAllProperties(BkgClusterFlag::stereo);

      if (cluster.getDistMethod() == BkgCluster::chi2) {
        _cchi2 = cluster.points().chisquared();
        _ccons = cluster.points().consistency();
        _cndof = cluster.points().nDOF();
      }



      //-- Spatial and temporal dimensions
      float sumEdep(0.), sumEcc(0.), sumwEcc(0.);
      float sqrSumDeltaX(0.), sqrSumDeltaY(0.), sqrSumDeltaR(0.), sqrSumDeltaT(0.), sqrSumQual(0.);
      float SumDeltaX(0.), SumDeltaY(0.), SumDeltaR(0.), SumDeltaT(0.);
      float rhmin(9999), rhmax(0),   phmin(100), phmax(-100), xmin(999), xmax(-999), ymin(999), ymax(-999);
      _nsh = _nactive = _nsha = _nbkg = 0;
      _chindices.clear();

      for (const auto& ich : hits) {
        const auto& ch   = _chcol->at(ich);
        const auto& shf  = _bkghitcol->at(ich).flag();
        auto wecc        = ch.nStrawHits();

        float dx         = ch.pos().X() - _cpos.X();
        float dy         = ch.pos().Y() - _cpos.Y();
        float dt         = ch.time() - _ctime;
        float dr         = sqrtf(ch.pos().perp2()) - sqrt(_cpos.perp2());

        sumEdep         += ch.energyDep()*ch.nStrawHits();
        sqrSumDeltaX    += dx*dx;
        sqrSumDeltaY    += dy*dy;
        sqrSumDeltaR    += dr*dr;
        sqrSumDeltaT    += dt*dt;
        SumDeltaX       += dx;
        SumDeltaY       += dy;
        SumDeltaR       += dr;
        SumDeltaT       += dt;
        sumEcc          += std::sqrt(1-(ch.vVar()/ch.uVar()))*wecc;
        sumwEcc         += wecc;

        float phi = ch.pos().phi() - _chcol->at(cluster.hits().at(0)).pos().phi();
        if (phi > M_PI)  phi -= M_PI;
        if (phi < -M_PI) phi += M_PI;
        phmin = std::min(phi,phmin);
        phmax = std::max(phi,phmax);
        rhmin = std::min(sqrtf(ch.pos().perp2()),rhmin);
        rhmax = std::max(sqrtf(ch.pos().perp2()),rhmax);
        xmin  = std::min(xmin,ch.pos().x());
        xmax  = std::max(xmax,ch.pos().x());
        ymin  = std::min(ymin,ch.pos().y());
        ymax  = std::max(ymax,ch.pos().y());

        _nsh += ch.nStrawHits();
        if (shf.hasAllProperties(StrawHitFlag::bkg)) _nbkg += ch.nStrawHits();

        if (shf.hasAllProperties(StrawHitFlag::active)){
          _nactive += ch.nStrawHits();
          if (StrawIdMask::station == ch._mask) _nsha += ch.nStrawHits();
        }
        _chindices.push_back(ich);
      }
      float varX = std::max(0.0f,sqrSumDeltaX/_nch - SumDeltaX*SumDeltaX/_nch/_nch);
      float varY = std::max(0.0f,sqrSumDeltaY/_nch - SumDeltaX*SumDeltaY/_nch/_nch);
      float varR = std::max(0.0f,sqrSumDeltaR/_nch - SumDeltaX*SumDeltaR/_nch/_nch);
      float varT = std::max(0.0f,sqrSumDeltaT/_nch - SumDeltaX*SumDeltaT/_nch/_nch);

      _avecedep = sumEdep/_nsh;
      _cecc     = sumEcc/sumwEcc;
      _rmscposx = sqrt(varX);
      _rmscposy = sqrt(varY);
      _rmscrho  = sqrt(varR);
      _rmsctime = sqrt(varT);
      _cdeltaR  = (rhmax-rhmin)/_nch;
      _cdeltaP  = std::min(phmax - phmin, 3.14159f - abs(phmax - phmin))/_nch;


      //-- Closet cluser - distance or time
      _mindt = _mindrho = 1.0e4;
      for (size_t jbkg = 0; jbkg < _bkgccol->size(); ++jbkg){
        if (ibkg == jbkg) continue;
        float dt   = fabs(_bkgccol->at(jbkg).time() - cluster.time());
        float drho = sqrt((_bkgccol->at(jbkg).pos() - cluster.pos()).Perp2());
        if (drho <50) _mindt = std::min(dt,_mindt);
        if (dt < 25)  _mindrho = std::min(drho,_mindrho);
      }

      //-- Number of core hits in the cluster
      _nCoreMax=0;
      for (size_t i=0;i<hits.size();++i) {
        const auto& chi = _chcol->at(hits[i]);
        int nCore(chi.nStrawHits());
        for (size_t j=i+1;j<hits.size();++j) {
          const auto& chj = _chcol->at(hits[j]);
          float dx  = chi.pos().x()-chj.pos().x();
          float dy  = chi.pos().y()-chj.pos().y();
          float dr2 = dx*dx+dy*dy;
          if (dr2 < _coreR2) nCore += chj.nStrawHits();
        }
        _nCoreMax = std::max(nCore, _nCoreMax);
      }


      //-- Cluster surface approximation
      _surface=0;
      if (hits.size()>1){
        const int nbinSurface(10);
        std::array<bool,nbinSurface*nbinSurface> occupied{};
        for (const auto& ich : hits) {
          const auto& ch   = _chcol->at(ich);
          int ibinX        = static_cast<int>((ch.pos().x()-xmin)/(xmax-xmin+1e-3)*nbinSurface);
          int ibinY        = static_cast<int>((ch.pos().y()-ymin)/(ymax-ymin+1e-3)*nbinSurface);
          occupied[ibinX*nbinSurface+ibinY] = true;
        }
        int nCells = std::count_if(occupied.begin(),occupied.end(),[](auto i){return i;});
        _surface = nCells*(xmax-xmin)*(ymax-ymin)/nbinSurface/nbinSurface;
      }


      //-- Z gaps and plane fraction
      std::vector<float> hz;
      std::array<bool,StrawId::_nplanes> hp{false};

      for (auto const& ich : cluster.hits()){
        hz.push_back(_chcol->at(ich).pos().Z());
        hp[_chcol->at(ich).strawId().plane()] = true;
      }

      std::sort(hz.begin(),hz.end());
      _zmin = hz.front();
      _zmax = hz.back();
      _zgap = 0.0;
      for (unsigned iz=1;iz<hz.size();++iz) _zgap=std::max(_zgap,hz[iz]-hz[iz-1]);

      _lp   = -1; // last plane in cluster
      _fp   = StrawId::_nplanes; // first plane in cluster
      _np   = 0;  //# of planes
      _pgap = 0;  // largest plane gap
      int lp(-1); // last plane seen
      for (int ip=0; ip < StrawId::_nplanes; ++ip){
        if (!hp[ip]) continue;
        _np++;
        if (lp > 0 && ip - lp -1 > _pgap) _pgap = ip - lp -1;
        if (ip > _lp) _lp = ip;
        if (ip < _fp) _fp = ip;
        lp = ip;
      }
      _pfrac = (_lp-_fp>0) ? static_cast<float>(_np)/static_cast<float>(_lp - _fp) : 0;




      //-- Tracklet information
      _nsth=0;
      float sumPitch(0.), sumwPitch(0.);
      for (auto const& ich : cluster.hits()){
        const auto& ch   = _chcol->at(ich);
        if (!ch.flag().hasAllProperties(StrawHitFlag::sline)) continue;

        ++_nsth;
        sqrSumQual     += std::pow(ch.qual(),2);

        sumPitch       += (ch.hcostVar() > 1e-6) ? ch.hDir().z() / ch.hcostVar() : 0;
        sumwPitch      += (ch.hcostVar() > 1e-6) ? 1.0 / ch.hcostVar() : 0;
        //sumPitch     += ch.hDir().z()*ch.nStrawHits();
        //sumwPitch    +=  ch.nStrawHits();
      }
      _cpitch   = sumwPitch > 1e-6 ? sumPitch/sumwPitch : 0.;
      _cqual    = std::sqrt(sqrSumQual/_nch); //average SLine fit quality of the cluster
      _csthqual = float(_nsth)/float(_nch);   //SLine ratio the cluster
      _ccomqual = _cqual + _csthqual;         //combined quality metric



      //-- Full bkg hit info
      _bkghinfo.clear();
      _bkghinfo.reserve(cluster.hits().size());
      _nrel = _nmain = _nsibling = _nprot =0;


      for (auto const& ich : cluster.hits()){
        const auto& ch   = _chcol->at(ich);
        const auto& bhit = _bkghitcol->at(ich);
        const auto& shf  = bhit.flag();

        BkgHitInfo bkghinfo;
        fillStrawHitInfo(ich,bkghinfo);

        auto psep        = ch.pos()-cluster.pos();
        auto pdir        = PerpVector(psep,GenVector::ZDir()).Unit();
        bkghinfo._rpos   = psep;
        bkghinfo._rerr   = std::max(float(2.5),ch.posRes(ComboHit::wire)*fabs(pdir.Dot(ch.uDir())));

        bkghinfo._active = shf.hasAllProperties(StrawHitFlag::active);
        bkghinfo._cbkg   = shf.hasAllProperties(StrawHitFlag::bkg);
        bkghinfo._gdist  = bhit.distance();
        bkghinfo._index  = ich;

        if (_mcdiag) {
          std::vector<StrawDigiIndex> dids;
          _chcol->fillStrawDigiIndices(ich,dids);
          const StrawDigiMC& mcdigi = _mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
          fillStrawHitInfoMC(mcdigi,mptr,bkghinfo);

          if (bkghinfo._mrel==MCRelationship::same)   _nmain    += ch.nStrawHits();
          if (bkghinfo._mrel>MCRelationship::sibling) _nsibling += ch.nStrawHits();
          if (bkghinfo._mrel>MCRelationship::none)    _nrel     += ch.nStrawHits();
          if (bkghinfo._mcpdg == PDGCode::proton)     _nprot    += ch.nStrawHits();
        }
        _bkghinfo.push_back(bkghinfo);
      }


      //-- Cluster - end cluster info
      _bcdiag->Fill();
      ++_cluIdx;
    }
    ++_iev;
  }


  //------------------------------------------------------------------------------------------------------------------------------
  bool BkgDiag2::findData(const art::Event& evt){
    _chcol = 0; _bkgccol = 0; _mcdigis = 0;
    auto chH = evt.getValidHandle(_chToken);
    _chcol = chH.product();
    auto bkgcH = evt.getValidHandle(_bkgcToken);
    _bkgccol = bkgcH.product();
    auto bkghH = evt.getValidHandle(_bkghToken);
    _bkghitcol = bkghH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle(_mcdigisToken);
      _mcdigis = mcdH.product();
      auto mcpH = evt.getValidHandle(_mcprimaryToken);
      _mcprimary = mcpH.product();
    }
    return _chcol != 0 && _bkgccol != 0 && ( (_mcdigis != 0 && _mcprimary != 0)  || !_mcdiag);
  }


  //------------------------------------------------------------------------------------------------------------------------------
  void BkgDiag2::findSimAncestor(const BkgCluster& cluster, art::Ptr<SimParticle>& mptr, StrawHitIndex& shid, std::vector<int>& icontrib) const
  {
    std::vector<StrawDigiIndex> dids;
    for (const auto& ich : cluster.hits()) _chcol->fillStrawDigiIndices(ich,dids);

    std::set<art::Ptr<SimParticle>> pp;
    for (auto id : dids) {
      const StrawDigiMC& mcdigi = _mcdigis->at(id);
      const art::Ptr<SimParticle>& spp = mcdigi.earlyStrawGasStep()->simParticle();
      if (spp.isNonnull()) pp.insert(spp);
    }

    // map these particles back to each other, to compress out particles generated inside the cluster
    // look for particles produced at the same point, like conversions. It's not enough to look for the same parent,
    // as that parent could produce multiple daughters at different times. Regardless of mechanism or genealogy,
    // call these 'the same' as they will contribute equally to the spiral

    std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle>> spmap;
    for (auto ipp=pp.begin(); ipp!=pp.end(); ++ipp){
      art::Ptr<SimParticle> sppi = *ipp;
      spmap[sppi] = sppi;
    }

    for (auto ipp=pp.begin();ipp!=pp.end();++ipp){
      if (!_mapSimAncestor) break; //exit if we don't remap ancestors (avoids big ugly if block)

      art::Ptr<SimParticle> sppi = *ipp;
      if (!sppi->genParticle().isNull()) continue;

      auto jpp=ipp;++jpp;
      for(;jpp!=pp.end();++jpp){
        art::Ptr<SimParticle> sppj = *jpp;
        if (!sppj->genParticle().isNull())continue;

        // call the particles 'the same' if they are related and were produced near each other
        MCRelationship rel(sppi,sppj);
        if (rel==MCRelationship::daughter || rel == MCRelationship::udaughter){
          spmap[sppi] = sppj;
          break;
        } else if (rel == MCRelationship::mother || rel == MCRelationship::umother){
          spmap[sppj] = sppi;
          break;
        } else if (rel == MCRelationship::sibling || rel == MCRelationship::usibling){
          double dist = (sppj->startPosition() - sppi->startPosition()).mag();
          if (dist < 10 && sppi->id().asInt() > sppj->id().asInt()) spmap[sppi] = sppj;
          if (dist < 10 && sppi->id().asInt() < sppj->id().asInt()) spmap[sppj] = sppi;
          if (dist < 10) break;
        }

      }
    }

    // check for remapping
    bool changed(true);
    while (_mapSimAncestor && changed){
      changed = false;
      for (auto im = spmap.begin();im!=spmap.end();++im){
        auto ifnd = spmap.find(im->second);
        if (ifnd->second != ifnd->first){
          changed = true;
          spmap[im->first] = ifnd->second;
        }
      }
    }

    // find the most likely SimParticle for this cluster
    std::map<int,int> mode;
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> spp = *ipp;
      spp = spmap[spp];
      ++mode[spp->id().asInt()];
    }

    int max(0);
    std::map<int,int>::iterator imax = mode.end();
    for (auto im=mode.begin();im!=mode.end();++im){
      icontrib.push_back(im->first);
      if (im->second>max){imax=im; max = im->second;}
    }

    unsigned pid = (imax != mode.end()) ? imax->first : 0;
    for (auto im = spmap.begin();im!=spmap.end();++im){
      if(im->first->id().asInt() != pid) continue;
      mptr = im->first;
      break;
    }

    // find the index of the mc digi for the ancestor SimParticle
    for (auto id : dids) {
      const auto& mcdigi = _mcdigis->at(id);
      const auto& sgsp   = mcdigi.earlyStrawGasStep();
      const auto& spp    = sgsp->simParticle();
      if (spp == mptr){shid = id; break;}
    }
  }



  //------------------------------------------------------------------------------------------------------------------------------
  void BkgDiag2::fillStrawHitInfoMC(StrawDigiMC const& mcdigi, art::Ptr<SimParticle>const& mptr, StrawHitInfo& shinfo) const {

    // use early step to define the MC match
    auto const& sgsp = mcdigi.earlyStrawGasStep();
    auto const& spp  = sgsp->simParticle();
    shinfo._mctime   = sgsp->time();
    shinfo._mcht     = mcdigi.wireEndTime(mcdigi.earlyEnd());
    shinfo._mcpdg    = spp->pdgId();
    shinfo._mcproc   = spp->creationCode();
    shinfo._mcedep   = mcdigi.energySum();
    shinfo._mcpos    = sgsp->position();
    shinfo._mcmom    = sqrt(sgsp->momentum().mag2());
    double cosd      = cos(sgsp->momentum().Theta());
    shinfo._mctd     = cosd/sqrt(1.0-cosd*cosd);


    // relationship to main particle
    shinfo._mrel = MCRelationship::none;
    if (!spp.isNonnull() || !mptr.isNonnull()) return;

    MCRelationship rel(spp,mptr);
    shinfo._mrel = rel.relationship();

    for (auto const& mcmptr : _mcprimary->primarySimParticles()){
      MCRelationship rel(spp,mcmptr);
      if (rel.relationship() == MCRelationship::none) continue;

      if(shinfo._prel > MCRelationship::none)
        shinfo._prel = std::min(shinfo._prel,(int)rel.relationship());
      else
        shinfo._prel = rel.relationship();
    }
  }


  //------------------------------------------------------------------------------------------------------------------------------
  void BkgDiag2::fillStrawHitInfo(size_t ich, StrawHitInfo& shinfo) const
  {
     auto const& ch  = _chcol->at(ich);
     auto const& shf = ch.flag();

     shinfo._stereo     = shf.hasAllProperties(StrawHitFlag::stereo);
     shinfo._tdiv       = shf.hasAllProperties(StrawHitFlag::tdiv);
     shinfo._esel       = shf.hasAllProperties(StrawHitFlag::energysel);
     shinfo._rsel       = shf.hasAllProperties(StrawHitFlag::radsel);
     shinfo._tsel       = shf.hasAllProperties(StrawHitFlag::timesel);
     shinfo._strawxtalk = shf.hasAllProperties(StrawHitFlag::strawxtalk);
     shinfo._elecxtalk  = shf.hasAllProperties(StrawHitFlag::elecxtalk);
     shinfo._isolated   = shf.hasAllProperties(StrawHitFlag::isolated);
     shinfo._bkg        = shf.hasAllProperties(StrawHitFlag::bkg);
     shinfo._bkgc       = shf.hasAllProperties(StrawHitFlag::bkgclust);

     shinfo._pos        = ch.pos();
     shinfo._time       = ch.correctedTime();
     shinfo._wdist      = ch.wireDist();
     shinfo._wres       = ch.posRes(ComboHit::wire);
     shinfo._tres       = ch.posRes(ComboHit::trans);
     shinfo._chisq      = ch.qual();
     shinfo._edep       = ch.energyDep();

     StrawId const& sid = ch.strawId();
     shinfo._plane      = sid.plane();
     shinfo._panel      = sid.panel();
     shinfo._layer      = sid.layer();
     shinfo._straw      = sid.straw();
     shinfo._stereo     = ch.flag().hasAllProperties(StrawHitFlag::stereo);
     shinfo._tdiv       = ch.flag().hasAllProperties(StrawHitFlag::tdiv);
  }
}

using mu2e::BkgDiag2;
DEFINE_ART_MODULE(BkgDiag2)
