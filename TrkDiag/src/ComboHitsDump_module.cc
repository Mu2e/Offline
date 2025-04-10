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

  class ComboHitsDump : public art::EDAnalyzer {
    public:

      struct Config {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>    ComboHitCollection       {Name("ComboHitCollection"),        Comment("ComboHit collection name")             };
        fhicl::Atom<art::InputTag>    StrawDigiMCCollection    {Name("StrawDigiMCCollection"),     Comment("StrawDigiMC collection name")          };
        fhicl::Atom<bool>             mcdiag                   {Name("MCDiag"),                    Comment("Monte Carlo Diag")                     };
     };

      explicit ComboHitsDump(const art::EDAnalyzer::Table<Config>& config);
      virtual ~ComboHitsDump();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);


    private:
      void fillStrawHitInfo(size_t ich, StrawHitInfo& bkghinfo) const;
      void fillStrawHitInfoMC(StrawDigiMC const& mcdigi, StrawHitInfo& shinfo) const;
      bool findData(const art::Event& e);


      art::ProductToken<ComboHitCollection>      _chToken;
      art::ProductToken<BkgClusterCollection>    _bkgcToken;
      art::ProductToken<BkgClusterHitCollection> _bkghToken;
      art::ProductToken<StrawDigiMCCollection>   _mcdigisToken;
      bool _mcdiag;

      // cache of event objects
      const ComboHitCollection      *_chcol;
      const BkgClusterHitCollection *_bkghitcol;
      const StrawDigiMCCollection   *_mcdigis;

      //Ntuple -  All hits variables
      TTree              *_hitDiag;
      Int_t              _iev, _nhits;
      std::vector<int>   _hitPdg,_hitproc,_hitnch,_hitnsh,_hitid;
      std::vector<float> _hitPosX,_hitPosY,_hitPosZ,_hitTime,_hiteDep,_hitDirX;
      std::vector<float> _hitDirY,_hitDirZ,_hitDirQ,_hitMCmomX,_hitMCmomY,_hitMCmomZ;
      std::vector<BkgHitInfo> _bkghinfo;
  };


  ComboHitsDump::ComboHitsDump(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _chToken       { consumes<ComboHitCollection>(config().ComboHitCollection() ) },
    _mcdigisToken  { consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection() ) },
    _mcdiag        ( config().mcdiag() )
  {}

  ComboHitsDump::~ComboHitsDump(){}

  void ComboHitsDump::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;


    _hitDiag = tfs->make<TTree>("hitDiag","hit diagnostics");
    _hitDiag->Branch("iev",     &_iev,       "iev/I");
    _hitDiag->Branch("nhits",   &_nhits,     "nhits/I");
    _hitDiag->Branch("posX",    &_hitPosX    );
    _hitDiag->Branch("posY",    &_hitPosY    );
    _hitDiag->Branch("posZ",    &_hitPosZ    );
    _hitDiag->Branch("dirX",    &_hitDirX    );
    _hitDiag->Branch("dirY",    &_hitDirY    );
    _hitDiag->Branch("dirZ",    &_hitDirZ    );
    _hitDiag->Branch("dirQ",    &_hitDirQ    );
    _hitDiag->Branch("time",    &_hitTime    );
    _hitDiag->Branch("eDep",    &_hiteDep    );
    _hitDiag->Branch("nch",     &_hitnch     );
    _hitDiag->Branch("nsh",     &_hitnsh     );
    _hitDiag->Branch("bkghinfo",&_bkghinfo   );

    if (_mcdiag){
      _hitDiag->Branch("mcpdg",  &_hitPdg     );
      _hitDiag->Branch("mcproc", &_hitproc    );
      _hitDiag->Branch("mcid",   &_hitid      );
      _hitDiag->Branch("mcmomX", &_hitMCmomX  );
      _hitDiag->Branch("mcmomY", &_hitMCmomY  );
      _hitDiag->Branch("mcmomZ", &_hitMCmomZ  );
    }

  }




  void ComboHitsDump::analyze(const art::Event& event ) {

    if (!findData(event))
      throw cet::exception("RECO")<<"mu2e::ComboHitsDump: data missing or incomplete"<< std::endl;



    // Ntuple - all hits
    //----------------------------------
      _hitPosX.clear(); _hitPosY.clear();_hitPosZ.clear();
      _hitDirX.clear(); _hitDirY.clear();_hitDirZ.clear();_hitDirQ.clear();
      _hitTime.clear(); _hiteDep.clear();_hitMCmomX.clear();_hitMCmomY.clear();_hitMCmomZ.clear();
      _hitnch.clear(); _hitnsh.clear();_hitPdg.clear();_hitproc.clear();_hitid.clear();
      _bkghinfo.clear();


      for (size_t ich=0;ich<_chcol->size();++ich){
        const auto& ch   = _chcol->at(ich);
        //const auto& shf  = _bkghitcol->at(ich).flag();

        _hitPosX.push_back(ch.pos().x());
        _hitPosY.push_back(ch.pos().y());
        _hitPosZ.push_back(ch.pos().z());
        _hitDirX.push_back(ch.hDir().x());
        _hitDirY.push_back(ch.hDir().y());
        _hitDirZ.push_back(ch.hDir().z());
        _hitDirQ.push_back(ch.qual());
        _hitTime.push_back(ch.time());
        _hiteDep.push_back(ch.energyDep());
        _hitnch.push_back(ch.nCombo());
        _hitnsh.push_back(ch.nStrawHits());

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


        BkgHitInfo bkghinfo;
        fillStrawHitInfo(ich,bkghinfo);

        if (_mcdiag) {
          std::vector<StrawDigiIndex> dids;
          _chcol->fillStrawDigiIndices(ich,dids);
          StrawDigiMC const& mcdigi = _mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
          fillStrawHitInfoMC(mcdigi,bkghinfo);
        }
        _bkghinfo.push_back(bkghinfo);


        ++_nhits;
      }

    _hitDiag->Fill();

    ++_iev;
  }



  bool ComboHitsDump::findData(const art::Event& evt){
    _chcol = 0; _bkghitcol = 0; _mcdigis = 0;
    auto chH = evt.getValidHandle(_chToken);
    _chcol = chH.product();
    auto bkghH = evt.getValidHandle(_bkghToken);
    _bkghitcol = bkghH.product();

    if(_mcdiag){
      auto mcdH = evt.getValidHandle(_mcdigisToken);
      _mcdigis = mcdH.product();
    }

    return _chcol != 0 && _bkghitcol != 0 && ( _mcdigis != 0 || !_mcdiag);
  }






  void ComboHitsDump::fillStrawHitInfoMC(StrawDigiMC const& mcdigi, StrawHitInfo& shinfo) const {

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
    shinfo._mrel     = MCRelationship::none;
  }


  void ComboHitsDump::fillStrawHitInfo(size_t ich, StrawHitInfo& shinfo) const
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

using mu2e::ComboHitsDump;
DEFINE_ART_MODULE(ComboHitsDump)
