//
//  Filter for selecting good seed (chisq) track fits: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
// mu2e
// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

//using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

namespace mu2e
{
  class KalSeedFilter : public art::EDFilter
  {
  public:
    struct KalSeedCutsConfig{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool>               requireCaloCluster  {     Name("requireCaloCluster"),      Comment("requireCaloCluster ")};
      fhicl::OptionalAtom<int>                fitparticle         {     Name("fitparticle"),             Comment("fitparticle       ") };
      fhicl::OptionalAtom<std::string>        fitdirection        {     Name("fitdirection"),            Comment("fitdirection (\"downstream\" or \"upstream\")") };
      fhicl::Atom<double>             minFitCons          {     Name("minFitCons"),              Comment("minFitCons        ") };
      fhicl::Atom<double>             minNHits            {     Name("minNStrawHits"),           Comment("minNStrawHits     ") };
      fhicl::Atom<double>             minMomentum         {     Name("minMomentum"),             Comment("minMomentum       ") };
      fhicl::Atom<double>             maxMomentum         {     Name("maxMomentum"),             Comment("maxMomentum       ") };
      fhicl::Atom<double>             minTanDip           {     Name("minTanDip"),               Comment("minTanDip         ") };
      fhicl::Atom<double>             maxTanDip           {     Name("maxTanDip"),               Comment("maxTanDip         ") };
      fhicl::Atom<double>             maxChi2DOF          {     Name("maxChi2DOF"),              Comment("maxChi2DOF        ") };
      fhicl::Atom<double>             maxMomErr           {     Name("maxMomErr"),               Comment("maxMomErr         ") };
      fhicl::Atom<double>             minD0               {     Name("minD0"),                   Comment("minD0             ") };
      fhicl::Atom<double>             maxD0               {     Name("maxD0"),                   Comment("maxD0             ") };
      fhicl::Atom<unsigned>           minNStereo          {     Name("minNStereo"),              Comment("Min number of 12 possible panel orientations on track"),0};
      fhicl::Atom<unsigned>           minNPlanes          {     Name("minNPlanes"),              Comment("Min number of planes hit "),0};
      fhicl::Sequence<std::string>    seedFitFlag         {     Name("seedFitFlag"),             Comment("seedFitFlag       ") , std::vector<std::string>{"SeedOK"}};
    };

    struct KalSeedCutsTool {
      KalSeedCutsTool(const KalSeedCutsConfig& config):
        _hascc     (config.requireCaloCluster()),
        _minfitcons(config.minFitCons()),
        _minnhits  (config.minNHits()),
        _minmom    (config.minMomentum()),
        _maxmom    (config.maxMomentum()),
        _mintdip   (config.minTanDip()),
        _maxtdip   (config.maxTanDip()),
        _maxchi2dof(config.maxChi2DOF()),
        _maxmomerr (config.maxMomErr()),
        _minD0     (config.minD0()),
        _maxD0     (config.maxD0()),
        _minnstereo(config.minNStereo()),
        _minnplanes(config.minNPlanes()),
        _goods     (config.seedFitFlag())
      {
        int tpart;
        _doParticleTypeCheck =config.fitparticle(tpart);
        if(_doParticleTypeCheck)_tpart = (PDGCode::type)tpart;
        std::string fdir;
        _doZPropDirCheck = config.fitdirection(fdir);
        if(_doZPropDirCheck)_fdir = TrkFitDirection::fitDirectionFromName(fdir);
      }

      KalSeedCutsTool() {}

      bool            _hascc; // Calo Cluster
      PDGCode::type   _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      double          _minfitcons;
      unsigned        _minnhits;
      double          _minmom, _maxmom, _mintdip, _maxtdip, _maxchi2dof, _maxmomerr;
      double          _minD0, _maxD0; // impact parameter limits
      unsigned        _minnstereo, _minnplanes;
      TrkFitFlag      _goods; // helix fit flag
      bool            _doParticleTypeCheck;
      bool            _doZPropDirCheck;
    };

    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Sequence<art::InputTag>          kalSeedCollections { Name("kalSeedCollections"),      Comment("kalSeedCollections ") };
      fhicl::Sequence<fhicl::Table<KalSeedCutsConfig>>  KalSeedCuts       { Name("KalSeedCuts"),            Comment("Cuts applied to the KalSeeds")};
      fhicl::Atom<int>                    debugLevel        { Name("debugLevel"),             Comment("debugLevel        ") , 0};
      fhicl::Atom<bool>                   noFilter          { Name("noFilter"),               Comment("don't apply any filter decision") , false};
      fhicl::Atom<bool>                   noInfo            { Name("noInfo"),                 Comment("don't create TriggerInfo object") , false};
      fhicl::Atom<unsigned>               minNTrks          { Name("minNTrks"),               Comment("minimum number of tracks passing the selection") , 1};
      fhicl::OptionalAtom<std::string>    surface           { Name("momSurface"),             Comment("Surface at which to compare fits. If unset use t0 segment")};
    };

    using Parameters = art::EDFilter::Table<Config>;

    explicit     KalSeedFilter(const Parameters& config);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;
    bool   checkKalSeed(const KalSeed&Ks, const KalSeedCutsTool&Cuts);

  private:
    std::vector<art::InputTag>   _ksTags;
    std::vector<KalSeedCutsConfig> _ksCutsConfig;
    std::vector<KalSeedCutsTool>   _ksCuts;
    int             _debug;
    bool            _noFilter, _noInfo;
    unsigned        _minNTrks;
    bool            _intmom; // check momentum at intersection? if not, use t0 segment
    SurfaceId       _momsid; //Surface for momentum test
    // counters
    unsigned        _nevt, _npass;
  };

  KalSeedFilter::KalSeedFilter(const Parameters& config):
    art::EDFilter{config},
    _ksTags      (config().kalSeedCollections()),
    _ksCutsConfig(config().KalSeedCuts()),
    _debug       (config().debugLevel()),
    _noFilter    (config().noFilter()),
    _noInfo      (config().noInfo()),
    _minNTrks    (config().minNTrks()),
    _nevt(0), _npass(0)
    {
      for (auto const&cf: _ksCutsConfig){
        _ksCuts.push_back(KalSeedCutsTool(cf));
      }
      if(!_noInfo)produces<TriggerInfo>();
      std::string momsurf;
      _intmom = config().surface(momsurf);
      if(_intmom)_momsid = SurfaceId(momsurf);
    }

  bool KalSeedFilter::filter(art::Event& evt){
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    unsigned nGoodTrks(0);
    // Loop over the collection
    for( auto kstag : _ksTags) {
      auto ksH = evt.getValidHandle<KalSeedCollection>(kstag);
      const KalSeedCollection* kscol = ksH.product();
      // loop over the collection: if any pass the selection, pass this event
      if(_debug > 1){
        std::cout << moduleDescription().moduleLabel() << " input " << kstag.encode() << " collection has " << kscol->size() << " tracks\n";
      }
      if(_debug > 2){
        if (kscol->size()>0) printf("[KalSeedFilter::filter]   nhits nst npl    mom     momErr    chi2ndof     fitCon   tanDip    d0      \n");
      }
      for(auto iks = kscol->begin(); iks != kscol->end(); ++iks) {
        auto const& ks = *iks;

        for (auto const&cuts : _ksCuts){
          if (checkKalSeed(ks, cuts)) {
            ++nGoodTrks;
            ++_npass;
            // Fill the trigger info object
            // associate to the helix which triggers.  Note there may be other helices which also pass the filter
            // but filtering is by event!
            size_t index = std::distance(kscol->begin(),iks);
            if(!_noInfo)triginfo->_tracks.push_back(art::Ptr<KalSeed>(ksH,index));

            if(_debug > 1){
              std::cout << moduleDescription().moduleLabel() << " --> accepted a track" << std::endl;
            }
            break;//no need to check the other ksCuts entries
          }
        }//end loop over the ksCuts
      }//end loop over the kalseeds
    }// end loop over KalSeed collections

    if(_debug > 1){
      std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
    }
    if(!_noInfo)evt.put(std::move(triginfo));

    if (!_noFilter){
      return (nGoodTrks >= _minNTrks);
    }else {
      return true;
    }
  }

  bool KalSeedFilter::checkKalSeed(const KalSeed&Ks, const KalSeedCutsTool&Cuts){
    if(_debug > 3){
      std::cout << "KalSeedFilter: in checkKalSeed status "<< Ks.status() << std::endl;
    }

    if( Ks.status().hasAllProperties(Cuts._goods) ){

      // extract test quantities from the fit segment at t0
      double t0;
      auto t0seg = Ks.t0Segment(t0);
      if(t0seg == Ks.segments().end()) return false;
      auto momvec = t0seg->momentum3();
      auto posvec = t0seg->position3();
      if(_intmom){
        auto kintercol = Ks.intersections(_momsid);
        for(auto kinter : kintercol) {
          if(kinter->momentum3().Z() > 0.0){
            momvec = kinter->momentum3();
            posvec = kinter->position3();
            break;
          }
        }
      }
      double td     = 1.0/tan(momvec.Theta());
      double mom    = momvec.R();
      double d0(0.0);
      if(Ks.loopHelixFit()){
        auto lhtraj = t0seg->loopHelix();
        d0 = std::copysign(lhtraj.minAxisDist(),posvec.Cross(momvec).Z());
      } else if (Ks.centralHelixFit()) {
        auto chtraj = t0seg->centralHelix();
        d0 = chtraj.d0();
      } else if (Ks.kinematicLineFit()){
        auto kltraj = t0seg->kinematicLine();
        d0 = kltraj.d0();
      } else {
        return false;
      }

      //check particle type and fitdirection
      if ( Cuts._doParticleTypeCheck){
        if (Ks.particle() != Cuts._tpart)  {
          if(_debug > 2)std::cout << "KalSeedFilter: particle cut failed" << std::endl;
          return false;
        }
      }
      if (Cuts._doZPropDirCheck){
        if (momvec.Z()*Cuts._fdir.dzdt() < 0) {
          if(_debug > 2)std::cout << "KalSeedFilter: dir cut failed" << std::endl;
          return false;
        }
      }
      unsigned nactive = Ks.nHits(true); //count active hits

      // compute number of planes and unique stereo orientations
      std::set<unsigned> stcount;
      std::set<unsigned> pcount;
      for ( auto const& hit : Ks.hits()){
        if(hit._flag.hasAllProperties(StrawHitFlag::active)){
          stcount.insert(hit._sid.stereoPanel());
          pcount.insert(hit._sid.plane());
        }
      }

      if(_debug > 2){
        printf("[KalSeedFilter::filter] %4d %4lu %4lu %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",
               nactive, stcount.size(), pcount.size(), mom, t0seg->momerr(),Ks.chisquared()/Ks.nDOF(), Ks.fitConsistency(), td, d0);
      }

      if( (!Cuts._hascc || Ks.caloCluster().isNonnull()) &&
          nactive >= Cuts._minnhits &&
          stcount.size() >= Cuts._minnstereo &&
          pcount.size() >= Cuts._minnplanes &&
          mom > Cuts._minmom && mom < Cuts._maxmom && t0seg->momerr() < Cuts._maxmomerr &&
          Ks.chisquared()/Ks.nDOF() < Cuts._maxchi2dof && // chisq/ndof isn't a statistically robust measure, this cut should be removed FIXME
          Ks.fitConsistency()       > Cuts._minfitcons &&
          td > Cuts._mintdip && td < Cuts._maxtdip &&
          d0 > Cuts._minD0   && d0 < Cuts._maxD0) {
        if(_debug > 1) std::cout << "Selected track " << std::endl;
        return true;
      } else if(_debug > 2){
        std::cout << "KalSeedFilter: parameter cuts failed: nactive " << std::endl;
      }
    } else {
      if(_debug > 2)std::cout << "KalSeedFilter: basic cuts failed: status "<< Ks.status() << " intersections " << Ks.intersections().size() << std::endl;
    }
    return false;
  }

  bool KalSeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::KalSeedFilter;
DEFINE_ART_MODULE(KalSeedFilter)
