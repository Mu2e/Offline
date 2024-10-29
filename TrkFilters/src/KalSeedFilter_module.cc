//
//  Filter for selecting good seed (chisq) track fits: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "fhiclcpp/types/Atom.h"
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
      fhicl::Atom<bool>               doParticleTypeCheck {     Name("doParticleTypeCheck"),     Comment("doParticleTypeCheck")};
      fhicl::Atom<int>                fitparticle         {     Name("fitparticle"),             Comment("fitparticle       ") };
      fhicl::Atom<bool>               doZPropDirCheck     {     Name("doZPropDirCheck"),         Comment("doZPropDirCheck   ") };
      fhicl::Atom<int>                fitdirection        {     Name("fitdirection"),            Comment("fitdirection      ") };
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
      fhicl::Atom<double>             minT0               {     Name("minT0"),                   Comment("minT0             ") };
      fhicl::Sequence<std::string>    seedFitFlag         {     Name("seedFitFlag"),             Comment("seedFitFlag       ") , std::vector<std::string>{"SeedOK"}};
    };

    struct KalSeedCutsTool {
      KalSeedCutsTool(const KalSeedCutsConfig& config):
        _hascc     (config.requireCaloCluster()),
        _tpart     ((PDGCode::type)config.fitparticle()),
        _fdir      ((TrkFitDirection::FitDirection)config.fitdirection()),
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
        _minT0     (config.minT0()),
        _goods     (config.seedFitFlag()),
        _doParticleTypeCheck(config.doParticleTypeCheck()),
        _doZPropDirCheck    (config.doZPropDirCheck()) {}

      KalSeedCutsTool() {}

      bool            _hascc; // Calo Cluster
      PDGCode::type   _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      double          _minfitcons;
      unsigned        _minnhits;
      double          _minmom, _maxmom, _mintdip, _maxtdip, _maxchi2dof, _maxmomerr;
      double          _minD0, _maxD0; // impact parameter limits
      double          _minT0;
      TrkFitFlag      _goods; // helix fit flag
      bool            _doParticleTypeCheck;
      bool            _doZPropDirCheck;
    };

    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>          kalSeedCollection { Name("kalSeedCollection"),      Comment("kalSeedCollection ") };
      fhicl::Sequence<fhicl::Table<KalSeedCutsConfig>>  KalSeedCuts       { Name("KalSeedCuts"),            Comment("Cuts applied to the KalSeeds")};
      fhicl::Atom<int>                    debugLevel        { Name("debugLevel"),             Comment("debugLevel        ") , 0};
      fhicl::Atom<bool>                   noFilter          { Name("noFilter"),               Comment("don't apply any filter decision") , 0};
      fhicl::Atom<unsigned>               minNTrks          { Name("minNTrks"),               Comment("minimum number of tracks passing the selection") , 1};
    };

    using Parameters = art::EDFilter::Table<Config>;

    explicit     KalSeedFilter(const Parameters& config);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;
    bool   checkKalSeed(const KalSeed&Ks, const KalSeedCutsTool&Cuts);

  private:
    art::InputTag   _ksTag;
    std::vector<KalSeedCutsConfig> _ksCutsConfig;
    std::vector<KalSeedCutsTool>   _ksCuts;
    int             _debug;
    bool            _noFilter;
    unsigned        _minNTrks;

    // counters
    unsigned        _nevt, _npass;
  };

  KalSeedFilter::KalSeedFilter(const Parameters& config):
    art::EDFilter{config},
    _ksTag       (config().kalSeedCollection()),
    _ksCutsConfig(config().KalSeedCuts()),
    _debug       (config().debugLevel()),
    _noFilter    (config().noFilter()),
    _minNTrks    (config().minNTrks()),
    _nevt(0), _npass(0)
    {
      for (auto const&cf: _ksCutsConfig){
        _ksCuts.push_back(KalSeedCutsTool(cf));
      }
      produces<TriggerInfo>();
    }

  bool KalSeedFilter::filter(art::Event& evt){
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    unsigned nGoodTrks(0);
    // find the collection
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    const KalSeedCollection* kscol = ksH.product();
    // loop over the collection: if any pass the selection, pass this event
    if(_debug > 2){
      if (kscol->size()>0) printf("[KalSeedFilter::filter]   nhits     mom     momErr    chi2ndof     fitCon   tanDip    d0      \n");
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
          triginfo->_tracks.push_back(art::Ptr<KalSeed>(ksH,index));

          if(_debug > 1){
            std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
          }
          break;//no need to check the other ksCuts entries
        }
      }//end loop over the ksCuts
    }//end loop over the kalseeds

    evt.put(std::move(triginfo));

    if (!_noFilter){
      return (nGoodTrks >= _minNTrks);
    }else {
      return true;
    }
  }

  bool KalSeedFilter::checkKalSeed(const KalSeed&Ks, const KalSeedCutsTool&Cuts){
    if( Ks.status().hasAllProperties(Cuts._goods) && Ks.intersections().size()>0){
      // get the first intersection
      auto const& kinter = Ks.intersections().front();
      // convert to LoopHeix; this is needed to get d0 cut (which shoudl be replaced FIXME)
      KinKal::LoopHelix lh(kinter.pstate_, KinKal::VEC3(kinter.bnom_));
      auto   momvec = kinter.momentum3();
      double td     = 1.0/tan(momvec.Theta());
      double mom    = momvec.R();
      double d0     = lh.minAxisDist();
      //check particle type and fitdirection
      if ( Cuts._doParticleTypeCheck){
        if (Ks.particle() != Cuts._tpart)  {
          return false;
        }
      }
      if (Cuts._doZPropDirCheck){
        if (momvec.Z()*Cuts._fdir.dzdt() < 0) {
          return false;
        }
      }
      unsigned nactive = Ks.nHits(true); //count active hits
      if(_debug > 2){
        printf("[KalSeedFilter::filter] %4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",
               nactive, mom, kinter.momerr(),Ks.chisquared()/Ks.nDOF(), Ks.fitConsistency(), td, d0);
      }
      if( (!Cuts._hascc || Ks.caloCluster().isNonnull()) &&
          nactive >= Cuts._minnhits &&
          mom > Cuts._minmom && mom < Cuts._maxmom && kinter.momerr() < Cuts._maxmomerr &&
          Ks.chisquared()/Ks.nDOF() < Cuts._maxchi2dof && // chisq/ndof isn't a statistically robust measure, this cut should be removed FIXME
          Ks.fitConsistency()       > Cuts._minfitcons &&
          td > Cuts._mintdip && td < Cuts._maxtdip &&
          d0 > Cuts._minD0   && d0 < Cuts._maxD0) { // d0 is not a global geometric parameter and signed cuts have particle-species dependence.   This should be replaced with a selection based on consistency with the stopping target FIXME

        return true;
      }
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
