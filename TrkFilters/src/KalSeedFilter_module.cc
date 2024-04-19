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

      struct Config{
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>      kalSeedCollection {     Name("kalSeedCollection"),       Comment("kalSeedCollection ") };
        fhicl::Atom<bool>               requireCaloCluster{     Name("requireCaloCluster"),      Comment("requireCaloCluster") };
        fhicl::Atom<int>                fitparticle       {     Name("fitparticle"),             Comment("fitparticle       ") };
        fhicl::Atom<int>                fitdirection      {     Name("fitdirection"),            Comment("fitdirection      ") };
        fhicl::Atom<double>             minFitCons        {     Name("minFitCons"),              Comment("minFitCons        ") };
        fhicl::Atom<double>             minNHits          {     Name("minNStrawHits"),           Comment("minNStrawHits     ") };
        fhicl::Atom<double>             minMomentum       {     Name("minMomentum"),             Comment("minMomentum       ") };
        fhicl::Atom<double>             maxMomentum       {     Name("maxMomentum"),             Comment("maxMomentum       ") };
        fhicl::Atom<double>             minTanDip         {     Name("minTanDip"),               Comment("minTanDip         ") };
        fhicl::Atom<double>             maxTanDip         {     Name("maxTanDip"),               Comment("maxTanDip         ") };
        fhicl::Atom<double>             maxChi2DOF        {     Name("maxChi2DOF"),              Comment("maxChi2DOF        ") };
        fhicl::Atom<double>             maxMomErr         {     Name("maxMomErr"),               Comment("maxMomErr         ") };
        fhicl::Atom<double>             minD0             {     Name("minD0"),                   Comment("minD0             ") };
        fhicl::Atom<double>             maxD0             {     Name("maxD0"),                   Comment("maxD0             ") };
        fhicl::Atom<double>             minT0             {     Name("minT0"),                   Comment("minT0             ") };
        fhicl::Sequence<std::string>    seedFitFlag       {     Name("seedFitFlag"),             Comment("seedFitFlag       ") , std::vector<std::string>{"SeedOK"}};
        fhicl::Atom<int>                debugLevel        {     Name("debugLevel"),             Comment("debugLevel        ") , 0};


      };

      using Parameters = art::EDFilter::Table<Config>;

      explicit     KalSeedFilter(const Parameters& config);
      virtual bool filter(art::Event& event) override;
      virtual bool endRun( art::Run& run ) override;

    private:
      art::InputTag   _ksTag;
      bool            _hascc; // Calo Cluster
      PDGCode::type     _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      double          _minfitcons;
      unsigned        _minnhits;
      double          _minmom, _maxmom, _mintdip, _maxtdip, _maxchi2dof, _maxmomerr;
      double          _minD0, _maxD0; // impact parameter limits
      double          _minT0;
      TrkFitFlag      _goods; // helix fit flag
      int             _debug;
      // counters
      unsigned        _nevt, _npass;
  };

  KalSeedFilter::KalSeedFilter(const Parameters& config):
    art::EDFilter{config},
    _ksTag     (config().kalSeedCollection()),
    _hascc     (config().requireCaloCluster()),
    _tpart     ((PDGCode::type)(config().fitparticle())),
    _fdir      ((TrkFitDirection::FitDirection)(config().fitdirection())),
    _minfitcons(config().minFitCons()),
    _minnhits  (config().minNHits()),
    _minmom    (config().minMomentum()),
    _maxmom    (config().maxMomentum()),
    _mintdip   (config().minTanDip()),
    _maxtdip   (config().maxTanDip()),
    _maxchi2dof(config().maxChi2DOF()),
    _maxmomerr (config().maxMomErr()),
    _minD0     (config().minD0()),
    _maxD0     (config().maxD0()),
    _minT0     (config().minT0()),
    _goods     (config().seedFitFlag()),
    _debug     (config().debugLevel()),
    _nevt(0), _npass(0)
    {
      produces<TriggerInfo>();
    }

  bool KalSeedFilter::filter(art::Event& evt){
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    const KalSeedCollection* kscol = ksH.product();
    // loop over the collection: if any pass the selection, pass this event
    if(_debug > 2){
      if (kscol->size()>0) printf("[KalSeedFilter::filter]   nhits     mom     momErr    chi2ndof     fitCon   tanDip    d0      \n");
    }
    for(auto iks = kscol->begin(); iks != kscol->end(); ++iks) {
      auto const& ks = *iks;
      if( ks.status().hasAllProperties(_goods) && ks.intersections().size()>0){
        // get the first intersection
        auto const& kinter = ks.intersections().front();
        // convert to LoopHeix; this is needed to get d0 cut (which shoudl be replaced FIXME)
        KinKal::LoopHelix lh(kinter.pstate_, KinKal::VEC3(kinter.bnom_));
        auto momvec = kinter.momentum3();
        double td = 1.0/tan(momvec.Theta());
        double mom = momvec.R();
        double d0 = lh.minAxisDist();
        //check particle type and fitdirection
        if ( (ks.particle() != _tpart) || (momvec.Z()*_fdir.dzdt() < 0))       continue;
        unsigned nactive = ks.nHits(true); //count active hits
        if(_debug > 2){
          printf("[KalSeedFilter::filter] %4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n", nactive, mom, kinter.momerr(),ks.chisquared()/ks.nDOF(), ks.fitConsistency(), td, d0);
        }
        if( (!_hascc || ks.caloCluster().isNonnull()) &&
            nactive >= _minnhits &&
            mom > _minmom && mom < _maxmom && kinter.momerr() < _maxmomerr &&
            ks.chisquared()/ks.nDOF() < _maxchi2dof && // chisq/ndof isn't a statistically robust measure, this cut should be removed FIXME
            ks.fitConsistency()  > _minfitcons &&
            td > _mintdip && td < _maxtdip &&
            d0 > _minD0 && d0 < _maxD0) { // d0 is not a global geometric parameter and signed cuts have particle-species dependence.   This should be replaced with a selection based on consistency with the stopping target FIXME
          retval = true;
          ++_npass;
          // Fill the trigger info object
          // associate to the helix which triggers.  Note there may be other helices which also pass the filter
          // but filtering is by event!
          triginfo->_tracks.push_back(KalSeed(ks));
          if(_debug > 1){
            std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
          }
        }
      }
    }
    evt.put(std::move(triginfo));
    return retval;
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
