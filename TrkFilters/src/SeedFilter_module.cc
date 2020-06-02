//
//  Filter for selecting good seed (chisq) track fits: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "fhiclcpp/ParameterSet.h"
// mu2e
// data
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
// BTrk
#include "BTrk/TrkBase/TrkParticle.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

namespace mu2e
{
  class SeedFilter : public art::EDFilter
  {
  public:
    explicit SeedFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag   _ksTag;
    bool            _hascc; // Calo Cluster
    TrkParticle     _tpart; // particle type being searched for
    TrkFitDirection _fdir;  // fit direction in search
    double          _minfitcons;
    unsigned        _minnhits;
    double          _minmom, _maxmom, _mintdip, _maxtdip, _maxchi2dof, _maxmomerr;
    double          _minD0, _maxD0; // impact parameter limits
    double          _minT0;
    TrkFitFlag      _goods; // helix fit flag
    std::string     _trigPath;
    int             _debug;
    // counters
    unsigned        _nevt, _npass;
  };

  SeedFilter::SeedFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _ksTag     (pset.get<art::InputTag>("kalSeedCollection","KSFDeM")),
    _hascc     (pset.get<bool>("requireCaloCluster",false)),
    _tpart     ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir      ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _minfitcons(pset.get<double>("minFitCons",-1.)),   //not used by default
    _minnhits  (pset.get<unsigned>("minNHits",15)),
    _minmom    (pset.get<double>("minMomentum",40.0)),
    _maxmom    (pset.get<double>("maxMomentum",200.0)) ,
    _mintdip   (pset.get<double>("minTanDip", 0.)),       //not used by default. 0.57735027
    _maxtdip   (pset.get<double>("maxTanDip", 100.)),     //not used by default. 1.5574077
    _maxchi2dof(pset.get<double>("maxChi2DOF",20.0)),
    _maxmomerr (pset.get<double>("maxMomErr",10)),
    _minD0     (pset.get<double>("minD0",-200.)),
    _maxD0     (pset.get<double>("maxD0", 200.)),
    _minT0     (pset.get<double>("minT0", 0.)),
    _goods     (pset.get<std::vector<std::string> >("seedFitFlag",std::vector<std::string>{"SeedOK"})),
    _trigPath  (pset.get<std::string>("triggerPath")),
    _debug     (pset.get<int>   ("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool SeedFilter::filter(art::Event& evt){
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    const KalSeedCollection* kscol = ksH.product();
    size_t trig_ind(0);
    // loop over the collection: if any pass the selection, pass this event
    for(auto iks = kscol->begin(); iks != kscol->end(); ++iks) {
      auto const& ks = *iks;
      //check particle type and fitdirection
      if ( (ks.particle() != _tpart) || (ks.fitDirection() != _fdir))       continue;

      // I should not be calculating NDOF here, this should be in an adapter, FIXME!!
      unsigned nactive(0);
      for(auto const& ish : ks.hits())
        if(ish.flag().hasAllProperties(StrawHitFlag::active))++nactive;
      float ndof = std::max(1.0,nactive - 5.0);
      // get the first segment
      KalSegment const& fseg = ks.segments().front();
      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << "status = " << ks.status() << " nactive = " << nactive << " mom = " << fseg.mom() << " chisq/dof = " << ks.chisquared()/ndof << std::endl;
      }
      if( ks.status().hasAllProperties(_goods) &&
          (!_hascc || ks.caloCluster().isNonnull()) &&
          nactive >= _minnhits &&
          fseg.mom() > _minmom && fseg.mom() < _maxmom && fseg.momerr() < _maxmomerr &&
          ks.chisquared()/ndof < _maxchi2dof &&
          ks.fitConsistency()  > _minfitcons &&
          fseg.helix().tanDip() > _mintdip && fseg.helix().tanDip() < _maxtdip &&
          fseg.helix().d0() > _minD0 && fseg.helix().d0() < _maxD0 ) {
        retval = true;
        ++_npass;
        // Fill the trigger info object
        if (trig_ind == 0){
	  triginfo->_triggerBits.merge(TriggerFlag::track);
	  triginfo->_triggerPath = _trigPath;
	}
        // associate to the helix which triggers.  Note there may be other helices which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(kscol->begin(),iks);
	triginfo->_tracks.push_back(art::Ptr<KalSeed>(ksH,index));
	++trig_ind;
        if(_debug > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
        }
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool SeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::SeedFilter;
DEFINE_ART_MODULE(SeedFilter);
