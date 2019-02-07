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

using namespace std;

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
    int             _debug;
    // counters
    unsigned        _nevt, _npass;
  };

  SeedFilter::SeedFilter(fhicl::ParameterSet const& pset) :
    _ksTag     (pset.get<art::InputTag>("KalSeedCollection","KSFDeM")),
    _hascc     (pset.get<bool>("RequireCaloCluster",false)),
    _tpart     ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir      ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _minfitcons(pset.get<double>("MinFitCons",-1.)),   //not used by default
    _minnhits  (pset.get<unsigned>("MinNHits",15)),
    _minmom    (pset.get<double>("MinMomentum",40.0)),
    _maxmom    (pset.get<double>("MaxMomentum",200.0)) ,
    _mintdip   (pset.get<double>("MinTanDip", 0.)),       //not used by default. 0.57735027
    _maxtdip   (pset.get<double>("MaxTanDip", 100.)),     //not used by default. 1.5574077
    _maxchi2dof(pset.get<double>("MaxChi2DOF",20.0)),
    _maxmomerr (pset.get<double>("MaxMomErr",10)),
    _minD0     (pset.get<double>("MinD0",-200.)),
    _maxD0     (pset.get<double>("MaxD0", 200.)),
    _minT0     (pset.get<double>("MinT0", 0.)),
    _goods     (pset.get<vector<string> >("SeedFitFlag",vector<string>{"SeedOK"})),
    _debug     (pset.get<int>("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool SeedFilter::filter(art::Event& evt){
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    const KalSeedCollection* kscol = ksH.product();
    // loop over the collection: if any pass the selection, pass this event
    for(auto iks = kscol->begin(); iks != kscol->end(); ++iks) {
      auto const& ks = *iks;
      //check particle type and fitdirection
      if ( (ks.particle() != _tpart) || (ks.fitDirection() != _fdir))       continue;

      // I should not be calculating NDOF here, this should be in an adapter, FIXME!!
      unsigned nactive(0);
      for(auto const& ish : ks.hits())
        if(ish.flag().hasAllProperties(StrawHitFlag::active))++nactive;
      float ndof = max(1.0,nactive - 5.0);
      // get the first segment
      KalSegment const& fseg = ks.segments().front();
      if(_debug > 2){
        cout << *currentContext()->moduleLabel() << "status = " << ks.status() << " nactive = " << nactive << " mom = " << fseg.mom() << " chisq/dof = " << ks.chisquared()/ndof << endl;
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
	triginfo->_triggerBits.merge(TriggerFlag::track);
	// associate to the helix which triggers.  Note there may be other helices which also pass the filter
	// but filtering is by event!
	size_t index = std::distance(kscol->begin(),iks);
	triginfo->_track = art::Ptr<KalSeed>(ksH,index);
	if(_debug > 1){
	  cout << *currentContext()->moduleLabel() << " passed event " << evt.id() << endl;
	}
	break;
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool SeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      cout << *currentContext()->moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }
    return true;
  }
}
using mu2e::SeedFilter;
DEFINE_ART_MODULE(SeedFilter);
