//
//  Filter for selecting good helices (pat. rec. output): this is part of the track trigger
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
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
// data
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
// mu2e
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
// helper function
#include "Offline/GeneralUtilities/inc/PhiPrescalingParams.hh"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
//#include "TrkFilters/inc/TrkFiltersHelpers.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

namespace mu2e
{
  class HelixFilter : public art::EDFilter
  {
    public:
      struct Config{
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>      helixSeedCollection  {     Name("helixSeedCollection"),     Comment("helixSeedCollection") };
        fhicl::Atom<bool>               requireCaloCluster   {     Name("requireCaloCluster"),      Comment("requireCaloCluster") };
        fhicl::Atom<bool>               doHelicityCheck      {     Name("doHelicityCheck"),         Comment("doHelicityCheck") };
        fhicl::Atom<int>                helicity             {     Name("helicity"),                Comment("helicity       ") };
        fhicl::Atom<int>                minNStrawHits        {     Name("minNStrawHits"),           Comment("minNStrawHits  ") };
        fhicl::Atom<double>             minHitRatio          {     Name("minHitRatio"),             Comment("minHitRatio    ") };
        fhicl::Atom<double>             minMomentum          {     Name("minMomentum"),             Comment("minMomentum    ") };
        fhicl::Atom<double>             maxMomentum          {     Name("maxMomentum"),             Comment("maxMomentum    ") };
        fhicl::Atom<double>             minPt                {     Name("minPt"),                   Comment("minPt          ") };
        fhicl::Atom<double>             maxChi2XY            {     Name("maxChi2XY"),               Comment("maxChi2XY      ") };
        fhicl::Atom<double>             maxChi2PhiZ          {     Name("maxChi2PhiZ"),             Comment("maxChi2PhiZ    ") };
        fhicl::Atom<double>             maxD0                {     Name("maxD0"),                   Comment("maxD0          ") };
        fhicl::Atom<double>             minD0                {     Name("minD0"),                   Comment("minD0          ") };
        fhicl::Atom<double>             maxAbsLambda         {     Name("maxAbsLambda"),            Comment("maxAbsLambda   ") };
        fhicl::Atom<double>             minAbsLambda         {     Name("minAbsLambda"),            Comment("minAbsLambda   ") };
        fhicl::Atom<double>             maxNLoops            {     Name("maxNLoops"),               Comment("maxNLoops      ") };
        fhicl::Atom<double>             minNLoops            {     Name("minNLoops"),               Comment("minNLoops      ") };
        fhicl::Sequence<std::string>    helixFitFlag         {     Name("helixFitFlag"),            Comment("helixFitFlag   "), std::vector<std::string>{"HelixOK"} };
        fhicl::Atom<int>                debugLevel           {     Name("debugLevel"),              Comment("debugLevel")     , 0 };
        fhicl::OptionalAtom<bool>       prescaleUsingD0Phi   {     Name("prescaleUsingD0Phi"),      Comment("prescaleUsingD0Phi") };
        fhicl::Table<PhiPrescalingParams::Config>             prescalerPar{     Name("prescalerPar"),      Comment("prescalerPar") };
      };

      using Parameters = art::EDFilter::Table<Config>;

      explicit HelixFilter(const Parameters& conf);
      virtual bool filter(art::Event& event) override;
      virtual bool beginRun(art::Run&   run   );
      virtual bool endRun( art::Run& run ) override;

    private:
      art::InputTag _hsTag;
      bool          _hascc; // Calo Cluster
      bool          _doHelicityCheck;
      int           _hel;
      int           _minnstrawhits;
      double        _minHitRatio;
      double        _minmom, _maxmom;
      double        _maxpT;
      double        _minpT;
      double        _maxchi2XY;
      double        _maxchi2PhiZ;
      double        _maxd0;
      double        _mind0;
      double        _maxlambda;
      double        _minlambda;
      double        _maxnloops;
      double        _minnloops;
      double        _bz0;
      const Tracker* _tracker;
      TrkFitFlag    _goodh; // helix fit flag
      std::string   _trigPath;
      bool          _prescaleUsingD0Phi;
      PhiPrescalingParams     _prescalerPar;
      int           _debug;
      // counters
      unsigned      _nevt, _npass;

      int evalIPAPresc(const float &phi0);
  };

  HelixFilter::HelixFilter(const Parameters& config):
    art::EDFilter{config},
    _hsTag             (config().helixSeedCollection()),
    _hascc             (config().requireCaloCluster()),
    _doHelicityCheck   (config().doHelicityCheck()),
    _hel               (config().helicity()),
    _minnstrawhits     (config().minNStrawHits()),
    _minHitRatio       (config().minHitRatio()),
    _minmom            (config().minMomentum()),
    _maxmom            (config().maxMomentum()),
    _minpT             (config().minPt()),
    _maxchi2XY         (config().maxChi2XY()),
    _maxchi2PhiZ       (config().maxChi2PhiZ()),
    _maxd0             (config().maxD0()),
    _mind0             (config().minD0()),
    _maxlambda         (config().maxAbsLambda()),
    _minlambda         (config().minAbsLambda()),
    _maxnloops         (config().maxNLoops()),
    _minnloops         (config().minNLoops()),
    _goodh             (config().helixFitFlag()),
    _debug             (config().debugLevel()),
    _nevt(0), _npass(0)
    {
      bool val;
      if (config().prescaleUsingD0Phi(val)) {
        _prescaleUsingD0Phi = val;
        if (_prescaleUsingD0Phi){
          _prescalerPar    = PhiPrescalingParams(config().prescalerPar());
        }
      }else {
        _prescaleUsingD0Phi = false;
      }
      produces<TriggerInfo>();
    }

  bool HelixFilter::beginRun(art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();
    return true;
  }

  int HelixFilter::evalIPAPresc(const float &phi0){
    //function defined by M. Whalen (m.whalen@yale.edu)
    // reference: docdb-xxxx
    int val= (_prescalerPar._amplitude - (_prescalerPar._amplitude-1)*sin(_prescalerPar._frequency*phi0 + _prescalerPar._phase));
    return val;
  }

  bool HelixFilter::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    const HelixSeedCollection* hscol = hsH.product();
    float mm2MeV = 3./10.*_bz0;
    // loop over the collection: if any pass the selection, pass this event
    for(auto ihs = hscol->begin();ihs != hscol->end(); ++ihs) {
      auto const& hs = *ihs;

      //check the helicity
      if (_doHelicityCheck && !(hs.helix().helicity() == Helicity(_hel)))        continue;

      HelixTool helTool(&hs, _tracker);
      // compute the helix momentum.  Note this is in units of mm!!!
      float hmom       = hs.helix().momentum()*mm2MeV;
      int   nstrawhits = helTool.nstrawhits();
      float hpT        = hs.helix().radius()*mm2MeV;
      float chi2XY     = hs.helix().chi2dXY();
      float chi2PhiZ   = hs.helix().chi2dZPhi();
      float d0         = hs.helix().rcent() - hs.helix().radius();
      float lambda     = std::fabs(hs.helix().lambda());
      float nLoops     = helTool.nLoops();
      float hRatio     = helTool.hitRatio();

      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << ": status = " << hs.status() << " nhits = " << nstrawhits << " mom = " << hmom << std::endl;
        std::cout << moduleDescription().moduleLabel() << ": chi2XY = " << chi2XY << " chi2ZPHI = " << chi2PhiZ << " d0 = " << d0 << " lambda = "<< lambda << " nLoops = " << nLoops << " hRatio = "<< hRatio << std::endl;
      }
      if( hs.status().hasAllProperties(_goodh) &&
          (!_hascc || hs.caloCluster().isNonnull()) &&
          nstrawhits >= _minnstrawhits &&
          hpT        >= _minpT &&
          chi2XY     <= _maxchi2XY &&
          chi2PhiZ   <= _maxchi2PhiZ &&
          d0         <= _maxd0 &&
          d0         >= _mind0 &&
          lambda     <= _maxlambda &&
          lambda     >= _minlambda &&
          nLoops     <= _maxnloops &&
          nLoops     >= _minnloops &&
          hmom       >= _minmom    &&
          hmom       <= _maxmom    &&
          hRatio     >= _minHitRatio ) {

        //now check if we want to prescake or not
        if (_prescaleUsingD0Phi) {
          float phiAtD0   = hs.helix().fcent();
          int   prescaler = evalIPAPresc(phiAtD0);
          if (_nevt % prescaler != 0)               continue;
        }
        retval = true;
        ++_npass;
        // Fill the trigger info object
        // associate to the helix which triggers.  Note there may be other helices which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(hscol->begin(),ihs);
        triginfo->_helixes.push_back(art::Ptr<HelixSeed>(hsH,index));
        if(_debug > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
        }
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool HelixFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " paassed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }

}
using mu2e::HelixFilter;
DEFINE_ART_MODULE(HelixFilter)
