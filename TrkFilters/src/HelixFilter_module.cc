//
//  Filter for selecting good helices (pat. rec. output): this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
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
    struct HelixCutsConfig{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool>                       configured           {     Name("configured"),              Comment("configured") };
      fhicl::OptionalAtom<bool>               requireCaloCluster   {     Name("requireCaloCluster"),      Comment("requireCaloCluster") };
      fhicl::OptionalAtom<bool>               doHelicityCheck      {     Name("doHelicityCheck"),         Comment("doHelicityCheck") };
      fhicl::OptionalAtom<int>                minNStrawHits        {     Name("minNStrawHits"),           Comment("minNStrawHits  ") };
      fhicl::OptionalAtom<double>             minHitRatio          {     Name("minHitRatio"),             Comment("minHitRatio    ") };
      fhicl::OptionalAtom<double>             minMomentum          {     Name("minMomentum"),             Comment("minMomentum    ") };
      fhicl::OptionalAtom<double>             maxMomentum          {     Name("maxMomentum"),             Comment("maxMomentum    ") };
      fhicl::OptionalAtom<double>             minPt                {     Name("minPt"),                   Comment("minPt          ") };
      fhicl::OptionalAtom<double>             maxChi2XY            {     Name("maxChi2XY"),               Comment("maxChi2XY      ") };
      fhicl::OptionalAtom<double>             maxChi2PhiZ          {     Name("maxChi2PhiZ"),             Comment("maxChi2PhiZ    ") };
      fhicl::OptionalAtom<double>             maxD0                {     Name("maxD0"),                   Comment("maxD0          ") };
      fhicl::OptionalAtom<double>             minD0                {     Name("minD0"),                   Comment("minD0          ") };
      fhicl::OptionalAtom<double>             maxAbsLambda         {     Name("maxAbsLambda"),            Comment("maxAbsLambda   ") };
      fhicl::OptionalAtom<double>             minAbsLambda         {     Name("minAbsLambda"),            Comment("minAbsLambda   ") };
      fhicl::OptionalAtom<double>             maxNLoops            {     Name("maxNLoops"),               Comment("maxNLoops      ") };
      fhicl::OptionalAtom<double>             minNLoops            {     Name("minNLoops"),               Comment("minNLoops      ") };
      fhicl::OptionalAtom<double>             slopeSigMin          {     Name("slopeSigMin"),             Comment("Minimum helix seed slope significance selection")};
      fhicl::OptionalAtom<double>             slopeSigMax          {     Name("slopeSigMax"),             Comment("Maximum helix seed slope significance selection")};
      fhicl::Sequence<std::string>            helixFitFlag         {     Name("helixFitFlag"),            Comment("helixFitFlag   "), std::vector<std::string>{"HelixOK"} };
      fhicl::OptionalAtom<bool>               prescaleUsingD0Phi   {     Name("prescaleUsingD0Phi"),      Comment("prescaleUsingD0Phi") };
      fhicl::Table<PhiPrescalingParams::Config>             prescalerPar{     Name("prescalerPar"),      Comment("prescalerPar") };
    };

    struct HelixCutsTool {
      HelixCutsTool(int helicity, const HelixCutsConfig& config):
        _configured  (config.configured()),
        _goodh       (config.helixFitFlag()){
        if (_configured){
          config.requireCaloCluster(_hascc);
          config.doHelicityCheck(_doHelicityCheck);
          _hel               = helicity;
          config.minNStrawHits(_minnstrawhits);
          config.minHitRatio(_minHitRatio);
          config.minMomentum(_minmom     );
          config.maxMomentum(_maxmom     );
          config.minPt(_minpT);
          config.maxChi2XY(_maxchi2XY);
          config.maxChi2PhiZ(_maxchi2PhiZ);
          config.maxD0(_maxd0);
          config.minD0(_mind0);
          config.maxAbsLambda(_maxlambda);
          config.minAbsLambda(_minlambda);
          config.maxNLoops(_maxnloops);
          config.minNLoops(_minnloops);
          _useSlopeSigMin = config.slopeSigMin(_slopeSigMin);
          _useSlopeSigMax = config.slopeSigMax(_slopeSigMax);
        }

        bool val;
        if (config.prescaleUsingD0Phi(val)) {
          _prescaleUsingD0Phi = val;
          if (_prescaleUsingD0Phi){
            _prescalerPar    = PhiPrescalingParams(config.prescalerPar());
          }
        }else {
          _prescaleUsingD0Phi = false;
        }
      };

      HelixCutsTool():_configured(false){};

      void setTrackerGeomHandle(const Tracker* TrackerGeom) { _myTracker = TrackerGeom; }

      int   evalIPAPresc(const float &phi0){
        int val= (_prescalerPar._amplitude - (_prescalerPar._amplitude-1)*sin(_prescalerPar._frequency*phi0 + _prescalerPar._phase));
        return val;
      }

      bool checkHelix(const HelixSeed&Helix, int NEvt, int Debug){
        if (!_configured) return true;

        //check the helicity
        if (_doHelicityCheck && !(Helix.helix().helicity() == Helicity(_hel)))  return false;
        float mm2MeV = 3./10.;//FIXME!
        HelixTool helTool(&Helix, _myTracker);
        // compute the helix momentum.  Note this is in units of mm!!!
        float hmom       = Helix.helix().momentum()*mm2MeV;
        int   nstrawhits = helTool.nstrawhits();
        float hpT        = Helix.helix().radius()*mm2MeV;
        float chi2XY     = Helix.helix().chi2dXY();
        float chi2PhiZ   = Helix.helix().chi2dZPhi();
        float d0         = Helix.helix().rcent() - Helix.helix().radius();
        float lambda     = std::fabs(Helix.helix().lambda());
        float nLoops     = helTool.nLoops();
        float hRatio     = helTool.hitRatio();
        const float slope    = Helix.recoDir().slope();
        const float slopeErr = std::fabs(Helix.recoDir().slopeErr());
        const float slopeSig = (slopeErr > 0.f) ? slope/slopeErr : 0.f;

        if(Debug > 2){
          std::cout << "[HelixFilter] : status = " << Helix.status() << " nhits = " << nstrawhits << " mom = " << hmom << std::endl;
          std::cout << "[HelixFilter] : chi2XY = " << chi2XY << " chi2ZPHI = " << chi2PhiZ << " d0 = " << d0 << " lambda = "<< lambda << " nLoops = " << nLoops << " hRatio = "<< hRatio << std::endl;
        }
        if( Helix.status().hasAllProperties(_goodh)      &&
            (!_hascc || Helix.caloCluster().isNonnull()) &&
            nstrawhits >= _minnstrawhits &&
            hpT        >= _minpT         &&
            chi2XY     <= _maxchi2XY     &&
            chi2PhiZ   <= _maxchi2PhiZ   &&
            d0         <= _maxd0         &&
            d0         >= _mind0         &&
            lambda     <= _maxlambda     &&
            lambda     >= _minlambda     &&
            nLoops     <= _maxnloops     &&
            nLoops     >= _minnloops     &&
            hmom       >= _minmom        &&
            hmom       <= _maxmom        &&
            (!_useSlopeSigMin || slopeSig > _slopeSigMin) &&
            (!_useSlopeSigMax || slopeSig < _slopeSigMax) &&
            hRatio     >= _minHitRatio ) {
          //now check if we want to prescake or not
          if (_prescaleUsingD0Phi) {
            float phiAtD0   = Helix.helix().fcent();
            int   prescaler = evalIPAPresc(phiAtD0);
            if (NEvt % prescaler != 0) return false;
          }
          return true;
        }
        return false;
      }
      bool          _configured;
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
      bool          _useSlopeSigMin;
      double        _slopeSigMin;
      bool          _useSlopeSigMax;
      double        _slopeSigMax;
      TrkFitFlag    _goodh; // helix fit flag
      bool          _prescaleUsingD0Phi;
      PhiPrescalingParams     _prescalerPar;
      const Tracker*_myTracker;
    };

    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>      helixSeedCollection  { Name("helixSeedCollection"),     Comment("helixSeedCollection") };
      fhicl::Table<HelixCutsConfig>   posHelicitySelection { Name("posHelicitySelection"),    Comment("cuts for the helices with positive helicity (like the conversion electron)")};
      fhicl::Table<HelixCutsConfig>   negHelicitySelection { Name("negHelicitySelection"),    Comment("cuts for the helices with negative helicity (like the conversion positron)")};
      fhicl::Atom<int>                debugLevel           { Name("debugLevel"),              Comment("debugLevel")     , 0 };
      fhicl::Atom<bool>               noFilter             { Name("noFilter"),                Comment("don't apply the filter decision"), false};
      fhicl::Atom<unsigned>           minNHelices          { Name("minNHelices"),             Comment("minimum number of helices passing the cuts"), 1};
      fhicl::Atom<float>              maxDt0               { Name("maxDt0"),                  Comment("maximum difference in time between helices"), -1};
    };


    using Parameters = art::EDFilter::Table<Config>;

    explicit HelixFilter(const Parameters& conf);
    virtual bool filter(art::Event& event) override;
    virtual void beginJob();
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag _hsTag;
    HelixCutsTool _posHelCuts;
    HelixCutsTool _negHelCuts;
    double        _bz0;
    const Tracker* _tracker;
    std::string   _trigPath;
    int           _debug;
    // counters
    unsigned      _nevt, _npass;
    bool          _noFilter;
    unsigned      _minNHelices;
    float         _maxDt0;

    bool  checkHelixFromHelicity(const HelixSeed&helix);
  };

  HelixFilter::HelixFilter(const Parameters& config):
    art::EDFilter{config},
    _hsTag             (config().helixSeedCollection()),
    _posHelCuts        ( 1, config().posHelicitySelection()),
    _negHelCuts        (-1, config().negHelicitySelection()),
    _debug             (config().debugLevel()),
    _nevt(0),
    _npass(0),
    _noFilter          (config().noFilter()),
    _minNHelices       (config().minNHelices()),
    _maxDt0            (config().maxDt0())
    {
      produces<TriggerInfo>();
      if(_minNHelices < 2 && _maxDt0 >= 0.) throw cet::exception("BADCONFIG") << "Requested a timing difference cut of " << _maxDt0 << " ns between helices but only " << _minNHelices << " helices required";
    }

  void HelixFilter::beginJob() {
    if ( (!_posHelCuts._configured) && (!_negHelCuts._configured)) {
      std::cout << moduleDescription().moduleLabel() << " NO HELIX CUT HAS BEEN SET. IF THAT'S NOT THE DESIRED BEHAVIOUR REVIEW YOUR CONFIGUREATION!" << std::endl;
    }
    if ( (!_posHelCuts._configured)){
      std::cout << moduleDescription().moduleLabel() << " NO HELIX CUT HAS BEEN SET FOR THE HELICES WITH POSITIVE HELICITY. IF THAT'S NOT THE DESIRED BEHAVIOUR REVIEW YOUR CONFIGUREATION!" << std::endl;
    }
    if ( (!_negHelCuts._configured)){
      std::cout << moduleDescription().moduleLabel() << " NO HELIX CUT HAS BEEN SET FOR THE HELICES WITH NEGATIVE HELICITY. IF THAT'S NOT THE DESIRED BEHAVIOUR REVIEW YOUR CONFIGUREATION!" << std::endl;
    }
  }

  bool HelixFilter::beginRun(art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();
    _posHelCuts.setTrackerGeomHandle(_tracker);
    _negHelCuts.setTrackerGeomHandle(_tracker);

    return true;
  }

  bool  HelixFilter::checkHelixFromHelicity(const HelixSeed&Helix){

    if (Helix.helix().helicity() == Helicity(_posHelCuts._hel)){
      return  _posHelCuts.checkHelix(Helix, _nevt, _debug);
    } else if (Helix.helix().helicity() == Helicity(_negHelCuts._hel)){
      return  _negHelCuts.checkHelix(Helix, _nevt, _debug);
    } else {
      return false;
    }
  }

  bool HelixFilter::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    // find the collection
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    const HelixSeedCollection* hscol = hsH.product();

    if(_debug > 1){
      std::cout << moduleDescription().moduleLabel() << " Input from: " << _hsTag << " NHelices = "<< hscol->size() << std::endl;
    }
    // loop over the collection: if any pass the selection, pass this event
    unsigned nGoodHelices(0);
    for(auto ihs = hscol->begin();ihs != hscol->end(); ++ihs) {
      auto const& hs = *ihs;

      if( checkHelixFromHelicity(hs) ) {
        ++_npass;
        ++nGoodHelices;
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

    bool dt_range = false;
    if (_maxDt0 < 0){
      dt_range = true;
    }
    else{
      for(auto ihs = hscol->begin();ihs != hscol->end(); ++ihs) {
        auto const&  hel0 = ihs;
        float  hel0_t0 =  hel0->t0()._t0;

        for(auto jhs = std::next(ihs); jhs != hscol->end(); ++jhs) {
          auto const&  hel1 =  jhs;
          float hel1_t0 = hel1->t0()._t0;
          float dt      = hel1_t0 - hel0_t0;
          if (dt < _maxDt0 && dt > -_maxDt0) {
            dt_range = true;
            break;
          }
        }
        if (dt_range) {
          break; // Break out of the outer loop
        }
      }
    }


    evt.put(std::move(triginfo));


    if(!_noFilter) {return (nGoodHelices >= _minNHelices) && dt_range;}
    else {return true;} //filtering is turned off

  }


  bool HelixFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){      std::cout << moduleDescription().moduleLabel() << " paassed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }

}
using mu2e::HelixFilter;
DEFINE_ART_MODULE(HelixFilter)
