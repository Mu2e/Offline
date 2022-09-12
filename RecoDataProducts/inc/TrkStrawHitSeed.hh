//
//  Persistent representation of a TrkStrawHit, used in the
//  persistent representation of Kalman Fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_TrkStrawHitSeed_HH
#define RecoDataProducts_TrkStrawHitSeed_HH
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include <Rtypes.h>
#include <functional>
namespace mu2e {
  struct TrkStrawHitSeed {
    // default constructor: initialization on declaration
    TrkStrawHitSeed() {}
    //KinKal constructor
    TrkStrawHitSeed(StrawHitIndex index, ComboHit const& chit,
        KinKal::ClosestApproachData const& rptca,
        KinKal::ClosestApproachData const& uptca,
        KinKal::Residual const& utresid, KinKal::Residual const& udresid,
        KinKal::Residual const& rtresid, KinKal::Residual const& rdresid,
        DriftInfo const& dinfo, WireHitState const& whs) :
      _index(index), _sid(chit.strawId()),_end(chit.driftEnd()),
      _flag(chit.flag()),
      _ambig(whs.state_), _algo(whs.algo_), _frozen(whs.frozen_),
      _edep(chit.energyDep()),_htime(chit.time()),_wdist(chit.wireDist()),_werr(chit.wireRes()), _tottdrift(chit.driftTime()),
      _ptoca(rptca.particleToca()),_stoca(rptca.sensorToca()),
      _rdoca(rptca.doca()),_rdocavar(rptca.docaVar()),
      _rdt(rptca.deltaT()), _rtocavar(rptca.tocaVar()),
      _udoca(uptca.doca()),_udocavar(uptca.docaVar()),
      _udt(uptca.deltaT()), _utocavar(uptca.tocaVar()),
      _rdrift(dinfo.driftDistance_),_rerr(dinfo.driftDistanceError_), _dvel(dinfo.driftVelocity_),_lang(dinfo.LorentzAngle_),
      _utresid(utresid.value()),_utresidmvar(utresid.measurementVariance()),_utresidpvar(utresid.parameterVariance()),
      _udresid(udresid.value()),_udresidmvar(udresid.measurementVariance()),_udresidpvar(udresid.parameterVariance()),
      _rtresid(rtresid.value()),_rtresidmvar(rtresid.measurementVariance()),_rtresidpvar(rtresid.parameterVariance()),
      _rdresid(rdresid.value()),_rdresidmvar(rdresid.measurementVariance()),_rdresidpvar(rdresid.parameterVariance()),
      _trklen(0),_hitlen(0),_stime(0.0)
    {
      // correct for end sign to return to Mu2e convention
      double endsign = 2.0*(chit.driftEnd()-0.5);
      _rupos = -endsign*rptca.sensorDirection().Dot(rptca.sensorPoca().Vect() - chit.centerPos());
      _uupos = -endsign*uptca.sensorDirection().Dot(uptca.sensorPoca().Vect() - chit.centerPos());
      _udoca *= endsign;
      _rdoca *= endsign;
      // correct flag
      _flag.merge(StrawHitFlag::track);
      if(whs.active())_flag.merge(StrawHitFlag::active);
    }

    //Legacy constructor for BTrk
    TrkStrawHitSeed(StrawHitIndex index, HitT0 const& t0, Float_t trklen, Float_t hitlen, Float_t rdrift,
        Float_t stime, Float_t upos, Float_t dt,
        Float_t wdoca, Int_t ambig, Float_t rerr, StrawHitFlag const& flag, ComboHit const& chit) :
      _index(index), _sid(chit.strawId()),_end(chit.driftEnd()),
      _flag(flag), _ambig(ambig), _algo(-10), _frozen(false),
      _edep(chit.energyDep()),_htime(chit.time()),_wdist(chit.wireDist()), _werr(chit.wireRes()),
      _tottdrift(chit.driftTime()),
      _ptoca(t0._t0),_stoca(chit.time()-stime),
      _rdoca(wdoca), _rdocavar(rerr*rerr), _rdt(dt), _rtocavar(t0._t0err*t0._t0err), _udoca(wdoca), _udocavar(rerr*rerr), _udt(dt), _utocavar(t0._t0err*t0._t0err),
      _rupos(upos),_uupos(upos),
      _rdrift(rdrift), _rerr(rerr), _dvel(0), _lang(0),
      _t0(t0), _trklen(trklen), _hitlen(hitlen),  _stime(stime){}

    // legacy interface
    auto index() const { return _index; }
    auto const&  strawId() const { return _sid; }
    auto const& flag() const { return _flag; }
    auto const& algorithm() const { return _algo; }
    auto strawHitState() const { return _ambig; }
    auto hitTime() const { return _htime; }
    auto energyDep() const { return _edep; }
    auto const& driftEnd() const { return _end; }
    auto wireDist() const { return _wdist; }
    auto wireRes() const { return _werr; }
    auto TOTDriftTime() const { return _tottdrift; }
    auto particleTOCA() const { return _ptoca; }
    auto sensorTOCA() const { return _stoca; }
    auto fitDOCA() const { return _udoca; }
    auto fitDOCAVar() const { return _udocavar; }
    auto fitDt() const { return _udt; }
    auto fitTOCAVar() const { return _utocavar; }
    auto refDOCA() const { return _rdoca; }
    auto refDOCAVar() const { return _rdocavar; }
    auto refDt() const { return _rdt; }
    auto reTOCAVar() const { return _rtocavar; }
    auto refPOCA_Upos() const { return _rupos; }
    auto driftRadius() const { return _rdrift; }
    auto radialErr() const { return _rerr; }
    HitT0 const&  t0() const { return _t0; }
    Float_t trkLen() const { return _trklen; }
    Float_t hitLen() const { return _hitlen; }
    Float_t signalTime() const { return _stime; }
    Float_t wireDOCA() const { return _rdoca; }
    Int_t ambig() const { return _ambig; }
    //
    //  Payload
    //
    StrawHitIndex   _index =0;       // index to the original straw (Combo) hit, and (for MC) MCDigi
    StrawId         _sid;   // which straw has the hit
    StrawEnd        _end;         // straw end used for hit time measurement
    StrawHitFlag    _flag;    // flag describing the status of this hit (active, ....)
    Int_t           _ambig =0;   // hit state, including LR ambiguity
    Int_t           _algo =0;     // hit updater algorithm
    Bool_t          _frozen =0; // hit state was frozen
    Float_t         _edep =0;        // reco energy deposition
    Float_t         _htime =0;   // raw hit time
    Float_t         _wdist =0;       // raw hit U position
    Float_t         _werr =0;    // raw hit U position error estimate
    Float_t         _tottdrift =0;   // drift time from TOT for this hit
    float_t         _ptoca =0;    // reference particle time of closest approach (TOCA)
    float_t         _stoca =0;    // reference sensor time of closest approach (TOCA)
    Float_t         _rdoca, _rdocavar =0;   // reference (biased) DOCA from the track to the wire, signed by the angular momentum WRT the wire and the measurement end (and variance)
    Float_t         _rdt, _rtocavar =0;   // reference (biased) time difference (and variance) at POCA
    Float_t         _udoca, _udocavar =0; // unbiaed DOCA (and variance)
    Float_t         _udt, _utocavar =0;   //unbiased dt and variance
    Float_t         _rupos =0; // reference POCA position along the straw WRT the straw middle
    Float_t         _uupos =0; // unbiased POCA position along the straw WRT the straw middle
    Float_t         _rdrift =0;  // drift radius for this hit
    Float_t         _rerr =0;    // intrinsic radial error
    Float_t         _dvel =0;  // instantaneous drift velocity
    Float_t         _lang =0; // Lorentz angle for EXB effects
    Float_t         _utresid=0, _utresidmvar=0, _utresidpvar =0; // unbiased time residual and associated measurement and parameter variances
    Float_t         _udresid=0, _udresidmvar=0, _udresidpvar =0; // unbiased distance residual and associated measurement and parameter variances
    Float_t         _rtresid=0, _rtresidmvar=0, _rtresidpvar =0; // reference time residual and associated measurement and parameter variances
    Float_t         _rdresid=0, _rdresidmvar=0, _rdresidpvar =0; // reference distance residual and associated measurement and parameter variances

    // BTrk legacy payload
    HitT0       _t0;     // time origin for this hit = track t0 + particle propagation time to this straw
    float_t     _trklen =0;    // track flightlength
    float_t     _hitlen =0;    // hit flightlength
    Float_t     _stime =0;   // signal propagation time for this hit, to the nearest end
  };
  // binary functor to sort TrkStrawHits by StrawHit index
  struct indexcompseed : public std::binary_function<TrkStrawHitSeed,TrkStrawHitSeed, bool> {
    bool operator()(const TrkStrawHitSeed& x,const TrkStrawHitSeed& y) { return x.index() < y.index(); }
  };
}
#endif
