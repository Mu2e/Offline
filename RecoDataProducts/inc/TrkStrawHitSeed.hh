//
//  Persistent representation of a TrkStrawHit, used in the
//  persistent representation of the BTrk Kalman Fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_TrkStrawHitSeed_HH
#define RecoDataProducts_TrkStrawHitSeed_HH
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <Rtypes.h>
#include <functional>
namespace mu2e {
  struct TrkStrawHitSeed {
    TrkStrawHitSeed() : _index(0), _ambig(0), _algo(-1), _edep(0),  _htime(0), _wdist(0), _werr(0),
    _ptoca(0),
    _wdoca(0), _wdocavar(0), _wdt(0), _wtocavar(0),
    _doca(0.0), _docavar(0), _dt(0), _tocavar(0),
    _upos(0),
    _tresid(0), _tresidmvar(0), _tresidpvar(0),
    _dresid(0), _dresidmvar(0), _dresidpvar(0),
    _trklen(0), _hitlen(0), _rdrift(0), _dtime(0), _stime(0), _rerr(0)    {}

    //KinKal constructor
    TrkStrawHitSeed(StrawHitIndex index, ComboHit const& chit,
        KinKal::ClosestApproachData const& refptca,
        KinKal::ClosestApproachData const& fitptca,
        KinKal::Residual utresid, KinKal::Residual udresid,
        int whstate, int algo) :
      _index(index), _sid(chit.strawId()),_end(chit.driftEnd()),
      _flag(chit.flag()), _ambig(whstate), _algo(algo),
      _edep(chit.energyDep()),_htime(chit.time()),_wdist(chit.wireDist()),_werr(chit.wireRes()),
      _ptoca(refptca.particleToca()),
      _wdoca(refptca.doca()),_wdocavar(refptca.docaVar()),
      _wdt(refptca.deltaT()), _wtocavar(refptca.tocaVar()),
      _doca(fitptca.doca()),_docavar(fitptca.docaVar()),
      _dt(fitptca.deltaT()), _tocavar(fitptca.tocaVar()),
      _tresid(utresid.value()),_tresidmvar(utresid.measurementVariance()),_tresidpvar(utresid.parameterVariance()),
      _dresid(udresid.value()),_dresidmvar(udresid.measurementVariance()),_dresidpvar(udresid.parameterVariance()),
      _rdrift(0.0), _dtime(0.0), _stime(0.0), _rerr(0.0)
    {
      // correct for end sign to return to Mu2e convention
      double endsign = 2.0*(chit.driftEnd()-0.5);
      _upos = -endsign*refptca.sensorDirection().Dot(refptca.sensorPoca().Vect() - chit.centerPos());
      _doca *= endsign;
      _wdoca *= endsign;
    }

    //Legacy constructor for BTrk
    TrkStrawHitSeed(StrawHitIndex index, HitT0 const& t0, Float_t trklen, Float_t hitlen, Float_t rdrift,
        Float_t stime, Float_t upos, Float_t dt,
        Float_t wdoca, Int_t ambig, Float_t rerr, StrawHitFlag const& flag, ComboHit const& chit) :
      _index(index), _sid(chit.strawId()),_end(chit.driftEnd()),
      _flag(flag), _ambig(ambig), _algo(-1),
      _edep(chit.energyDep()),_htime(chit.time()),_wdist(chit.wireDist()), _werr(chit.wireRes()), _ptoca(t0._t0),
       _wdoca(wdoca), _wdocavar(rerr*rerr), _wdt(dt), _wtocavar(t0._t0err*t0._t0err), _doca(wdoca), _docavar(rerr*rerr), _dt(dt), _tocavar(t0._t0err*t0._t0err),_upos(upos),
      _t0(t0), _trklen(trklen), _hitlen(hitlen), _rdrift(rdrift), _dtime(chit.driftTime()), _stime(stime), _rerr(rerr){}

    // accessors
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
    auto particleTOCA() const { return _ptoca; }
    auto fitDOCA() const { return _doca; }
    auto fitDOCAVar() const { return _docavar; }
    auto fitDt() const { return _dt; }
    auto fitTOCAVar() const { return _tocavar; }
    auto refDOCA() const { return _wdoca; }
    auto refDOCAVar() const { return _wdocavar; }
    auto refDt() const { return _wdt; }
    auto reTOCAVar() const { return _wtocavar; }
    auto refPOCA_Upos() const { return _upos; }
    auto timeResidual() const { return _tresid; }
    auto timeResidMeasurementVariance() const { return _tresidmvar; }
    auto timeResidParameterVariance() const { return _tresidpvar; }
    auto distResidual() const { return _dresid; }
    auto distResidMeasurementVariance() const { return _dresidmvar; }
    auto distResidParameterVariance() const { return _dresidpvar; }

    // legacy BTrk interface
    HitT0 const&  t0() const { return _t0; }
    Float_t trkLen() const { return _trklen; }
    Float_t hitLen() const { return _hitlen; }
    Float_t driftRadius() const { return _rdrift; }
    Float_t signalTime() const { return _stime; }
    Float_t TOTDriftTime() const { return _dtime; }
    Float_t radialErr() const { return _rerr; }
    Float_t wireDOCA() const { return _wdoca; }
    Int_t ambig() const { return _ambig; }
    //
    //  Payload
    //
    StrawHitIndex   _index;       // index to the original straw (Combo) hit, and (for MC) MCDigi
    StrawId         _sid;   // which straw has the hit
    StrawEnd        _end;         // straw end used for hit time measurement
    StrawHitFlag    _flag;    // flag describing the status of this hit (active, ....)
    Int_t           _ambig;   // hit state, including LR ambiguity
    Int_t           _algo;     // hit updater algorithm
    Float_t         _edep;        // reco energy deposition
    Float_t         _htime;   // raw hit time
    Float_t         _wdist;       // raw hit U position
    Float_t         _werr;    // raw hit U position error estimate
    float_t         _ptoca;    // reference particle time of closest approach (TOCA)
    Float_t         _wdoca, _wdocavar;   // reference (biased) DOCA from the track to the wire, signed by the angular momentum WRT the wire and the measurement end (and variance)
    Float_t         _wdt, _wtocavar;   // reference (biased) time difference (and variance) at POCA
    Float_t         _doca, _docavar; // unbiaed DOCA (and variance)
    Float_t         _dt, _tocavar;   // fit (unbiased) dt and variance
    Float_t         _upos; // POCA position along the straw WRT the straw middle
    Float_t         _tresid, _tresidmvar, _tresidpvar; // unbiased time residual and associated measurement and parameter variances
    Float_t         _dresid, _dresidmvar, _dresidpvar; // unbiased distance residual and associated measurement and parameter variances

    // BTrk legacy payload
    HitT0       _t0;     // time origin for this hit = track t0 + particle propagation time to this straw
    float_t     _trklen;    // track flightlength
    float_t     _hitlen;    // hit flightlength
    Float_t     _rdrift;  // drift radius for this hit
    Float_t     _dtime;   // drift time from TOT for this hit
    Float_t     _stime;   // signal propagation time for this hit, to the nearest end
    Float_t     _rerr;    // intrinsic radial error
  };
  // binary functor to sort TrkStrawHits by StrawHit index
  struct indexcompseed : public std::binary_function<TrkStrawHitSeed,TrkStrawHitSeed, bool> {
    bool operator()(const TrkStrawHitSeed& x,const TrkStrawHitSeed& y) { return x.index() < y.index(); }
  };
}
#endif
