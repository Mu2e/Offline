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
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <functional>
namespace mu2e {
  struct TrkStrawHitSeed {
    // default constructor: initialization on declaration
    TrkStrawHitSeed() {}
    //KinKal constructor
    TrkStrawHitSeed(StrawHitIndex index, ComboHit const& chit,
        KinKal::ClosestApproachData const& rptca,
        KinKal::ClosestApproachData const& uptca,
        KinKal::Residual const& utresid, KinKal::Residual const& udresid, KinKal::Residual const& ulresid,
        KinKal::Residual const& rtresid, KinKal::Residual const& rdresid, KinKal::Residual const& rlresid,
        DriftInfo const& dinfo,
        WireHitState const& whs,
        Straw const& straw) :
      _index(index), _sid(chit.strawId()),_eend(chit.earlyEnd()),
      _flag(chit.flag()),_kkshflag(whs.flag_),
      _ambig(whs.state_), _algo(whs.algo_), _frozen(whs.frozen_),
      _bkgqual(whs.quality(WireHitState::bkg)),
      _signqual(whs.quality(WireHitState::sign)),
      _driftqual(whs.quality(WireHitState::drift)),
      _chi2qual(whs.quality(WireHitState::chi2)),
      _edep(chit.energyDep()),_etime(chit.endTimes()),_wdist(chit.wireDist()),_werr(chit.wireRes()), _tottdrift(chit.driftTime()),
      _tot(chit.TOTs()),
      _ptoca(rptca.particleToca()),_stoca(rptca.sensorToca()),
      _rdoca(rptca.doca()),_rdocavar(rptca.docaVar()),
      _rdt(rptca.deltaT()), _rtocavar(rptca.tocaVar()),
      _udoca(uptca.doca()),_udocavar(uptca.docaVar()),
      _udt(uptca.deltaT()), _utocavar(uptca.tocaVar()),
      _rdrift(dinfo.rDrift_),_cdrift(dinfo.cDrift_ ),
      _sderr(dinfo.signedDriftError_),_uderr(dinfo.unsignedDriftError_), _dvel(dinfo.driftVelocity_),_lang(dinfo.LorentzAngle_),
      _utresid(utresid.value()),_utresidmvar(utresid.measurementVariance()),_utresidpvar(utresid.parameterVariance()),
      _udresid(udresid.value()),_udresidmvar(udresid.measurementVariance()),_udresidpvar(udresid.parameterVariance()),
      _ulresid(ulresid.value()),_ulresidmvar(ulresid.measurementVariance()),_ulresidpvar(ulresid.parameterVariance()),
      _rtresid(rtresid.value()),_rtresidmvar(rtresid.measurementVariance()),_rtresidpvar(rtresid.parameterVariance()),
      _rdresid(rdresid.value()),_rdresidmvar(rdresid.measurementVariance()),_rdresidpvar(rdresid.parameterVariance()),
      _rlresid(rlresid.value()),_rlresidmvar(rlresid.measurementVariance()),_rlresidpvar(rlresid.parameterVariance()),
      _upoca(XYZVectorF(uptca.particlePoca().Vect()))
    {
      // compute position along wire according to Mu2e convention
      double endsign = chit.earlyEnd().endSign();
      _rupos = -endsign*rptca.sensorDirection().Dot(rptca.sensorPoca().Vect() - chit.centerPos());
      _uupos = -endsign*uptca.sensorDirection().Dot(uptca.sensorPoca().Vect() - chit.centerPos());
      // this signing makes the values consistent with MC truth, but needs to be accomodated when rebuilding the fit
      _udoca *= endsign;
      _rdoca *= endsign;
      // correct flag
      _flag.merge(StrawHitFlag::track);
      if(whs.active())_flag.merge(StrawHitFlag::active);

      // calculate the doca and phi relative to the straw envelope at POCA to wire
      auto ppoca = XYZVectorF(uptca.particlePoca().Vect());
      static XYZVectorF zdir(0.0,0.0,1.0); // relative to Z
      auto wmid = XYZVectorF(straw.wirePosition());
      auto wdir = XYZVectorF(straw.wireDirection());
      auto delta = ppoca - wmid;
      float dw = delta.Dot(wdir);
      XYZVectorF cperp = delta - dw*wdir; // just perp part
      auto raddir = wdir.Cross(zdir);
      if (raddir.Dot(wmid) < 0.0) raddir *= -1.0; // sign radially outwards
      _uwirephi = atan2(cperp.Dot(raddir),cperp.Dot(zdir));

      auto smid = XYZVectorF(straw.strawPosition());
      auto sdir = XYZVectorF(straw.strawDirection());
      delta = ppoca - smid; // particle poca to wire WRT straw middle
      dw = delta.Dot(sdir);
      cperp = delta - dw*sdir; // just perp part
      raddir = sdir.Cross(zdir);
      if(raddir.Dot(smid) < 0.0) raddir *= -1.0; // sign radially outwards
      _ustrawphi = atan2(cperp.Dot(raddir),cperp.Dot(zdir)); // angle around wire WRT z axis in range -pi,pi
      _ustrawdist = sqrt(cperp.mag2());

      _wdot = uptca.particleDirection().Dot(straw.wireDirection());
    }

    //Legacy constructor for BTrk
    TrkStrawHitSeed(StrawHitIndex index, HitT0 const& t0, float trklen, float hitlen, float rdrift,
        float stime, float upos, float dt,
        float wdoca, int ambig, float rerr, StrawHitFlag const& flag, ComboHit const& chit) :
      _index(index), _sid(chit.strawId()),_eend(chit.earlyEnd()),
      _flag(flag), _ambig(ambig), _algo(-10), _frozen(false),
      _edep(chit.energyDep()),_wdist(chit.wireDist()), _werr(chit.wireRes()),
      _tottdrift(chit.driftTime()), _tot{0.0,0.0},
      _ptoca(t0._t0),_stoca(chit.time()-stime),
      _rdoca(wdoca), _rdocavar(rerr*rerr), _rdt(dt), _rtocavar(t0._t0err*t0._t0err), _udoca(wdoca), _udocavar(rerr*rerr), _udt(dt), _utocavar(t0._t0err*t0._t0err),
      _rupos(upos),_uupos(upos),
      _rdrift(rdrift),_cdrift(rdrift), _sderr(rerr), _dvel(0), _lang(0),
      _ustrawdist(0), _ustrawphi(0), _uwirephi(0), _wdot(0),
      _t0(t0), _trklen(trklen), _hitlen(hitlen),  _stime(stime){
        _etime[chit.earlyEnd()] = chit.time();
      }

    // legacy interface
    auto index() const { return _index; }
    auto const&  strawId() const { return _sid; }
    auto const& flag() const { return _flag; }
    auto const& algorithm() const { return _algo; }
    auto strawHitState() const { return _ambig; }
    auto time() const { return _etime[_eend]; }
    auto energyDep() const { return _edep; }
    auto const& earlyEnd() const { return _eend; }
    StrawEnd lateEnd() const { return _eend.otherEnd(); }
    auto wireDist() const { return _wdist; }
    auto wireRes() const { return _werr; }
    auto TOTDriftTime() const { return _tottdrift; }
    auto particleToca() const { return _ptoca; }
    auto sensorToca() const { return _stoca; }
    auto fitDOCA() const { return _udoca; }
    auto fitDOCAVar() const { return _udocavar; }
    auto fitDt() const { return _udt; }
    auto fitTocaVar() const { return _utocavar; }
    auto refDOCA() const { return _rdoca; }
    auto refDOCAVar() const { return _rdocavar; }
    auto refDt() const { return _rdt; }
    auto reTocaVar() const { return _rtocavar; }
    auto refPOCA_Upos() const { return _rupos; }
    auto driftRadius() const { return _rdrift; }
    auto radialErr() const { return _sderr; }
    HitT0 const&  t0() const { return _t0; }
    float trkLen() const { return _trklen; }
    float hitLen() const { return _hitlen; }
    float signalTime() const { return _stime; }
    float wireDOCA() const { return _rdoca; }
    int ambig() const { return _ambig; }
    // return a true WireHitState
    WireHitState wireHitState() const {
      return WireHitState(static_cast<WireHitState::State>(_ambig),static_cast<StrawHitUpdaters::algorithm>(_algo),_kkshflag);
    }
    //
    //  Payload
    //
    StrawHitIndex   _index =0;       // index to the original straw (Combo) hit, and (for MC) MCDigi
    StrawId         _sid;   // which straw has the hit
    StrawEnd        _eend;         // straw end with the earliest TDC reading
    StrawHitFlag    _flag;    // flag describing the status of this hit (active, ....)
    KKSHFlag        _kkshflag; // flag from KinKal fit
    int           _ambig =0;   // hit state, including LR ambiguity
    int           _algo =0;     // hit updater algorithm
    bool          _frozen =0; // hit state was frozen
    float         _bkgqual =0; // hit background rejection quality
    float         _signqual =0; // hit LR ambiguity sign quality
    float         _driftqual =0; // hit drift quality
    float         _chi2qual =0; // hit chi2 quality
    float         _edep =0;        // reco energy deposition
    TrkTypes::TDCTimes _etime = {0,0};   // raw hit times, by end
    float         _wdist =0;       // raw hit U position
    float         _werr =0;    // raw hit U position error estimate
    float         _tottdrift =0;   // drift time from TOT for this hit
    TrkTypes::TOTTimes _tot = {0,0};   // TOT times in ns from each end
    float         _ptoca =0;    // reference particle time of closest approach (TOCA)
    float         _stoca =0;    // reference sensor time of closest approach (TOCA)
    float         _rdoca, _rdocavar =0;   // reference (biased) DOCA from the track to the wire, signed by the angular momentum WRT the wire and the measurement end (and variance)
    float         _rdt, _rtocavar =0;   // reference (biased) time difference (and variance) at POCA
    float         _udoca, _udocavar =0; // unbiaed DOCA (and variance)
    float         _udt, _utocavar =0;   //unbiased dt and variance
    float         _rupos =0; // reference POCA position along the straw WRT the straw middle
    float         _uupos =0; // unbiased POCA position along the straw WRT the straw middle
    float         _rdrift =0;  // drift radius for this hit
    float         _cdrift =0;  // cluster drift radius for this hit
    float         _sderr =0;    //signed drift radius error
    float         _uderr =0;    // unsigned drift radius error
    float         _dvel =0;  // instantaneous drift velocity
    float         _lang =0; // Lorentz angle for EXB effects
    float         _utresid=0, _utresidmvar=0, _utresidpvar =0; // unbiased time residual and associated measurement and parameter variances
    float         _udresid=0, _udresidmvar=0, _udresidpvar =0; // unbiased distance residual and associated measurement and parameter variances
    float         _ulresid=0, _ulresidmvar=0, _ulresidpvar =0; // unbiased longitudinal residual and associated measurement and parameter variances
    float         _rtresid=0, _rtresidmvar=0, _rtresidpvar =0; // reference time residual and associated measurement and parameter variances
    float         _rdresid=0, _rdresidmvar=0, _rdresidpvar =0; // reference distance residual and associated measurement and parameter variances
    float         _rlresid=0, _rlresidmvar=0, _rlresidpvar =0; // reference longitudinal residual and associated measurement and parameter variances
    float         _ustrawdist = 0;
    float         _ustrawphi = 0;
    float         _uwirephi = 0;
    float         _wdot = 0;
    XYZVectorF    _upoca;

    // BTrk legacy payload
    HitT0       _t0;     // time origin for this hit = track t0 + particle propagation time to this straw
    float     _trklen =0;    // track flightlength
    float     _hitlen =0;    // hit flightlength
    float     _stime =0;   // signal propagation time for this hit, to the nearest end
  };
  // binary functor to sort TrkStrawHits by StrawHit index
  struct indexcompseed {
    bool operator()(const TrkStrawHitSeed& x,const TrkStrawHitSeed& y) { return x.index() < y.index(); }
  };
}
#endif
