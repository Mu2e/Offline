//
//  Persistent representation of a TrkCaloHit, used in the
//  persistent representation of the BTrk Kalman Fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_TrkCaloHitSeed_HH
#define RecoDataProducts_TrkCaloHitSeed_HH
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <Rtypes.h>
namespace mu2e {
  struct TrkCaloHitSeed {
    TrkCaloHitSeed() {} // initialization on declaration

    // KinKal constructor
    TrkCaloHitSeed(art::Ptr<CaloCluster>const& cluster, StrawHitFlag const& flag,
        Float_t cdepth,
        KinKal::ClosestApproachData const& rptca,
        KinKal::ClosestApproachData const& uptca,
        KinKal::Residual const& tresid,
        KinKal::VEC3 const& tmom) : _cluster(cluster),_flag(flag),
    _cdepth(cdepth),
    _rdoca(rptca.doca()), _rptoca(rptca.particlePoca().T()), _rtocavar(rptca.tocaVar()), _rdt(rptca.deltaT()),
    _udoca(uptca.doca()), _uptoca(uptca.particlePoca().T()), _utocavar(uptca.tocaVar()), _udt(uptca.deltaT()),
    _tresid(tresid.value()), _tresidmvar(tresid.measurementVariance()), _tresidpvar(tresid.parameterVariance()),
    _cpos(rptca.sensorPoca().Vect()), _tmom(tmom), _trklen(0), _hitlen(0),  _time(0), _terr(0), _rerr(0)
    {}

    // Legacy constructor for BTrk
    TrkCaloHitSeed(HitT0 const& t0, Float_t trklen, Float_t hitlen, Float_t cdoca, Float_t rerr,
        Float_t time, Float_t terr, XYZVectorF const& cpos, XYZVectorF const& tmom, StrawHitFlag const& flag) :
      _flag(flag), _cdepth(0),
      _rdoca(cdoca), _rptoca(0), _rtocavar(0), _rdt(0),
      _udoca(cdoca), _uptoca(0), _utocavar(0), _udt(0),
      _tresid(0), _tresidmvar(0), _tresidpvar(0),
      _cpos(cpos), _tmom(tmom), _trklen(trklen), _hitlen(hitlen), _time(time), _terr(terr),
      _t0(t0), _rerr(rerr)
    {}
    // accessors
    art::Ptr<CaloCluster> const& caloCluster() const { return _cluster; }

    // legacy functions
    HitT0 const&  t0() const { return _t0; }
    Float_t        trkLen() const { return _trklen; }
    Float_t        hitLen() const { return _hitlen; }
    Float_t        time() const { return _time; }
    Float_t        timeErr() const { return _terr; }
    Float_t        transverseErr() const { return _rerr; }
    //
    // Payload
    //
    art::Ptr<CaloCluster> _cluster; // cluster this hit is based on
    StrawHitFlag          _flag;          // flag describing the status of this hit (active, ....)
    Float_t               _cdepth =0;   // depth along the particle from the disk front face
    Float_t               _rdoca =0;          // reference  DOCA from the track to the cluster axis, signed by the angular momentum WRT the wire
    Float_t               _rptoca =0;          // reference  particle TOCA at POCA
    Float_t               _rtocavar =0;      // reference  variance on TOCA
    Float_t               _rdt =0;          // reference  Delta T (=cluster TOCA - particle TOCA)
    Float_t               _udoca =0;          // unbiased  DOCA from the track to the cluster axis, signed by the angular momentum WRT the wire
    Float_t               _uptoca =0;          // unbiased  particle TOCA at POCA
    Float_t               _utocavar =0;      // unbiased  variance on TOCA
    Float_t               _udt =0;          // unbiased  Delta T at POCA
    Float_t               _tresid =0, _tresidmvar=0, _tresidpvar =0; // unbiased time residual and associated measurement and parameter variances
    XYZVectorF            _cpos;       // cluster position at POCA
    XYZVectorF            _tmom;       // track momentum at POCA
    //
    // legacy payload
    //
    Float_t               _trklen =0;          // track Length from the front face to the POCA of this hit
    Float_t               _hitlen =0;          // Length along the cluster axis to the POCA of this hit
    Float_t               _time =0;          // time of this hit, = cluster time at the SIPM
    Float_t               _terr =0;          // time error assigned to this hit
    HitT0                 _t0;          // time origin for this hit
    Float_t               _rerr =0;        // intrinsic radial error
  };
}
#endif
