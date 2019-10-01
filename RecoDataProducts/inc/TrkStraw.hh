//
//  Persistent representation recording which straws were hit
//  by a track, and where, and how much material was traversed
//  this info is part of the persistent representation of the BTrk Kalman Fit
//  Original author: Dave Brown (LBNL) 20 Feb 2017
//
#ifndef RecoDataProducts_TrkStraw_HH
#define RecoDataProducts_TrkStraw_HH
#include "DataProducts/inc/StrawId.hh"
namespace mu2e {
  struct TrkStraw {
    TrkStraw() : _doca(-1.0), _trklen(0.0), _wirelen(0.0), _slen(0.0), _radlen(0.0)  {}
    TrkStraw(StrawId const& id, double doca, double trklen, double wirelen, double slen, double radlen, double pfrac, bool active) : _straw(id), _doca(doca), _trklen(trklen), 
    _wirelen(wirelen), _slen(slen), _radlen(radlen), _pfrac(pfrac), _active(active) {}
   
    StrawId const& straw() const { return _straw; }
    Float_t doca() const { return _doca; }
    Float_t trkLen() const { return _trklen; }
    Float_t wireLen() const { return _wirelen; }
    Float_t strawLen() const { return _slen; }
    Float_t radLen() const { return _radlen; }
    Float_t pfrac() const { return _pfrac; }
    Bool_t active() const { return _active; }
  
    StrawId   _straw; // which straw was traversed
    Float_t   _doca; // DOCA from the track to the wire of this straw
    Float_t   _trklen; // length along the track from z=0 to track-wire POCA
    Float_t   _wirelen; // length along the wire from the wire middle to the track-wire POCA
    Float_t   _slen; // path length through this straw
    Float_t   _radlen; // total radiation length of material traversed in this straw
    Float_t   _pfrac; // fractional momentum change due to energy loss
    Bool_t    _active; // was this material active in the fit?
  };
}
#endif
