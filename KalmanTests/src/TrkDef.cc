//
// Track definition object
//
// $Id: TrkDef.cc,v 1.7 2012/07/25 20:56:57 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/25 20:56:57 $
//
// Original author David Brown, LBNL
//
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"

namespace mu2e 
{
    // dummy variables
  HepVector TrkDef::_dpar(5,0);
  HepSymMatrix TrkDef::_dcov(5,0);
  TrkParticle TrkDef::_eminus(TrkParticle::e_minus);
  TrkFitDirection TrkDef::_downstream(TrkFitDirection::downstream);

  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
    const HelixTraj& helix, TrkParticle const& tpart, TrkFitDirection const& fdir,
    double t0, double t0err) :
    _straws(strawcollection), _indices(strawhits),_h0(helix),_ptraj(0),_tpart(tpart),_fdir(fdir),_t0(t0,t0err)
  {}
    
  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
    const HepVector& parvec, const HepSymMatrix& covar,
    TrkParticle const& tpart, TrkFitDirection const& fdir,
    double t0, double t0err) :
    _straws(strawcollection),_indices(strawhits),_h0(parvec,covar),_ptraj(0),_tpart(tpart),_fdir(fdir),_t0(t0,t0err)
  {}

  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
    TrkParticle const& tpart, TrkFitDirection const& fdir) :
    _straws(strawcollection), _indices(strawhits),_h0(_dpar,_dcov),_ptraj(0),_tpart(tpart),_fdir(fdir),_t0(0.0,-1.0)
  {}

  TrkDef::TrkDef(const StrawHitCollection* strawcollection,TrkParticle const& tpart, TrkFitDirection const& fdir) :
  _straws(strawcollection),_h0(_dpar,_dcov),_ptraj(0),_tpart(tpart),_fdir(fdir),_t0(0.0,-1.0)
  {}

  TrkDef::TrkDef(TrkParticle const& tpart, TrkFitDirection const& fdir) :
  _straws(0),_h0(_dpar,_dcov),_ptraj(0),_tpart(tpart),_fdir(fdir),_t0(0.0,-1.0)
  {}
  
  TrkDef::TrkDef(const TrkDef& other ) : _straws(other._straws),_indices(other._indices),
    _h0(other._h0), _ptraj(other._ptraj), _tpart(other._tpart),
    _fdir(other._fdir), _t0(other._t0)
  {}
  
  TrkDef&
  TrkDef::operator = (const TrkDef& other) {
    if(this != &other){
      _straws = other._straws;
      _indices = other._indices;
      _h0 = other._h0;
      _ptraj = other._ptraj;
      _tpart = other._tpart;
      _fdir = other._fdir;
      _t0 = other._t0;
    }
    return *this;
  }
    
  TrkDef::~TrkDef(){}
}
