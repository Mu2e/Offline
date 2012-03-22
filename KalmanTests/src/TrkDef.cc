//
// Track definition object
//
// $Id: TrkDef.cc,v 1.6 2012/03/22 22:30:44 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/03/22 22:30:44 $
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
  
  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
    const HelixTraj& helix, double t0, double t0err) :
    _straws(strawcollection), _indices(strawhits),_h0(helix),_ptraj(0),_t0(t0,t0err)
  {}
    
  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
    const HepVector& parvec, const HepSymMatrix& covar, double t0, double t0err) :
    _straws(strawcollection),_indices(strawhits),_h0(parvec,covar),_ptraj(0),_t0(t0,t0err)
  {}

  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits) :
  _straws(strawcollection), _indices(strawhits),_h0(_dpar,_dcov),_ptraj(0),_t0(0.0,-1.0)
  {}

  TrkDef::TrkDef(const StrawHitCollection* strawcollection) :
  _straws(strawcollection),_h0(_dpar,_dcov),_ptraj(0),_t0(0.0,-1.0)
  {}
  
  TrkDef::TrkDef() : _straws(0),_h0(_dpar,_dcov),_ptraj(0),_t0(0.0,-1.0)
  {}
  
  TrkDef::TrkDef(const TrkDef& other ) : _straws(other._straws),_indices(other._indices),
    _h0(other._h0), _ptraj(other._ptraj), _t0(other._t0)
  {}
  
  TrkDef&
  TrkDef::operator = (const TrkDef& other) {
    if(this != &other){
      _straws = other._straws;
      _indices = other._indices;
      _h0 = other._h0;
      _ptraj = other._ptraj;
      _t0 = other._t0;
    }
    return *this;
  }
    
  TrkDef::~TrkDef(){}
}
