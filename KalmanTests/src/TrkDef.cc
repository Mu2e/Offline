//
// Track definition object
//
// $Id: TrkDef.cc,v 1.1 2011/06/02 00:00:05 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/02 00:00:05 $
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
  
  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<size_t>& strawhits,
    const HelixTraj& helix, double t0, double t0err) :
    _straws(strawcollection), _indices(strawhits),_h0(helix),_t0(t0),_t0err(t0err)
  {}
    
  TrkDef::TrkDef(const StrawHitCollection* strawcollection, const std::vector<size_t>& strawhits,
    const HepVector& parvec, const HepSymMatrix& covar, double t0, double t0err) :
    _straws(strawcollection),_indices(strawhits),_h0(parvec,covar),_t0(t0),_t0err(t0err)
  {}
  
  TrkDef::TrkDef() : _straws(0),_h0(_dpar,_dcov),_t0(0),_t0err(-1.0)
  {}
  
  TrkDef::TrkDef(const TrkDef& other ) : _straws(other._straws),_indices(other._indices),
    _h0(other._h0), _t0(other._t0), _t0err(other._t0err)
  {}
  
  TrkDef&
  TrkDef::operator = (const TrkDef& other) {
    if(this != &other){
      _straws = other._straws;
      _indices = other._indices;
      _h0 = other._h0;
      _t0 = other._t0;
      _t0err = other._t0err;
    }
    return *this;
  }
    
  TrkDef::~TrkDef(){}
}
