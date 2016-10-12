//
// Track definition object
//
// $Id: TrkDefHack.cc,v 1.9 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
// Original author David Brown, LBNL
//
#include "BTrk/BaBar/BaBar.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e 
{
    // dummy variables
  HepVector TrkDefHack::_dpar(5,0);
  HepSymMatrix TrkDefHack::_dcov(5,0);
  TrkParticle TrkDefHack::_eminus(TrkParticle::e_minus);
  TrkFitDirection TrkDefHack::_downstream(TrkFitDirection::downstream);

  TrkDefHack::TrkDefHack(TimeCluster const& tclust, const HelixTraj& helix,
      TrkParticle const& tpart, TrkFitDirection const& fdir) : _shcol(0),
    _timeCluster(tclust),_h0(helix),_tpart(tpart),_fdir(fdir)
  {}

  TrkDefHack::TrkDefHack(const StrawHitCollection* shcol, std::vector<StrawHitIndex> const& hits,
      HelixTraj const& htraj,
      TrkParticle const& tpart, TrkFitDirection const& fdir) : _shcol(shcol),
  _h0(htraj), _tpart(tpart),_fdir(fdir)
  {
    setIndices(hits);
  }

  TrkDefHack::TrkDefHack(const StrawHitCollection* shcol, std::vector<StrawHitIndex> const& hits,
      TrkParticle const& tpart, TrkFitDirection const& fdir) : _shcol(shcol),
  _h0(_dpar,_dcov), _tpart(tpart),_fdir(fdir)
  {
    setIndices(hits);
  }

  TrkDefHack::TrkDefHack(TrkParticle const& tpart, TrkFitDirection const& fdir) :
   _h0(_dpar,_dcov),_tpart(tpart),_fdir(fdir)
  {}
  
  TrkDefHack::TrkDefHack(const TrkDefHack& other ) : 
    _shcol(other._shcol),
    _timeCluster(other._timeCluster),
    _h0(other._h0), _tpart(other._tpart),
    _fdir(other._fdir)
  {}
  
  TrkDefHack&
  TrkDefHack::operator = (const TrkDefHack& other) {
    if(this != &other){
      _shcol = other._shcol;
      _timeCluster = other._timeCluster;
      _h0 = other._h0;
      _tpart = other._tpart;
      _fdir = other._fdir;
    }
    return *this;
  }
    
  TrkDefHack::~TrkDefHack(){}
}
