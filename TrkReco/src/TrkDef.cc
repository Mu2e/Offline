//
// Track definition object
//
// $Id: TrkDef.cc,v 1.9 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
// Original author David Brown, LBNL
//
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkDef.hh"
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e 
{
    // dummy variables
  HepVector TrkDef::_dpar(5,0);
  HepSymMatrix TrkDef::_dcov(5,0);
  TrkParticle TrkDef::_eminus(TrkParticle::e_minus);
  TrkFitDirection TrkDef::_downstream(TrkFitDirection::downstream);

  TrkDef::TrkDef(TimeCluster const& tclust, const HelixTraj& helix,
      TrkParticle const& tpart, TrkFitDirection const& fdir) :
    _timeCluster(tclust),_h0(helix),_tpart(tpart),_fdir(fdir)
  {}

  TrkDef::TrkDef(TrkParticle const& tpart, TrkFitDirection const& fdir) :
   _h0(_dpar,_dcov),_tpart(tpart),_fdir(fdir)
  {}
  
  TrkDef::TrkDef(const TrkDef& other ) : 
    _timeCluster(other._timeCluster),
    _h0(other._h0), _tpart(other._tpart),
    _fdir(other._fdir)
  {}
  
  TrkDef&
  TrkDef::operator = (const TrkDef& other) {
    if(this != &other){
      _timeCluster = other._timeCluster;
      _h0 = other._h0;
      _tpart = other._tpart;
      _fdir = other._fdir;
    }
    return *this;
  }
    
  TrkDef::~TrkDef(){}
}
