//
// Track definition object
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
  TrkDef::TrkDef(TimeCluster const& tclust, const HelixTraj& helix,
      TrkParticle const& tpart, TrkFitDirection const& fdir) :
    _timeCluster(tclust),_h0(helix),_tpart(tpart),_fdir(fdir)
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
