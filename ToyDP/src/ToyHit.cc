//
// A hit in a CTracker.
// The CTracker is a toy device that is present only to illustrate the
// framework.
//
// $Id: ToyHit.cc,v 1.2 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke

#include "ToyDP/inc/ToyHit.hh"

namespace mu2e {

  ToyHit::ToyHit(): _position(0.,0.,0.), _pulseheight(0.) { }

  ToyHit::ToyHit(double x, double y, double z, double ph):
    _position(x,y,z), _pulseheight(ph) { }

  ToyHit::ToyHit( const CLHEP::Hep3Vector& pos, double ph):
    _position(pos), _pulseheight(ph) { }

  ToyHit::~ToyHit(){
  }

}

