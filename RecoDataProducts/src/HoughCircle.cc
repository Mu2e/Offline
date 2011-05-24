//
// Describe a found HoughCircle
//
// $Id: HoughCircle.cc,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Peter Shanahan

#include "RecoDataProducts/inc/HoughCircle.hh"

namespace mu2e {

  HoughCircle::HoughCircle(): _center(0.,0.), _radius(0.), _nStraws(0) { }

  HoughCircle::HoughCircle(double x, double y, double radius, unsigned nstraws):
    _center(x,y), _radius(radius), _nStraws(nstraws) { }

  HoughCircle::HoughCircle( const CLHEP::Hep2Vector& center, double radius, unsigned nstraws):
    _center(center), _radius(radius), _nStraws(nstraws) { }

  HoughCircle::~HoughCircle(){
  }

}

