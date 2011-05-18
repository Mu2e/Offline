//
// Describe a found HoughCircle
//
// $Id: HoughCircle.cc,v 1.2 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Peter Shanahan

#include "ToyDP/inc/HoughCircle.hh"

namespace mu2e {

  HoughCircle::HoughCircle(): _center(0.,0.), _radius(0.), _nStraws(0) { }

  HoughCircle::HoughCircle(double x, double y, double radius, uint32_t nstraws):
    _center(x,y), _radius(radius), _nStraws(nstraws) { }

  HoughCircle::HoughCircle( const CLHEP::Hep2Vector& center, double radius, uint32_t nstraws):
    _center(center), _radius(radius), _nStraws(nstraws) { }

  HoughCircle::~HoughCircle(){
  }

}

