//
// Describe a found HoughCircle
//
// $Id: HoughCircle.cc,v 1.1 2010/04/12 22:44:51 shanahan Exp $
// $Author: shanahan $
// $Date: 2010/04/12 22:44:51 $
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

