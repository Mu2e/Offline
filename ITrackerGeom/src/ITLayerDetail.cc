//
// $Id: ITLayerDetail.cc,v 1.4 2012/09/25 10:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:29 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/ITLayerDetail.hh"

using namespace std;

namespace mu2e {

  ITLayerDetail::ITLayerDetail():
          _center_radius_ringIn(0.0),
          _center_radius_ringOut(0.0),
          _epsilonIn(0.0),
          _epsilonOut(0.0),
          _halfLength(0.0),
          _materialNames("")
  {
  }

  ITLayerDetail::ITLayerDetail(  double center_radius_ringIn, double center_radius_ringOut, double epsilonIn,
                     double epsilonOut, double halfLength, std::string materialNames
                     ):
     _center_radius_ringIn(center_radius_ringIn),
     _center_radius_ringOut(center_radius_ringOut),
    _epsilonIn(epsilonIn),
    _epsilonOut(epsilonOut),
    _halfLength(halfLength),
    _materialNames(materialNames)
  {
  }

  ITLayerDetail::~ITLayerDetail (){
  }

} // namespace mu2e
