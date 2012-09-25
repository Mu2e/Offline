//
// $Id: WireDetail.cc,v 1.5 2012/09/25 10:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:29 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/WireDetail.hh"

using namespace std;

namespace mu2e {

  WireDetail::WireDetail( std::vector<double> & thicknesses,
                   std::vector<std::string> & materialNames,
               double halfLength
                           ):
    _materialNames(materialNames),
    _shellsThicknesses(thicknesses),
    _halfLength(halfLength)
  {
          _radius=0.;
          vector<double>::iterator ithick = thicknesses.begin();
          while(ithick!= thicknesses.end()){
                  _radius += *ithick;
                  ++ithick;
          }
  }

  WireDetail::~WireDetail (){
  }

} // namespace mu2e
