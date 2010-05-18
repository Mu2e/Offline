#include "ITrackerGeom/inc/WireDetail.hh"

using namespace std;

namespace mu2e {

  WireDetail::WireDetail( std::vector<double> & thicknesses,
                   std::vector<std::string> & materialNames,
               double halfLength
                           ):
    _shellsThicknesses(thicknesses),
    _materialNames(materialNames),
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
