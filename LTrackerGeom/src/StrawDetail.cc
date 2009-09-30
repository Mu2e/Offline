//
// Details common to many straws.
//
//
// $Id: StrawDetail.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "LTrackerGeom/inc/StrawDetail.hh"

using namespace std;

namespace mu2e {

  StrawDetail::StrawDetail( vector<string> const& materialNames,
			    double radius,
			    double thickness,
			    double halfLength,
			    double rwire
			   ):
    _materialNames(materialNames),
    _radius(radius),
    _thickness(thickness),
    _halfLength(halfLength),
    _rwire(rwire)
  {
  }
  
  StrawDetail::~StrawDetail (){
  }

} // namespace mu2e
