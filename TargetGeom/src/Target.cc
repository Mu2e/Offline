//
// Class to represent the system of target foils.
//
// $Id: Target.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke

#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

Target::Target( const SimpleConfig& c ){

  const double z0            = c.getDouble("target.z0");
  const double deltaZ        = c.getDouble("target.deltaZ");
  const double halfThickness = c.getDouble("target.halfThickness");
  vector<double> radii;
  c.getVectorDouble("target.radii", radii);

  const double rIn(0.);

  // Tools for computing the z position of each foil.
  const int n0 = radii.size()/2;
  const double offset = ( radii.size()%2 == 1) ?
    z0 : z0 + deltaZ/2.;

  for ( vector<double>::size_type i=0;
	i<radii.size(); ++i){

    // z position of the center of the foil.
    const double z = offset + (int(i)-n0)*deltaZ;

    _foils.push_back( TargetFoil( i, 
				  Hep3Vector(0.,0.,z),
				  radii[i],
				  rIn, 
				  halfThickness
				  )
		      );
  }
}


Target::~Target(){}

}

