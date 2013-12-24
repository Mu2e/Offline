//
// Original author David Norvil Brown

#include <string>
#include <sstream>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxesMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxes.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldCendBoxes>  ExtNeutShieldCendBoxesMaker::make(const SimpleConfig& c,
                                                                             double solenoidOffsetInX )
  {

    const int nB = c.getInt("ExtNeutShieldCendBoxes.numberOfBoxes");

    std::unique_ptr<ExtNeutShieldCendBoxes> shield ( new ExtNeutShieldCendBoxes );
    shield->_dimensions     .reserve(nB);
    shield->_materialNames  .reserve(nB);
    shield->_centerPositions.reserve(nB);
    shield->_hasHole        .reserve(nB);
    shield->_holeIndexes    .reserve(nB);
    shield->_holeLocations  .reserve(nB);
    shield->_holeRadius     .reserve(nB);
    shield->_holeHalfLength .reserve(nB);

    std::string baseName="ExtNeutShieldCendBoxes.b";
    int nHoles = 0;

    // Loop through all boxes
    for ( int icb = 1; icb <= nB; icb++ ) {
      // A stream for reading in data - easy to use with loop
      std::ostringstream variableStream;
      variableStream << baseName << icb;

      shield->_materialNames.push_back(c.getString(variableStream.str()+"materialName") );

      // Get the length, width, and height of each box.  Make half dim for G4
      // width is x dimension
      // height is y dimension
      // length is z dimension
      std::vector<double> tempdim;
      tempdim.push_back(c.getDouble(variableStream.str()+"width" )*CLHEP::mm/2.);
      tempdim.push_back(c.getDouble(variableStream.str()+"height")*CLHEP::mm/2.);
      tempdim.push_back(c.getDouble(variableStream.str()+"length")*CLHEP::mm/2.);
      shield->_dimensions.push_back(tempdim);

      // Location of the center of the box relative to the parent volume
      // (the detector hall).
      CLHEP::Hep3Vector loc( c.getDouble(variableStream.str()+"centerX")*CLHEP::mm,
			     c.getDouble(variableStream.str()+"centerY")*CLHEP::mm,
			     c.getDouble(variableStream.str()+"centerZ")*CLHEP::mm);

      shield->_centerPositions.push_back( loc );
      shield->_holeLocations  .push_back( loc );

      // Check on whether box has hole
      shield->_hasHole.push_back(c.getBool(variableStream.str()+"hasHole"));

      if ( shield->_hasHole.back() ) {
	// The index of this hole in the vectors of radius and length
	shield->_holeIndexes.push_back(nHoles);
        shield->_holeLocations.back().setY( 0. );
	nHoles++;

	// Read in radius and length of box.
	shield->_holeRadius    .push_back( c.getDouble(variableStream.str()+"holeRadius"    )*CLHEP::mm );
	shield->_holeHalfLength.push_back( c.getDouble(variableStream.str()+"holeHalfLength")*CLHEP::mm );

        // Translate position of hole wrt beamline
        const double deltaR   = c.getDouble( variableStream.str()+"holeDeltaR"   );
        const double deltaPhi = c.getDouble( variableStream.str()+"holeDeltaPhi" );

        if ( deltaR > 0 ) {
          // (-) in front of solenoidOffsetInX, which is positive by
          // definition, since we're in the DS portion of Mu2e
          shield->_holeLocations.back().setX( deltaR*cos(deltaPhi*CLHEP::degree) - solenoidOffsetInX );
          shield->_holeLocations.back().setY( deltaR*sin(deltaPhi*CLHEP::degree) );
        }

      } else {
	shield->_holeIndexes.push_back(-1);
      }

    }
			       
    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCendBoxes.verbosityLevel") > 0) {
      std::cout<<*shield<<std::endl;
    }

    return shield;
  }

} // namespace mu2e
