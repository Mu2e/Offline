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

  std::unique_ptr<ExtNeutShieldCendBoxes>  ExtNeutShieldCendBoxesMaker::make(const 
							      SimpleConfig& c)
  {


    int nB = c.getInt("ExtNeutShieldCendBoxes.numberOfBoxes");

    std::vector<std::vector<double>> dims;
    dims.reserve(nB);
    std::vector<std::string> mats;
    mats.reserve(nB);
    std::vector<CLHEP::Hep3Vector> sites;
    sites.reserve(nB);
    std::vector<bool> holey;
    holey.reserve(nB);
    std::vector<int> holeID;
    holeID.reserve(nB);
    // For now, only one hole, so don't need to reserve extra space in vector
    std::vector<double> rads;
    std::vector<double> hlengs;

    std::string baseName="ExtNeutShieldCendBoxes.b";
    int nHoles = 0;

    // Loop through all boxes
    for ( int icb = 1; icb <= nB; icb++ ) {
      // A stream for reading in data - easy to use with loop
      std::ostringstream variableStream;

      // Get the material name for the box
      variableStream << baseName << icb << "materialName";
      std::string amat = c.getString(variableStream.str());
      mats.push_back(amat);
      variableStream.str("");

      // Get the length, width, and height of each box.  Make half dim for G4
      // width is x dimension
      // height is y dimension
      // length is z dimension
      std::vector<double> tempdim;
      variableStream << baseName << icb << "width";
      tempdim.push_back(c.getDouble(variableStream.str())*CLHEP::mm/2.);
      variableStream.str("");
      variableStream << baseName << icb << "height";
      tempdim.push_back(c.getDouble(variableStream.str())*CLHEP::mm/2.);
      variableStream.str("");
      variableStream << baseName << icb << "length";
      tempdim.push_back(c.getDouble(variableStream.str())*CLHEP::mm/2.);
      dims.push_back(tempdim);

      // Location of the center of the box relative to the parent volume
      // (the detector hall).
      variableStream.str("");
      variableStream << baseName << icb << "centerX";
      double tempX = c.getDouble(variableStream.str());
      variableStream.str("");
      variableStream << baseName << icb << "centerY";
      double tempY = c.getDouble(variableStream.str());
      variableStream.str("");
      variableStream << baseName << icb << "centerZ";
      double tempZ = c.getDouble(variableStream.str());

      CLHEP::Hep3Vector loc(tempX * CLHEP::mm,
			    tempY * CLHEP::mm,
			    tempZ * CLHEP::mm);
      sites.push_back(loc);

      // Check on whether box has hole
      variableStream.str("");
      variableStream << baseName << icb << "hasHole";
      bool hasHole = c.getBool(variableStream.str());
      holey.push_back(hasHole);

      if ( hasHole ) {
	// The index of this hole in the vectors of radius and length
	holeID.push_back(nHoles);
	nHoles++;

	// Read in radius and length of box.
	variableStream.str("");
	variableStream << baseName << icb << "holeRadius";
	rads.push_back(c.getDouble(variableStream.str())*CLHEP::mm);
	variableStream.str("");
	variableStream << baseName << icb << "holeHalfLength";
	hlengs.push_back(c.getDouble(variableStream.str())*CLHEP::mm);

      } else {
	holeID.push_back(-1);
      }

    }
			       
    std::unique_ptr<ExtNeutShieldCendBoxes> res(new ExtNeutShieldCendBoxes(
								     dims,
								     mats,
								     sites,
								     holey,
								     holeID,
								     rads,
								     hlengs)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCendBoxes.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
