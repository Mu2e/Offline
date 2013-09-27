//
// Original author David Norvil Brown

#include <string>
#include <sstream>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCryoBoxesMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCryoBoxes.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldCryoBoxes>  ExtNeutShieldCryoBoxesMaker::make(const 
							      SimpleConfig& c)
  {

    int nB = c.getInt("ExtNeutShieldCryoBoxes.numberOfBoxes");

    std::vector<std::vector<double>> dims;
    dims.reserve(nB);
    std::vector<std::string> mats;
    mats.reserve(nB);
    std::vector<CLHEP::Hep3Vector> sites;
    sites.reserve(nB);

    std::string baseName="ExtNeutShieldCryoBoxes.b";
    std::string amat;
    for ( int icb = 1; icb <= nB; icb++ ) {
      // A stream for reading in data - easy to use with loop
      std::ostringstream variableStream;

      // Get the material name for the box
      variableStream << baseName << icb << "materialName";
      std::string amat = c.getString(variableStream.str());
      mats.push_back(amat);
      variableStream.str("");

    // Get the length, width, and height of each box.  Make half dim for G4
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
    }



			   
    
    std::unique_ptr<ExtNeutShieldCryoBoxes> res(new ExtNeutShieldCryoBoxes(
								     dims,
								     mats,
								     sites)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCryoBoxes.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
