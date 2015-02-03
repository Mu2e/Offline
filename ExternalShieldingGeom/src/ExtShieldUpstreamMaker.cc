//
// Original author David Norvil Brown, University of Louisville
//
// This file is the implementation of the ExtShieldUpstreamMaker class.
// It builds the ExtShieldUpstream using SimpleConfig.

#include <string>
#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalShieldingGeom/inc/ExtShieldUpstreamMaker.hh"
#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtShieldUpstream>  ExtShieldUpstreamMaker::make(const 
							      SimpleConfig& c)
  {

    // We use just a few basic types, so simplify config file by
    // only specifying dimensions, tolerances, and material once per type.
    // Locations and orientations (?) will have to be entered for each box.

    const int nType = c.getInt("ExtShieldUpstream.numberOfBoxTypes");

    std::vector<int>                       nBoxesOfType;
    std::vector<std::vector<double> >      dimsOfType;
    std::vector<std::vector<double> >      tolsOfType;
    std::vector<std::string>               materialOfType;
    nBoxesOfType                           .reserve(nType);
    dimsOfType                             .reserve(nType);
    tolsOfType                             .reserve(nType);
    materialOfType                         .reserve(nType);

    std::string typeBaseName     = "ExtShieldUpstream.nBoxType";
    std::string dimsBaseName     = "ExtShieldUpstream.dims";
    std::string tolsBaseName     = "ExtShieldUpstream.tols";
    std::string materialBaseName = "ExtShieldUpstream.materialType";
    std::string centerBaseName   = "ExtShieldUpstream.center";
    std::string orientBaseName   = "ExtShieldUpstream.orientation";

    // Some temporary holders
    std::vector<double> tempDoubleVec;

    // Loop over the various box types and fill the vectors with 
    // information needed for the types.
    for ( int iType = 1; iType <= nType; iType++ ) {
      // A stream trick for creating variable variable names!
      std::ostringstream bTypeNumberVarName;
      bTypeNumberVarName << typeBaseName << iType;
      nBoxesOfType.push_back(c.getInt(bTypeNumberVarName.str()));

      // Get the dimensions for this type box - u,v,w are the self-coordinates
      std::ostringstream bdimsUVarName;
      bdimsUVarName << dimsBaseName << "UType" << iType;
      std::ostringstream bdimsVVarName;
      bdimsVVarName << dimsBaseName << "VType" << iType;
      std::ostringstream bdimsWVarName;
      bdimsWVarName << dimsBaseName << "WType" << iType;
      // Divide dimensions by 2.0 because G4 boxes are created in terms of 
      // half-lengths.
      tempDoubleVec.push_back(c.getDouble(bdimsUVarName.str())*CLHEP::mm/2.0);
      tempDoubleVec.push_back(c.getDouble(bdimsVVarName.str())*CLHEP::mm/2.0);
      tempDoubleVec.push_back(c.getDouble(bdimsWVarName.str())*CLHEP::mm/2.0);
      dimsOfType.push_back(tempDoubleVec);
      tempDoubleVec.clear();  // So it can be re-used

      // Get the tolerances for this type box - du,dv,dw
      std::ostringstream btolsUVarName;
      btolsUVarName << tolsBaseName << "UType" << iType;
      std::ostringstream btolsVVarName;
      btolsVVarName << tolsBaseName << "VType" << iType;
      std::ostringstream btolsWVarName;
      btolsWVarName << tolsBaseName << "WType" << iType;
      double du = c.getDouble(btolsUVarName.str());
      double dv = c.getDouble(btolsVVarName.str());
      double dw = c.getDouble(btolsWVarName.str());
      if ( du > 10.0 || dv > 10.0 || dw > 10.0 || du < -20.0 || dv < -20.0 || dw < -20.0 ) {
	// Throw if tolerances out of tolerance.
	throw cet::exception("GEOM")
	  << "Tolerance on ExternalShielding Upstream outside limits. "
	  << "\nTolerances must be between -20 and 10 mm.";
      }
      tempDoubleVec.push_back(du*CLHEP::mm/2.0);
      tempDoubleVec.push_back(dv*CLHEP::mm/2.0);
      tempDoubleVec.push_back(dw*CLHEP::mm/2.0);
      tolsOfType.push_back(tempDoubleVec);
      tempDoubleVec.clear();  // So it can be re-used

      // Get the material name for this type box
      std::ostringstream bmatVarName;
      bmatVarName << materialBaseName << iType;
      materialOfType.push_back( c.getString(bmatVarName.str()));

    }

    // Get total number of boxes from info collected above
    int nBoxesTot = 0;
    for ( unsigned int ib = 0; ib < nBoxesOfType.size(); ib++ ) {
      nBoxesTot += nBoxesOfType[ib];
    }

    // Set up the actual vectors used to build the upstream shielding
    std::vector<std::vector<double> > dims;
    std::vector<std::vector<double> > tols;
    std::vector<std::string>          mats;
    std::vector<CLHEP::Hep3Vector>    sites;
    std::vector<std::string>          orients;
    dims.reserve(nBoxesTot);
    tols.reserve(nBoxesTot);
    mats.reserve(nBoxesTot);
    sites.reserve(nBoxesTot);
    orients.reserve(nBoxesTot);

    // Loop over all the boxes and fill the vectors used to build them
    for ( int it = 0; it < nType; it++ ) {
      for ( int iboxt = 0; iboxt < nBoxesOfType[it]; iboxt++ ) {
	dims.push_back(dimsOfType[it]);
	tols.push_back(tolsOfType[it]);
	mats.push_back(materialOfType[it]);


	// Location of the center of the box in Mu2e coords
	// Use our now-familiar trick for variable names
	std::ostringstream bCentXVarName;
	bCentXVarName << centerBaseName << "XType" << it+1 << "Box" << iboxt+1;
	std::ostringstream bCentYVarName;
	bCentYVarName << centerBaseName << "YType" << it+1 << "Box" << iboxt+1;
	std::ostringstream bCentZVarName;
	bCentZVarName << centerBaseName << "ZType" << it+1 << "Box" << iboxt+1;
	CLHEP::Hep3Vector boxCenter(c.getDouble(bCentXVarName.str())
			   *CLHEP::mm,
				    c.getDouble(bCentYVarName.str())
			   *CLHEP::mm,
				    c.getDouble(bCentZVarName.str())
			   *CLHEP::mm);
	sites.push_back(boxCenter);			   
    
	std::ostringstream bOrientVarName;
	bOrientVarName << orientBaseName << "Type" << it+1 << "Box" << iboxt+1;

	std::string orientation( c.getString(bOrientVarName.str(),"000"));
	if ( orientation.length() != 3 ) {
	  throw cet::exception("GEOM")
	    << "Orientation must be specified with a three-digit number"
	    << "\nentered as a string.  You specified: " << orientation;
	}

	orients.push_back(orientation);

      } // end loop over boxes of type...
    } // end loop over types...

    // Now make the pointer to the object itself.
    std::unique_ptr<ExtShieldUpstream> res(new ExtShieldUpstream(
								     dims,
								     tols,
								     mats,
								     sites,
								     orients)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtShieldUpstream.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
