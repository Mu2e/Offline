//
// Original author David Norvil Brown, University of Louisville
//
// This file is the implementation of the ExtShieldDownstreamMaker class.
// It builds the ExtShieldDownstream using SimpleConfig.

#include <iostream>
#include <string>
#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalShieldingGeom/inc/ExtShieldDownstreamMaker.hh"
#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtShieldDownstream>  ExtShieldDownstreamMaker::make(const 
							      SimpleConfig& c)
  {

    // We use just a few basic types, so simplify config file by
    // only specifying dimensions, tolerances, and material once per type.
    // Locations and orientations (?) will have to be entered for each box.

    const int nType = c.getInt("ExtShieldDownstream.numberOfBoxTypes");

    std::vector<int>                       nBoxesOfType;
    std::vector<std::vector<std::vector<double> > >      outlineOfType;
    std::vector<double>                    lengthOfType;
    std::vector<std::vector<double> >      tolsOfType;
    std::vector<std::string>               materialOfType;
    nBoxesOfType                           .reserve(nType);
    outlineOfType                          .reserve(nType);
    lengthOfType                           .reserve(nType);
    tolsOfType                             .reserve(nType);
    materialOfType                         .reserve(nType);

    std::string typeBaseName     = "ExtShieldDownstream.nBoxType";
    std::string vertBaseName     = "ExtShieldDownstream.nVertType";
    std::string outlineBaseName  = "ExtShieldDownstream.outline";
    std::string lengthBaseName   = "ExtShieldDownstream.length";
    std::string tolsBaseName     = "ExtShieldDownstream.tols";
    std::string materialBaseName = "ExtShieldDownstream.materialType";
    std::string centerBaseName   = "ExtShieldDownstream.center";
    std::string orientBaseName   = "ExtShieldDownstream.orientation";
    std::string hasHoleBaseName  = "ExtShieldDownstream.hasHole";
    std::string hasNotchBaseName = "ExtShieldDownstream.hasNotch";
    std::string hCentBaseName    = "ExtShieldDownstream.holeCenter";
    std::string hRadBaseName     = "ExtShieldDownstream.holeRadius";
    std::string hLenBaseName     = "ExtShieldDownstream.holeLength";
    std::string hOrBaseName      = "ExtShieldDownstream.holeOrientation";
    std::string nCentBaseName    = "ExtShieldDownstream.notchCenter";
    std::string nDimBaseName     = "ExtShieldDownstream.notchDim";

    // Some temporary holders
    std::vector<double> tempDoubleVec;

    // Loop over the various box types and fill the vectors with 
    // information needed for the types.
    for ( int iType = 1; iType <= nType; iType++ ) {
      // A stream trick for creating variable variable names!
      std::ostringstream bTypeNumberVarName;
      bTypeNumberVarName << typeBaseName << iType;
      nBoxesOfType.push_back(c.getInt(bTypeNumberVarName.str()));

      std::ostringstream bnVertsVarName;
      bnVertsVarName << vertBaseName << iType;
      int nVert = c.getInt(bnVertsVarName.str());

      // Get the vertices for this type box - u,v,w are the self-coordinates
      // Specify shape in u-v plane and extrude in w
      // Loop over all the vertices
      std::vector<std::vector<double> > outlineThisType;
      outlineThisType.reserve(nVert);
      for ( int iVert = 1; iVert <= nVert; iVert++ ) {
	std::ostringstream bOutlineUVarName;
	bOutlineUVarName << outlineBaseName << "Type" << iType << "UVert" << iVert;
	std::ostringstream bOutlineVVarName;
	bOutlineVVarName << outlineBaseName << "Type" << iType << "VVert" << iVert;
	tempDoubleVec.push_back(c.getDouble(bOutlineUVarName.str())*CLHEP::mm);
	tempDoubleVec.push_back(c.getDouble(bOutlineVVarName.str())*CLHEP::mm);
	outlineThisType.push_back(tempDoubleVec);
	tempDoubleVec.clear();  // So it can be re-used
      }
      outlineOfType.push_back(outlineThisType);

      // Get the length of this type of box - the w component
      std::ostringstream bLengthVarName;
      bLengthVarName << lengthBaseName << "Type" << iType;
      lengthOfType.push_back(c.getDouble(bLengthVarName.str())*CLHEP::mm/2.0);

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
      if ( du >10.0 ||dv > 10.0 || dw> 10.0 || du < -20.0 ||dv < -20.0 || dw < -20.0 ) {
        // Throw if tolerances out of tolerance.
        throw cet::exception("GEOM")
          << "Tolerance on ExternalShielding Downstream outside limits. "
          << "\nTolerances mustbe between -20 and 10 mm.";
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

    // ***************************************************
    // Above we set up the type-by-type information
    // Below we set up the specific box-by-box information
    // ***************************************************

    // Set up the actual vectors used to build the downstream shielding
    std::vector<std::vector<std::vector<double> > > outlines;
    std::vector<double>               lengths;
    std::vector<std::vector<double> > tols;
    std::vector<std::string>          mats;
    std::vector<CLHEP::Hep3Vector>    sites;
    std::vector<std::string>          orients;
    std::vector<bool>                 holeYes;
    std::vector<bool>                 notchYa;
    std::vector<int>                  holeID;
    std::vector<int>                  notchID;
    std::vector<CLHEP::Hep3Vector>    holeLoc;
    std::vector<double>               holeRad;
    std::vector<double>               holeLen;
    std::vector<std::string>          holeOri;
    std::vector<CLHEP::Hep3Vector>    notchLoc;
    std::vector<std::vector<double> > notchDim;
    outlines                         .reserve(nBoxesTot);
    lengths                          .reserve(nBoxesTot);
    tols                             .reserve(nBoxesTot);
    mats                             .reserve(nBoxesTot);
    sites                            .reserve(nBoxesTot);
    orients                          .reserve(nBoxesTot);
    holeYes                          .reserve(nBoxesTot);
    notchYa                          .reserve(nBoxesTot);
    holeID                           .reserve(nBoxesTot);
    notchID                          .reserve(nBoxesTot);

    // Helper variables for building lists of holes and notches
    int holeIdx = 0;
    int notchIdx = 0;

    // Loop over all the boxes and fill the vectors used to build them
    for ( int it = 0; it < nType; it++ ) {
      for ( int iboxt = 0; iboxt < nBoxesOfType[it]; iboxt++ ) {
	outlines.push_back(outlineOfType[it]);
	lengths.push_back(lengthOfType[it]);
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
            << "\nentered as a string.You specified: " << orientation;
	}
	orients.push_back(orientation);

	// ********************************************
	// Now handle special stuff - holes and notches
	// ********************************************

	std::ostringstream bHasHoleVarName;
	bHasHoleVarName << hasHoleBaseName << "Type" << it+1 << "Box" << iboxt+1;

	std::ostringstream bHasNotchVarName;
	bHasNotchVarName << hasNotchBaseName << "Type" << it+1 << "Box" << iboxt+1;
	//	std::cout << "DNB** About to get bools for holes/notches for type " << it+1 << ", Box " << iboxt+1  << std::endl;
	bool hasHole = c.getBool(bHasHoleVarName.str(),false);
	bool hasNotch = c.getBool(bHasNotchVarName.str(),false);
	holeYes.push_back(hasHole);
	notchYa.push_back(hasNotch);

	if ( hasHole ) { 
	  holeID.push_back(holeIdx);
	  holeIdx++;
	  // Location of the center of the hole in block coords
	  // Use our now-familiar trick for variable names
	  std::ostringstream hCentUVarName;
	  hCentUVarName << hCentBaseName << "UType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream hCentVVarName;
	  hCentVVarName << hCentBaseName << "VType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream hCentWVarName;
	  hCentWVarName << hCentBaseName << "WType" << it+1 << "Box" << iboxt+1;
	  //	  std::cout << "DNB** About to read hole center for type " << it+1 << ", Box " << iboxt+1 << std::endl;

	  CLHEP::Hep3Vector holeCenter(c.getDouble(hCentUVarName.str())
				      *CLHEP::mm,
				      c.getDouble(hCentVVarName.str())
				      *CLHEP::mm,
				      c.getDouble(hCentWVarName.str())
				      *CLHEP::mm);
	  holeLoc.push_back(holeCenter);			   

	  // Get the hole radius 
	  std::ostringstream hRadVarName;
	  hRadVarName << hRadBaseName << "Type" << it+1 << "Box" << iboxt+1;
	  double tempDouble = c.getDouble(hRadVarName.str());
	  holeRad.push_back(tempDouble*CLHEP::mm);
	  // Get the hole length (remember G4 uses half-length)
	  std::ostringstream hLenVarName;
	  hLenVarName << hLenBaseName << "Type" << it+1 << "Box" << iboxt+1;
	  tempDouble = c.getDouble(hLenVarName.str());
	  holeLen.push_back(tempDouble*CLHEP::mm/2.0);
	  // Get the hole orientation (uses same convention as blocks)
	  std::ostringstream hOriVarName;
	  hOriVarName << hOrBaseName << "Type" << it+1 << "Box" << iboxt+1;
	  std::string orio( c.getString(hOriVarName.str(),"000"));
	  if ( orio.length() != 3 ) {
	    throw cet::exception("GEOM")
	      << "Hole orientation must be specified with a three-digit number"
	      << "\nentered as a string.You specified: " << orio;
	  }
	  holeOri.push_back(orio);

	} else {
	  holeID.push_back(-1);
	}

	if ( hasNotch ) { 
	  notchID.push_back(notchIdx);
	  notchIdx++;
	  // Location of the center of the notch in block coords
	  // Use our now-familiar trick for variable names
	  std::ostringstream nCentUVarName;
	  nCentUVarName << nCentBaseName << "UType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream nCentVVarName;
	  nCentVVarName << nCentBaseName << "VType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream nCentWVarName;
	  nCentWVarName << nCentBaseName << "WType" << it+1 << "Box" << iboxt+1;

	  CLHEP::Hep3Vector notchCenter(c.getDouble(nCentUVarName.str())
					*CLHEP::mm,
					c.getDouble(nCentVVarName.str())
					*CLHEP::mm,
					c.getDouble(nCentWVarName.str())
					*CLHEP::mm);
	  notchLoc.push_back(notchCenter);			   

	  // And now get the dimensions of the box
	  std::ostringstream nDimUVarName;
	  nDimUVarName << nDimBaseName << "UType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream nDimVVarName;
	  nDimVVarName << nDimBaseName << "VType" << it+1 << "Box" << iboxt+1;
	  std::ostringstream nDimWVarName;
	  nDimWVarName << nDimBaseName << "WType" << it+1 << "Box" << iboxt+1;

	  std::vector<double> tempDubs;
	  tempDubs.push_back(c.getDouble(nDimUVarName.str())*CLHEP::mm/2.0);
	  tempDubs.push_back(c.getDouble(nDimVVarName.str())*CLHEP::mm/2.0);
	  tempDubs.push_back(c.getDouble(nDimWVarName.str())*CLHEP::mm/2.0);

	  notchDim.push_back(tempDubs);			   

	} else {
	  notchID.push_back(-1);
	}


      } // end loop over boxes of type...
    } // end loop over types...

    // Now make the pointer to the object itself.
    std::unique_ptr<ExtShieldDownstream> res(new ExtShieldDownstream(
								     outlines,
								     lengths,
								     tols,
								     mats,
								     sites,
								     orients,
								     holeYes,
								     notchYa,
								     holeID,
								     holeLoc,
								     holeRad,
								     holeLen,
								     holeOri,
								     notchID,
								     notchLoc,
								     notchDim)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtShieldDownstream.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
