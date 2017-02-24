//
// Original author David Norvil Brown, University of Louisville
//
// This file is the implementation of the ElectronicRackMaker class.
// It builds the ElectronicRack using SimpleConfig.

#include <string>
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/ElectronicRackMaker.hh"
#include "ServicesGeom/inc/ElectronicRack.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ElectronicRack>  ElectronicRackMaker::make(const
                                                              SimpleConfig& c)
  {

    // We use just a few basic types, so simplify config file by
    // only specifying dimensions, tolerances, and material once per type.
    // Locations and orientations (?) will have to be entered for each rack.

    const int nType = c.getInt("ElectronicRack.numberOfRackTypes",0);

    std::vector<int>                       nRacksOfType;
    std::vector<std::vector<double> >      dimsOfType;
    std::vector<std::string>               materialOfType;
    nRacksOfType                           .reserve(nType);
    dimsOfType                             .reserve(nType);
    materialOfType                         .reserve(nType);

    std::string typeBaseName     = "ElectronicRack.nRackType";
    std::string dimsBaseName     = "ElectronicRack.dims";
    std::string materialBaseName = "ElectronicRack.materialType";
    std::string centerBaseName   = "ElectronicRack.center";
    std::string orientBaseName   = "ElectronicRack.orientation";

    // Some temporary holders
    std::vector<double> tempDoubleVec;

    // Loop over the various rack types and fill the vectors with
    // information needed for the types.
    for ( int iType = 1; iType <= nType; iType++ ) {
      // A stream trick for creating variable variable names!
      std::ostringstream bTypeNumberVarName;
      bTypeNumberVarName << typeBaseName << iType;
      nRacksOfType.push_back(c.getInt(bTypeNumberVarName.str()));

      // Get the dimensions for this type rack - u,v,w are the self-coordinates
      std::ostringstream bdimsVarName;
      bdimsVarName << dimsBaseName << "Type" << iType;
      // Divide dimensions by 2.0 because G4 racks are created in terms of
      // half-lengths.
      c.getVectorDouble(bdimsVarName.str(),tempDoubleVec,3);
      for ( int itmp = 0; itmp < 3; itmp++ ) tempDoubleVec[itmp]*=(CLHEP::mm/2.0);
      dimsOfType.push_back(tempDoubleVec);
      tempDoubleVec.clear();  // So it can be re-used

      // Get the material name for this type rack
      std::ostringstream bmatVarName;
      bmatVarName << materialBaseName << iType;
      materialOfType.push_back( c.getString(bmatVarName.str()));
    }


    // Get total number of racks from info collected above
    int nRacksTot = 0;
    for ( unsigned int ib = 0; ib < nRacksOfType.size(); ib++ ) {
      nRacksTot += nRacksOfType[ib];
    }

    // Set up the actual vectors used to build the upstream shielding
    std::vector<std::vector<double> > dims;
    std::vector<std::string>          mats;
    std::vector<CLHEP::Hep3Vector>    sites;
    std::vector<std::string>          orients;
    dims.reserve(nRacksTot);
    mats.reserve(nRacksTot);
    sites.reserve(nRacksTot);
    orients.reserve(nRacksTot);

    // Loop over all the racks and fill the vectors used to build them
    for ( int it = 0; it < nType; it++ ) {
      for ( int irackt = 0; irackt < nRacksOfType[it]; irackt++ ) {
        dims.push_back(dimsOfType[it]);
        mats.push_back(materialOfType[it]);


        // Location of the center of the rack in Mu2e coords
        // Use our now-familiar trick for variable names
        std::ostringstream bCentVarName;
        bCentVarName << centerBaseName << "Type" << it+1 << "Rack" << irackt+1;
        CLHEP::Hep3Vector rackCenter = c.getHep3Vector(bCentVarName.str());
        rackCenter *= CLHEP::mm;
        sites.push_back(rackCenter);

        std::ostringstream bOrientVarName;
        bOrientVarName << orientBaseName << "Type" << it+1 << "Rack" << irackt+1;

        std::string orientation( c.getString(bOrientVarName.str(),"000"));
        if ( orientation.length() != 3 ) {
          throw cet::exception("GEOM")
            << "Orientation must be specified with a three-digit number"
            << "\nentered as a string.  You specified: " << orientation;
        }

        orients.push_back(orientation);

      } // end loop over racks of type...
    } // end loop over types...

      // Now make the pointer to the object itself.
    std::unique_ptr<ElectronicRack> res(new ElectronicRack(
							   dims,
							   mats,
							   sites,
							   orients)
					);
    
    //----------------------------------------------------------------
    if(c.getInt("ElectronicRack.verbosityLevel",0) > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }
} // namespace mu2e
