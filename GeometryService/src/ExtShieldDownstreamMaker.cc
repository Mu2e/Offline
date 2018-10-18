//
// Original author David Norvil Brown, University of Louisville
//
// This file is the implementation of the ExtShieldDownstreamMaker class.
// It builds the ExtShieldDownstream using SimpleConfig.

#include <iostream>
#include <string>
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/ExtShieldDownstreamMaker.hh"
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
    std::string nHolesBaseName   = "ExtShieldDownstream.nHoles";
    std::string nNotchesBaseName = "ExtShieldDownstream.nNotches";
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
      std::ostringstream bOutlineUVarName;
      bOutlineUVarName << outlineBaseName << "Type" << iType << "UVerts";
      std::ostringstream bOutlineVVarName;
      bOutlineVVarName << outlineBaseName << "Type" << iType << "VVerts";
      std::vector<double> uVerts;
      std::vector<double> vVerts;
      c.getVectorDouble(bOutlineUVarName.str(),uVerts,nVert);
      c.getVectorDouble(bOutlineVVarName.str(),vVerts,nVert);
      for ( int iVert = 0; iVert < nVert; iVert++ ) {
        tempDoubleVec.push_back(uVerts[iVert]*CLHEP::mm);
        tempDoubleVec.push_back(vVerts[iVert]*CLHEP::mm);
        outlineThisType.push_back(tempDoubleVec);
        tempDoubleVec.clear();  // So it can be re-used
      }
      outlineOfType.push_back(outlineThisType);

      // Get the length of this type of box - the w component
      std::ostringstream bLengthVarName;
      bLengthVarName << lengthBaseName << "Type" << iType;
      lengthOfType.push_back(c.getDouble(bLengthVarName.str())*CLHEP::mm/2.0);

      // Get the tolerances for this type box - du,dv,dw
      std::ostringstream btolsVarName;
      btolsVarName << tolsBaseName << "Type" << iType;
      c.getVectorDouble(btolsVarName.str(),tempDoubleVec,3);
      for ( int itmp = 0; itmp < 3; itmp++ ) {
        if ( tempDoubleVec[itmp] > 10.0 || tempDoubleVec[itmp] < -20.0 ) {
          // Throw if tolerances out of tolerance.
          throw cet::exception("GEOM")
            << "Tolerance on ExternalShielding Downstream outside limits. "
            << "\nTolerances mustbe between -20 and 10 mm.";
        }
        tempDoubleVec[itmp] *= (CLHEP::mm/2.0);
      }
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
    std::vector<int>                  nHoles;
    std::vector<int>                  nNotches;
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
    nHoles                           .reserve(nBoxesTot);
    nNotches                         .reserve(nBoxesTot);
    holeID                           .reserve(nBoxesTot);
    notchID                          .reserve(nBoxesTot);

    // Helper variables for building lists of holes and notches
    int holeIdx = 0;
    int notchIdx = 0;

    // Loop over all the boxes and fill the vectors used to build them
    for ( int it = 0; it < nType; it++ ) {
      for ( int iboxt = 0; iboxt < nBoxesOfType[it]; iboxt++ ) {
        outlines.push_back(outlineOfType[it]); // This information comes
        lengths.push_back(lengthOfType[it]);   // from the type-specified
        tols.push_back(tolsOfType[it]);        // info filled above
        mats.push_back(materialOfType[it]);


        // Location of the center of the box in Mu2e coords
        // Use our now-familiar trick for variable names
        std::ostringstream bCentVarName;
        bCentVarName << centerBaseName << "Type" << it+1 << "Box" << iboxt+1;
        CLHEP::Hep3Vector boxCenter = c.getHep3Vector(bCentVarName.str());
        boxCenter *= CLHEP::mm;
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

        std::ostringstream bNHolesVarName;
        bNHolesVarName << nHolesBaseName << "Type" << it+1 << "Box" << iboxt+1;

        std::ostringstream bNNotchesVarName;
        bNNotchesVarName << nNotchesBaseName << "Type" << it+1 << "Box" << iboxt+1;

        int nHinput = c.getInt(bNHolesVarName.str(),0);
        int nNinput = c.getInt(bNNotchesVarName.str(),0);
        nHoles.push_back(nHinput);
        nNotches.push_back(nNinput);


        if ( nHinput > 0 ) {
          holeID.push_back(holeIdx);  // Keeps track of where in the list of
          // holes this block's holes begin
          for ( int iHole = 0; iHole < nHinput; iHole++ ) {
            holeIdx++; // Keep track of number of holes total
            // Location of the center of the hole in block coords
            // Use our now-familiar trick for variable names
            std::ostringstream hCentVarName;
            hCentVarName << hCentBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Hole" << iHole+1;


            CLHEP::Hep3Vector holeCenter = c.getHep3Vector(hCentVarName.str());
            holeCenter *= CLHEP::mm;


            holeLoc.push_back(holeCenter);

            // Get the hole radius
            std::ostringstream hRadVarName;
            hRadVarName << hRadBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Hole" << iHole+1;
            double tempDouble = c.getDouble(hRadVarName.str());
            holeRad.push_back(tempDouble*CLHEP::mm);
            // Get the hole length (remember G4 uses half-length)
            std::ostringstream hLenVarName;
            hLenVarName << hLenBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Hole" << iHole+1;
            tempDouble = c.getDouble(hLenVarName.str());
            holeLen.push_back(tempDouble*CLHEP::mm/2.0);
            // Get the hole orientation (uses same convention as blocks)
            std::ostringstream hOriVarName;
            hOriVarName << hOrBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Hole" << iHole+1;
            std::string orio( c.getString(hOriVarName.str(),"000"));
            if ( orio.length() != 3 ) {
              throw cet::exception("GEOM")
                << "Hole orientation must be specified with a three-digit number"
                << "\nentered as a string.You specified: " << orio;
            }
            holeOri.push_back(orio);
          } // End of loop over holes for this block

        } else {
          holeID.push_back(-1);  // Indicates this block has no holes
        } // End checking for holes for this block

        if ( nNinput > 0 ) {
          notchID.push_back(notchIdx); // Record starting point in list of
          // notches for the notches belonging to this block.
          for ( int iNotch = 0; iNotch<nNinput; iNotch++ ) {
            notchIdx++; // Keep  track of the number of notches

            // Location of the center of the notch in block coords
            // Use our now-familiar trick for variable names
            std::ostringstream nCentVarName;
            nCentVarName << nCentBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Notch" << iNotch+1;

            CLHEP::Hep3Vector notchCenter = c.getHep3Vector(nCentVarName.str());
            notchCenter *= CLHEP::mm;
            notchLoc.push_back(notchCenter);

            // And now get the dimensions of the box that defines the notch
            std::ostringstream nDimVarName;
            nDimVarName << nDimBaseName << "Type" << it+1 << "Box" << iboxt+1 << "Notch" << iNotch+1;

            std::vector<double> tempDubs;
            c.getVectorDouble(nDimVarName.str(),tempDubs,3);
            for ( int idummy = 0; idummy < 3; idummy++) tempDubs[idummy] *= CLHEP::mm/2.0;
            notchDim.push_back(tempDubs);
          } // End of loop over notches for this block
        } else {
          notchID.push_back(-1); // Indicates no notches for this block
        } // End of checking for notches for this block


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
                                                                     nHoles,
                                                                     nNotches,
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
