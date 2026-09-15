//
// Original author David Norvil Brown, University of Louisville
//
// This file is the implementation of the PipeMaker class.
// It builds the Pipe using SimpleConfig.

#include <iostream>
#include <string>
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/GeometryService/inc/PipeMaker.hh"
#include "Offline/ServicesGeom/inc/Pipe.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<Pipe>  PipeMaker::make(const
                                             SimpleConfig& c)
  {

    // We use a few types of pipe, with several copies of each allowed.
    // A pipe contains components, each of which is also a pipe.  The
    // "components" are embedded within the first pipe, so you get a pipe
    // full of pipes.  Component positions are specified with respect to
    // the first pipe and then the first pipe is placed in the world.

    const int version = c.getInt("Pipe.version",1);

    const int nType = c.getInt("Pipe.numberOfPipeTypes");

    // Qualities of types
    std::vector<int>                       nComponentsOfType;
    std::vector<int>                       nPipesOfType;
    std::vector<Pipe::Flavor>              flavorOfType;
    std::vector<std::string>               fillOfType;
    nComponentsOfType                      .reserve(nType);
    nPipesOfType                           .reserve(nType);
    flavorOfType                           .reserve(nType);
    fillOfType                             .reserve(nType);

    // Qualities of Pipes, by type
    std::vector<std::vector<CLHEP::Hep3Vector> >         sites;
    std::vector<std::vector<std::string> >               orients;
    std::vector<std::vector<double> >                    lengths;
    sites                                  .reserve(nType);
    orients                                .reserve(nType);
    lengths                                .reserve(nType);

    // Collections by type of qualities of components.  Filled below
    std::vector<std::vector<double> >      rInOfCompByType;
    std::vector<std::vector<double> >      rOutOfCompByType;
    std::vector<std::vector<std::string> > materialOfCompByType;
    std::vector<std::vector<double> >      uOffsetsOfCompByType;
    std::vector<std::vector<double> >      vOffsetsOfCompByType;
    rInOfCompByType                        .reserve(nType);
    rOutOfCompByType                       .reserve(nType);
    materialOfCompByType                   .reserve(nType);
    uOffsetsOfCompByType                   .reserve(nType);
    vOffsetsOfCompByType                   .reserve(nType);

    // Base names of type qualities
    std::string nCompBaseName     = "Pipe.nComponentsInType";
    std::string nPipeBaseName     = "Pipe.nPipeType";
    std::string lengthBaseName    = "Pipe.lengthType";
    std::string flavBaseName      = "Pipe.flavorType";
    std::string fillBaseName      = "Pipe.fillMaterialType";

    // Base names for qualities of pipes
    std::string centerBaseName   = "Pipe.centerType";
    std::string orientBaseName   = "Pipe.orientationType";

    // Base names for qualities of components
    std::string rInBaseName      = "Pipe.rInType";
    std::string rOutBaseName     = "Pipe.rOutType";
    std::string materialBaseName = "Pipe.materialType";
    std::string uOffsetBaseName  = "Pipe.uOffsetType";
    std::string vOffsetBaseName  = "Pipe.vOffsetType";


    // *********************
    // Loop over the various pipe _types_ and fill the vectors with
    // information needed for the _types_.
    for ( int iType = 1; iType <= nType; iType++ ) {
      // *** For this Type ***
      // *** Get the number of these pipes ***
      // A stream trick for creating variable variable names!
      std::ostringstream bTypeNumberVarName;
      bTypeNumberVarName << nPipeBaseName << iType;
      nPipesOfType.push_back(c.getInt(bTypeNumberVarName.str()));

      // *** For this Type ***
      // *** Get the number of components in these pipes ***
      std::ostringstream bCompNumberVarName;
      bCompNumberVarName << nCompBaseName << iType;
      int nComp = c.getInt(bCompNumberVarName.str(),1);
      nComponentsOfType.push_back(nComp);

      // Get the flavor - straight or bend for this type
      std::ostringstream pFlavVarName;
      pFlavVarName << flavBaseName << iType;
      std::string theFlav = c.getString(pFlavVarName.str(),"straight");
      if ( theFlav == "straight" ) {
        flavorOfType.push_back(Pipe::Flavor::straight);
      } else if ( theFlav == "bend" ) {
        flavorOfType.push_back(Pipe::Flavor::bend);
      } else {
        throw cet::exception("GEOM")
          << pFlavVarName.str() << " must be \"straight\" or \"bend\"."
          << "\nYou specified: " << theFlav;
      }

      // Get the fill material for the type of pipe
      std::ostringstream pFillVarName;
      pFillVarName << fillBaseName << iType;
      std::string theFill = c.getString(pFillVarName.str());
      fillOfType.push_back(theFill);

    } // End of collecting Type information

    // ****************
    // Now collect individual pipe information
    // ****************
    std::vector<CLHEP::Hep3Vector> tmpVecHep3V;
    std::vector<std::string> tmpVecOri;
    std::vector<double> tmpVecLen;
    for ( int it = 0; it < nType; it++ ) {
      tmpVecHep3V.clear();
      tmpVecOri.clear();
      tmpVecLen.clear();
      for ( int ip = 0; ip < nPipesOfType[it]; ip++ ) {
        // Location of the center of the pipe in Mu2e coords
        // Use our now-familiar trick for variable names
        std::ostringstream bCentVarName;
        bCentVarName << centerBaseName << it+1 << "Pipe" << ip+1;
        CLHEP::Hep3Vector pipeCenter = c.getHep3Vector(bCentVarName.str());
        pipeCenter *= CLHEP::mm;
        tmpVecHep3V.push_back(pipeCenter);

        std::ostringstream bOrientVarName;
        bOrientVarName << orientBaseName << it+1 << "Pipe" << ip+1;

        std::string orientation( c.getString(bOrientVarName.str(),"000"));
        if ( orientation.length() != 3 ) {
          throw cet::exception("GEOM")
            << "Orientation must be specified with a three-digit number"
            << "\nentered as a string.You specified: " << orientation;
        }
        tmpVecOri.push_back(orientation);

        // Full length of the pipe - will be cut in half when making Tubs,
        // not Tori.  Version 1 gives one length per type.  Version 2 allows
        // a length per pipe, falling back to the length of the type.
        std::ostringstream pipeLengthVarName;
        pipeLengthVarName << lengthBaseName << it+1 << "Pipe" << ip+1;
        std::ostringstream typeLengthVarName;
        typeLengthVarName << lengthBaseName << it+1;

        double pipeLength = 0.0;
        if ( version > 1 && c.hasName(pipeLengthVarName.str()) ) {
          pipeLength = c.getDouble(pipeLengthVarName.str());
        } else {
          pipeLength = c.getDouble(typeLengthVarName.str());
        }
        tmpVecLen.push_back(pipeLength*CLHEP::mm);

      }
      orients.push_back(tmpVecOri);
      sites.push_back(tmpVecHep3V);
      lengths.push_back(tmpVecLen);
    }

    // ****************
    // Now get the information for each component.


    // ***************************************************
    // Above we set up the type-by-type information
    // Below we set up the specific component-by-component information
    // ***************************************************

    // Set up the actual vectors used to build the Pipe systems
    // Qualities of individual components


    std::vector<double> tmpVecRin;
    std::vector<double> tmpVecRout;
    std::vector<std::string> tmpVecMat;
    std::vector<double> tmpVecuOff;
    std::vector<double> tmpVecvOff;

    // Loop over all the tubes and fill the vectors used to build them
    for ( int it = 0; it < nType; it++ ) {
      tmpVecRin.clear();
      tmpVecRout.clear();
      tmpVecMat.clear();
      tmpVecuOff.clear();
      tmpVecvOff.clear();
      for ( int ipipet = 0; ipipet < nComponentsOfType[it]; ipipet++ ) {
        // holder variable
        double tmpDub = 0.0;

        // Get the inner radii of components
        std::ostringstream pRInVarName;
        pRInVarName << rInBaseName << it+1 << "Comp" << ipipet+1;
        tmpDub = c.getDouble(pRInVarName.str())*CLHEP::mm;
        tmpVecRin.push_back(tmpDub);

        // Get the outer radii of components
        std::ostringstream pROutVarName;
        pROutVarName << rOutBaseName << it+1 << "Comp" << ipipet+1;
        tmpDub = c.getDouble(pROutVarName.str())*CLHEP::mm;
        tmpVecRout.push_back(tmpDub);

        // Get the "u" offsets of centers of components
        std::ostringstream puOffVarName;
        puOffVarName << uOffsetBaseName << it+1 << "Comp" << ipipet+1;
        tmpDub = c.getDouble(puOffVarName.str())*CLHEP::mm;
        tmpVecuOff.push_back(tmpDub);

        // Get the "v" offsets of centers of components
        std::ostringstream pvOffVarName;
        pvOffVarName << vOffsetBaseName << it+1 << "Comp" << ipipet+1;
        tmpDub = c.getDouble(pvOffVarName.str())*CLHEP::mm;
        tmpVecvOff.push_back(tmpDub);

        // Get the material of each component
        std::ostringstream pMatVarName;
        pMatVarName << materialBaseName << it+1 << "Comp" << ipipet+1;
        tmpVecMat.push_back(c.getString(pMatVarName.str()));

      } // end loop over components of type...

      rInOfCompByType        .push_back(tmpVecRin);
      rOutOfCompByType       .push_back(tmpVecRout);
      materialOfCompByType   .push_back(tmpVecMat);
      uOffsetsOfCompByType   .push_back(tmpVecuOff);
      vOffsetsOfCompByType   .push_back(tmpVecvOff);

    } // end loop over types...

    // Now make the pointer to the object itself.
    std::unique_ptr<Pipe> res(new Pipe( nComponentsOfType,
                                        nPipesOfType,
                                        flavorOfType,
                                        fillOfType,
                                        sites,
                                        orients,
                                        lengths,
                                        rInOfCompByType,
                                        rOutOfCompByType,
                                        materialOfCompByType,
                                        uOffsetsOfCompByType,
                                        vOffsetsOfCompByType)
                              );

    //----------------------------------------------------------------
    if(c.getInt("Pipe.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
