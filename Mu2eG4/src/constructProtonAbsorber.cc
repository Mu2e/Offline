//
// Free function to create Proton Absorber
//
// $Id: constructProtonAbsorber.cc,v 1.10 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/27 06:05:35 $
//
// Original author KLG based on Mu2eWorld constructProtonAbs
//
// Notes:
// Construct the  Proton Absorber

// C++ includes
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestCons.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4BooleanSolid.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;

namespace mu2e {

  void constructProtonAbsorber( SimpleConfig const * const _config
                                ){

    if( !_config->getBool("hasProtonAbsorber", true) ) return;

    int  const verbosityLevel           = _config->getInt("protonabsorber.verbosityLevel", 0);

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    VolumeInfo const & parent1Info  = _helper->locateVolInfo("ToyDS2Vacuum");
    VolumeInfo const & parent2Info  = _helper->locateVolInfo("ToyDS3Vacuum");

    // smaller and larger outer radii
    double pabs1rOut0   = _config->getDouble("protonabsorber.OutRadius0");
    double pabs2rOut1   = _config->getDouble("protonabsorber.OutRadius1");
    double pabsZHalfLen = _config->getDouble("protonabsorber.halfLength");
    double thick        = _config->getDouble("protonabsorber.thickness");

    // adding virtual detector before and after target
    double vdHL = 0.;
    art::ServiceHandle<GeometryService> geom;
    if( geom->hasElement<VirtualDetector>() ) {
      GeomHandle<VirtualDetector> vdg;
      if( vdg->nDet()>0 ) vdHL = vdg->getHalfLength();
    }

    // subtract virtual detector thickness from the larger outer
    // radius of the proton absorber

    pabs2rOut1 -= 2.*vdHL;

    double pabs1rIn0  = pabs1rOut0 - thick;
    double pabs2rIn1  = pabs2rOut1 - thick;

    MaterialFinder materialFinder(*_config);
    G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");


    double z0DSup   = parent1Info.centerInMu2e().z();

    // the target info should be gotten from the geometry service... not as done here

    vector<double> targetRadius;  _config->getVectorDouble("target.radii", targetRadius);

    double numoftf = (targetRadius.size()-1.0)*0.5;

    double foilwid=_config->getDouble("target.deltaZ"); 

    // we add space for the virtual detector here
    double taglen =(foilwid*numoftf) + 5.0 + 2.*vdHL; // what/why is 5.0 hardcoded here?

    double z0valt =_config->getDouble("target.z0");
    double tagoff =z0valt - z0DSup + 12000.0;

    double targetEnd = tagoff + taglen;
    double ds2HalfLen = _config->getDouble("toyDS2.halfLength");
    double pabs1len = ds2HalfLen - targetEnd;
    G4ThreeVector  pabs1Offset(0.0, 0.0, (pabs1len*0.5) + targetEnd);

    if ( verbosityLevel > 0) {
      cout << __func__ <<
        " z0DSup                            : " << z0DSup << endl;
      cout << __func__ <<
        " tagoff                            : " << tagoff << endl;
      cout << __func__ << 
	" ds2HalfLen                        : " << ds2HalfLen << endl;
      cout << __func__ << 
	" pabs1len                          : " << pabs1len << endl;
    }

    // interpolating the outer radius of the DS2 part

    double pabs1rOut1 = ((pabs2rOut1 - pabs1rOut0)*(pabs1len/(2.0*pabsZHalfLen))) + pabs1rOut0;
    double pabs1rIn1  = pabs1rOut1 - thick;

    G4Tubs const * parent1solid0 = static_cast<G4Tubs*>(parent1Info.solid);
    G4Tubs const * parent2solid0 = static_cast<G4Tubs*>(parent2Info.solid->GetConstituentSolid(0));
    double ds3HalfLen = _config->getDouble("toyDS3.halfLength");

    if ( verbosityLevel > 0) {
      cout << __func__ << 
	" ds3HalfLen                        : " << ds3HalfLen << endl;
      cout << __func__ << 
	" ds3HalfLen from ConstituentSolid(0): " <<
        parent2solid0->GetZHalfLength()<< endl; 
    }

    double pabs2len  = (2.0*pabsZHalfLen) - pabs1len;

    // double pabs2ZOffset = (pabs2len*0.5) - ds3HalfLen; 

    // we get the ds3HalfLen from the G4 physical volume itself:
    // ConstituentSolid(0) of a boolean solid, here a subtraction solid

    double pabs2ZOffset = (pabs2len*0.5) - parent2solid0->GetZHalfLength(); 

    // protonabs2 should touch protonabs1 and both of them should touch the ds2/ds3 boundary
    // it looks like the boolean solid center is in the center of ConstituentSolid(0)
    // note that the subtraction solid may be shorter, depending how the subtraction is done

    if ( verbosityLevel > 0) {
      double theZ  = parent1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
      double theHL = parent1solid0->GetZHalfLength();
      cout << __func__ << " " << parent1Info.name << " Z offset in Mu2e    : " <<
	theZ << endl;
      cout << __func__ << " " << parent1Info.name << " Z extent in Mu2e    : " <<
        theZ - theHL << ", " << theZ + theHL << endl;
    }

    if ( verbosityLevel > 0) {
      double theZ  = parent2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
      double theHL = parent2solid0->GetZHalfLength();
      cout << __func__ << " " << parent2Info.name << " Z offset in Mu2e    : " <<
	theZ << endl;
      cout << __func__ << " " << parent2Info.name << " Z extent in Mu2e    : " <<
        theZ - theHL << ", " << theZ + theHL << endl;
    }

    G4ThreeVector pabs2Offset(0.0, 0.0, pabs2ZOffset);

    // proton absorber in DS2
    double pabs1Param[7] = { pabs1rIn0, pabs1rOut0, pabs1rIn1, pabs1rOut1, pabs1len*0.5,
			     0.0, 360.0*CLHEP::degree };
    bool pabsIsVisible = _config->getBool("protonabsorber.visible",true);
    bool pabsIsSolid   = _config->getBool("protonabsorber.solid",true);

    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    mf::LogInfo log("GEOM");
    log << "Constructing Proton Absorber -- \n";
    log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
    log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
    log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
    log << "halflength: "<< pabs1len*0.5 <<"\n";
    log << "Proton Abs Offset in DS3:  " << pabs2Offset <<"\n";
    log << "rIn,  rOut (-z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
    log << "rIn,  rOut (+z): "<< pabs2rIn1 <<"  "<< pabs2rOut1<<"  ";
    log << "halflength: "<< pabs2len*0.5 <<"\n";

    VolumeInfo protonabs1Info = nestCons( "protonabs1",
                                          pabs1Param,
                                          pabsMaterial,
                                          0,
                                          pabs1Offset,
                                          parent1Info,
                                          0,
                                          pabsIsVisible,
                                          G4Color::White(),
                                          pabsIsSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    if ( verbosityLevel > 0) {
      double pzhl   = static_cast<G4Cons*>(protonabs1Info.solid)->GetZHalfLength();
      double pabs1Z = protonabs1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
      cout << __func__ << " protonabs1 Z offset in Mu2e    : " <<
        pabs1Z << endl;
      cout << __func__ << " protonabs1 Z extent in Mu2e    : " <<
        pabs1Z - pzhl  << ", " <<  pabs1Z + pzhl  << endl;
    }


    // proton absorber in DS3
    double pabs2Param[7] = { pabs1rIn1, pabs1rOut1, pabs2rIn1, pabs2rOut1, pabs2len*0.5, 
                             0.0, 360.0*CLHEP::degree };

    VolumeInfo protonabs2Info = nestCons( "protonabs2",
                                          pabs2Param,
                                          pabsMaterial,
                                          0,
                                          pabs2Offset,
                                          parent2Info,
                                          0,
                                          pabsIsVisible,
                                          G4Color::Yellow(),
                                          pabsIsSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );
    if ( verbosityLevel > 0) {
      double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();
      double pabs2Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
      cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
        pabs2Z << endl;
      cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
        pabs2Z - pzhl  << ", " <<  pabs2Z + pzhl  << endl;

      // we also check how the offsets are handled

      cout << __func__ << " " << protonabs2Info.name << " local input offset in G4                  : " <<
        pabs2Offset << endl;
      cout << __func__ << " " << protonabs2Info.name << " local GetTranslation()       offset in G4 : " <<
        protonabs2Info.physical->GetTranslation() << endl; // const &

    }

  }

}
