//
// Free function to create Proton Absorber
//
// $Id: constructProtonAbsorber.cc,v 1.4 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
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
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestCons.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"

using namespace std;

namespace mu2e {

  void constructProtonAbsorber( SimpleConfig const * const _config
                                ){
    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    VolumeInfo const & parent1  = _helper->locateVolInfo("ToyDS2Vacuum");
    VolumeInfo const & parent2  = _helper->locateVolInfo("ToyDS3Vacuum");
    VolumeInfo const & hallInfo = _helper->locateVolInfo("HallAir");

    double pabs1rOut0 = _config->getDouble("protonabsorber.OutRadius0");

    double pabs2rOut1 = _config->getDouble("protonabsorber.OutRadius1");
    double zLen       = _config->getDouble("protonabsorber.halfLength");
    double thick      = _config->getDouble("protonabsorber.thickness");

    double pabs1rIn0  = pabs1rOut0 - thick;
    double pabs2rIn1  = pabs2rOut1 - thick;

    MaterialFinder materialFinder(*_config);
    G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");

    // Add virtual detector before and after target
    double vdHL = 0.;
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<VirtualDetector>() ) return;
    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()>0 ) vdHL = vdg->getHalfLength();

    double z0DSup   = parent1.centerInWorld.z()+hallInfo.centerInMu2e().z();
    vector<double> targetRadius;  _config->getVectorDouble("target.radii", targetRadius);
    double numoftf = (targetRadius.size()-1)/2.0;
    double foilwid=_config->getDouble("target.deltaZ"); double taglen =(foilwid*numoftf) + 5.0 + vdHL;
    double z0valt =_config->getDouble("target.z0");     double tagoff =z0valt - z0DSup + 12000.0;
    double targetEnd = tagoff + taglen;
    double ds2halflen = _config->getDouble("toyDS2.halfLength");
    double pabs1len = ds2halflen - targetEnd;
    G4ThreeVector  pabs1Offset(0.0, 0.0, (pabs1len/2.0) + targetEnd);

    double pabs1rOut1 = ((pabs2rOut1 - pabs1rOut0)*(pabs1len/(2.0*zLen))) + pabs1rOut0;
    double pabs1rIn1  = pabs1rOut1 - thick;
    double ds3halflen = _config->getDouble("toyDS3.halfLength");
    double pabs2len  = (2.0*zLen) - pabs1len;
    double pabs2zoff = (pabs2len/2.0) - ds3halflen;
    G4ThreeVector  pabs2Offset(0.0, 0.0, pabs2zoff);

    // proton absorber in DS2
    double pabs1Param[7] = { pabs1rIn0, pabs1rOut0, pabs1rIn1, pabs1rOut1, pabs1len/2.0, 0.0, 360.0*CLHEP::degree };
    bool pabsVisible = _config->get<bool>("protonabsorber.visible",true);
    bool pabsSolid   = _config->get<bool>("protonabsorber.solid",true);

    bool forceAuxEdgeVisible = _config->get<bool>("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->get<bool>("g4.doSurfaceCheck",false);
    bool const placePV       = true;


    if( _config->get<bool>("hasProtonAbsorber", true) ){

      mf::LogInfo log("GEOM");
      log << "Constructing Proton Absorber -- \n";
      log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
      log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
      log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
      log << "halflength: "<< pabs1len/2.0 <<"\n";
      log << "Proton Abs Offset in DS3:  " << pabs2Offset <<"\n";
      log << "rIn,  rOut (-z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
      log << "rIn,  rOut (+z): "<< pabs2rIn1 <<"  "<< pabs2rOut1<<"  ";
      log << "halflength: "<< pabs2len/2.0 <<"\n";

      VolumeInfo protonabs1Info = nestCons( "protonabs1",
                                            pabs1Param,
                                            pabsMaterial,
                                            0,
                                            pabs1Offset,
                                            parent1,
                                            0,
                                            pabsVisible,
                                            G4Color::White(),
                                            pabsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

      // proton absorber in DS3
      double pabs2Param[7] = { pabs1rIn1, pabs1rOut1, pabs2rIn1, pabs2rOut1, pabs2len/2.0, 0.0, 360.0*CLHEP::degree };

      VolumeInfo protonabs2Info = nestCons( "protonabs2",
                                            pabs2Param,
                                            pabsMaterial,
                                            0,
                                            pabs2Offset,
                                            parent2,
                                            0,
                                            pabsVisible,
                                            G4Color::Magenta(),
                                            pabsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    }
  }

}
