// $Id: constructExternalNeutronShielding.cc,v 1.4 2014/06/05 21:08:24 genser Exp $
// $Author: genser $
// $Date: 2014/06/05 21:08:24 $
// David Norvil Brown, August 2013
//
//
// Modified by K.L.Genser to make the windows using G4SubtractionSolid

#include "Mu2eG4/inc/constructExternalNeutronShielding.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Include each shield here...
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1a.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1b.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream2.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamTop.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamBottom.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRight.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRightb.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexLeft.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRoof.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLAbove.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLCeiling.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCryoBoxes.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxes.hh"
// etc...
#include "GeometryService/inc/GeomHandle.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
//#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"

#include <vector>
#include <sstream>
#include <algorithm>


namespace mu2e {

  //================================================================

  void constructExternalNeutronShielding(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ExtNeutShieldUpstream1a> ensu1a;
    // Do similar thing for each of the sub-shields...
    GeomHandle<ExtNeutShieldUpstream1b> ensu1b;
    GeomHandle<ExtNeutShieldUpstream2> ensu2;
    GeomHandle<ExtNeutShieldUpstreamTop> ensutop;
    GeomHandle<ExtNeutShieldUpstreamBottom> ensubot;
    GeomHandle<ExtNeutShieldCavexRight> enscer;
    GeomHandle<ExtNeutShieldCavexRightb> enscerb;
    GeomHandle<ExtNeutShieldCavexLeft> enscel;
    GeomHandle<ExtNeutShieldCavexRoof> enscero;
    GeomHandle<ExtNeutShieldLAbove> ensla;
    GeomHandle<ExtNeutShieldLCeiling> enslc;
    GeomHandle<ExtNeutShieldCryoBoxes> enscb;
    GeomHandle<ExtNeutShieldCendBoxes> enscendb;


    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    //----------------------------------------------------------------
    // ExtNeutShieldUpstream1a <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensu1aOutline; 
    std::copy(ensu1a->externalShieldOutline().begin(),
	      ensu1a->externalShieldOutline().end(),
	      std::back_inserter(ensu1aOutline));

    VolumeInfo ensu1aShield("ExtNeutShieldUpstream1a",
                            ensu1a->centerOfShield()-parent.centerInMu2e(),
                            parent.centerInWorld);

    ensu1aShield.solid = new G4ExtrudedSolid(ensu1aShield.name,
                                             ensu1aOutline,
                                             ensu1a->getLength()/2.,
                                             G4TwoVector(0,0), 1.,
                                             G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensu1aRotat(CLHEP::HepRotation::IDENTITY);
    ensu1aRotat.rotateX(ensu1a->getRotation());

    finishNesting(ensu1aShield,
 		  findMaterialOrThrow(ensu1a->materialName()),
 		  &ensu1aRotat,
 		  ensu1aShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldUpstream1a.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldUpstream1a.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);


    //----------------------------------------------------------------
    // ExtNeutShieldUpstream1b <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensu1bOutline; 
    std::copy(ensu1b->externalShieldOutline().begin(),
	      ensu1b->externalShieldOutline().end(),
	      std::back_inserter(ensu1bOutline));

    VolumeInfo ensu1bShield("ExtNeutShieldUpstream1b",
                            ensu1b->centerOfShield()-parent.centerInMu2e(),
                            parent.centerInWorld);

    ensu1bShield.solid = new G4ExtrudedSolid(ensu1bShield.name,
                                             ensu1bOutline,
                                             ensu1b->getLength()/2.,
                                             G4TwoVector(0,0), 1.,
                                             G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensu1bRotat(CLHEP::HepRotation::IDENTITY);
    ensu1bRotat.rotateX(ensu1b->getRotation());

    finishNesting(ensu1bShield,
 		  findMaterialOrThrow(ensu1b->materialName()),
 		  &ensu1bRotat,
 		  ensu1bShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldUpstream1b.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldUpstream1b.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

  
    //----------------------------------------------------------------
    // ExtNeutShieldUpstream2 <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensu2Outline; 
    std::copy(ensu2->externalShieldOutline().begin(),
	      ensu2->externalShieldOutline().end(),
	      std::back_inserter(ensu2Outline));

    VolumeInfo ensu2Shield("ExtNeutShieldUpstream2",
                           ensu2->centerOfShield()-parent.centerInMu2e(),
                           parent.centerInWorld);

    ensu2Shield.solid = new G4ExtrudedSolid(ensu2Shield.name,
                                            ensu2Outline,
                                            ensu2->getLength()/2.,
                                            G4TwoVector(0,0), 1.,
                                            G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensu2Rotat(CLHEP::HepRotation::IDENTITY);
    ensu2Rotat.rotateX(ensu2->getRotation());

    finishNesting(ensu2Shield,
 		  findMaterialOrThrow(ensu2->materialName()),
 		  &ensu2Rotat,
 		  ensu2Shield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldUpstream2.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldUpstream2.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

  
    //----------------------------------------------------------------
    // ExtNeutShieldUpstreamTop <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensutopOutline; 
    std::copy(ensutop->externalShieldOutline().begin(),
	      ensutop->externalShieldOutline().end(),
	      std::back_inserter(ensutopOutline));

    VolumeInfo ensutopShield("ExtNeutShieldUpstreamTop",
                             ensutop->centerOfShield()-parent.centerInMu2e(),
                             parent.centerInWorld);

    ensutopShield.solid = new G4ExtrudedSolid(ensutopShield.name,
                                              ensutopOutline,
                                              ensutop->getLength()/2.,
                                              G4TwoVector(0,0), 1.,
                                              G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensutopRotat(CLHEP::HepRotation::IDENTITY);
    ensutopRotat.rotateX(ensutop->getRotation());

    finishNesting(ensutopShield,
 		  findMaterialOrThrow(ensutop->materialName()),
 		  &ensutopRotat,
 		  ensutopShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldUpstreamTop.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldUpstreamTop.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

  
    //----------------------------------------------------------------
    // ExtNeutShieldUpstreamBottom <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensubotOutline; 
    std::copy(ensubot->externalShieldOutline().begin(),
	      ensubot->externalShieldOutline().end(),
	      std::back_inserter(ensubotOutline));

    VolumeInfo ensubotShield("ExtNeutShieldUpstreamBottom",
                             ensubot->centerOfShield()-parent.centerInMu2e(),
                             parent.centerInWorld);

    ensubotShield.solid = new G4ExtrudedSolid(ensubotShield.name,
                                              ensubotOutline,
                                              ensubot->getLength()/2.,
                                              G4TwoVector(0,0), 1.,
                                              G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensubotRotat(CLHEP::HepRotation::IDENTITY);
    ensubotRotat.rotateX(ensubot->getRotation());

    finishNesting(ensubotShield,
 		  findMaterialOrThrow(ensubot->materialName()),
 		  &ensubotRotat,
 		  ensubotShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldUpstreamBottom.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldUpstreamBottom.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

  
    //----------------------------------------------------------------
    // ExtNeutShieldCavexRight <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> enscerOutline; 
    std::copy(enscer->externalShieldOutline().begin(),
	      enscer->externalShieldOutline().end(),
	      std::back_inserter(enscerOutline));

    VolumeInfo enscerShield("ExtNeutShieldCavexRight",
                            enscer->centerOfShield()-parent.centerInMu2e(),
                            parent.centerInWorld);

    enscerShield.solid = new G4ExtrudedSolid(enscerShield.name,
                                             enscerOutline,
                                             enscer->getLength()/2.,
                                             G4TwoVector(0,0), 1.,
                                             G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation enscerRotat(CLHEP::HepRotation::IDENTITY);
    enscerRotat.rotateX(enscer->getRotation());

    finishNesting(enscerShield,
 		  findMaterialOrThrow(enscer->materialName()),
 		  &enscerRotat,
 		  enscerShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldCavexRight.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldCavexRight.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // ExtNeutShieldCavexRightb <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> enscerbOutline; 
    std::copy(enscerb->externalShieldOutline().begin(),
	      enscerb->externalShieldOutline().end(),
	      std::back_inserter(enscerbOutline));

    VolumeInfo enscerbShield("ExtNeutShieldCavexRightb",
                             enscerb->centerOfShield()-parent.centerInMu2e(),
                             parent.centerInWorld);

    enscerbShield.solid = new G4ExtrudedSolid(enscerbShield.name,
                                              enscerbOutline,
                                              enscerb->getLength()/2.,
                                              G4TwoVector(0,0), 1.,
                                              G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation enscerbRotat(CLHEP::HepRotation::IDENTITY);
    enscerbRotat.rotateX(enscerb->getRotation());

    finishNesting(enscerbShield,
 		  findMaterialOrThrow(enscerb->materialName()),
 		  &enscerbRotat,
 		  enscerbShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldCavexRightb.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldCavexRightb.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // ExtNeutShieldCavexLeft <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> enscelOutline; 
    std::copy(enscel->externalShieldOutline().begin(),
	      enscel->externalShieldOutline().end(),
	      std::back_inserter(enscelOutline));

    VolumeInfo enscelShield("ExtNeutShieldCavexLeft",
                            enscel->centerOfShield()-parent.centerInMu2e(),
                            parent.centerInWorld);

    enscelShield.solid = new G4ExtrudedSolid(enscelShield.name,
                                             enscelOutline,
                                             enscel->getLength()/2.,
                                             G4TwoVector(0,0), 1.,
                                             G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation enscelRotat(CLHEP::HepRotation::IDENTITY);
    enscelRotat.rotateX(enscel->getRotation());

    finishNesting(enscelShield,
 		  findMaterialOrThrow(enscel->materialName()),
 		  &enscelRotat,
 		  enscelShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldCavexLeft.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldCavexLeft.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // ExtNeutShieldCavexRoof <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> ensceroOutline; 
    std::copy(enscero->externalShieldOutline().begin(),
	      enscero->externalShieldOutline().end(),
	      std::back_inserter(ensceroOutline));

    VolumeInfo ensceroShield("ExtNeutShieldCavexRoof",
                             enscero->centerOfShield()-parent.centerInMu2e(),
                             parent.centerInWorld);

    ensceroShield.solid = new G4ExtrudedSolid(ensceroShield.name,
                                              ensceroOutline,
                                              enscero->getLength()/2.,
                                              G4TwoVector(0,0), 1.,
                                              G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation ensceroRotat(CLHEP::HepRotation::IDENTITY);
    ensceroRotat.rotateX(enscero->getRotation());

    finishNesting(ensceroShield,
 		  findMaterialOrThrow(enscero->materialName()),
 		  &ensceroRotat,
 		  ensceroShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldCavexRoof.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldCavexRoof.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // ExtNeutShieldLAbove <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> enslaOutline; 
    std::copy(ensla->externalShieldOutline().begin(),
	      ensla->externalShieldOutline().end(),
	      std::back_inserter(enslaOutline));

    VolumeInfo enslaShield("ExtNeutShieldLAbove",
                           ensla->centerOfShield()-parent.centerInMu2e(),
                           parent.centerInWorld);

    enslaShield.solid = new G4ExtrudedSolid(enslaShield.name,
                                            enslaOutline,
                                            ensla->getLength()/2.,
                                            G4TwoVector(0,0), 1.,
                                            G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation enslaRotat(CLHEP::HepRotation::IDENTITY);
    enslaRotat.rotateX(ensla->getRotation());

    finishNesting(enslaShield,
 		  findMaterialOrThrow(ensla->materialName()),
 		  &enslaRotat,
 		  enslaShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldLAbove.visible"),
 		  G4Colour::Blue(), //G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldLAbove.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // ExtNeutShieldLCeiling <=====
    //----------------------------------------------------------------

    // convert internally held std::vector of Hep2Vectors and convert to
    // G4TwoVectors for G4

    std::vector<G4TwoVector> enslcOutline; 
    std::copy(enslc->externalShieldOutline().begin(),
	      enslc->externalShieldOutline().end(),
	      std::back_inserter(enslcOutline));

    VolumeInfo enslcShield("ExtNeutShieldLCeiling",
                           enslc->centerOfShield()-parent.centerInMu2e(),
                           parent.centerInWorld);

    enslcShield.solid = new G4ExtrudedSolid(enslcShield.name,
                                            enslcOutline,
                                            enslc->getLength()/2.,
                                            G4TwoVector(0,0), 1.,
                                            G4TwoVector(0,0), 1.);

    // rotation, needed for next call
    static CLHEP::HepRotation enslcRotat(CLHEP::HepRotation::IDENTITY);
    enslcRotat.rotateX(enslc->getRotation());

    finishNesting(enslcShield,
 		  findMaterialOrThrow(enslc->materialName()),
 		  &enslcRotat,
 		  enslcShield.centerInParent,
 		  parent.logical,
 		  0,
 		  config.getBool("ExtNeutShieldLCeiling.visible"),
 		  G4Colour::Magenta(),
 		  config.getBool("ExtNeutShieldLCeiling.solid"),
 		  forceAuxEdgeVisible,
 		  placePV,
 		  doSurfaceCheck);

    //----------------------------------------------------------------
    // DS-Cryo boxes <=====
    //----------------------------------------------------------------

    std::vector<std::vector<double>> dims = enscb->getDimensions();
    std::vector<std::string> mats = enscb->materialNames();
    std::vector<CLHEP::Hep3Vector> sites = enscb->centersOfBoxes();

    int nBox = dims.size();

    for(int i = 0; i < nBox; i++)
      {

	std::vector<double> lwhs = dims[i];

	//  Make the name of the box
	std::ostringstream name;
	name << "ExtNeutShieldCryoBox_" << i+1 ;

	// Make a non-rotating rotation
	static CLHEP::HepRotation fakeRotat(CLHEP::HepRotation::IDENTITY);

	// Build each box here

	nestBox( name.str(), lwhs, findMaterialOrThrow(mats[i]),
		 &fakeRotat, sites[i]-parent.centerInMu2e(),
		 parent.logical,
		 0,
		 config.getBool("ExtNeutShieldCryoBoxes.visible"),
		 G4Colour::Magenta(),
		 config.getBool("ExtNeutShieldCryoBoxes.solid"),
		 forceAuxEdgeVisible,
		 placePV,
		 doSurfaceCheck);

      }

    //----------------------------------------------------------------
    // DS-Cend boxes <=====
    //----------------------------------------------------------------

    std::vector<std::vector<double>> dimes = enscendb->getDimensions();
    std::vector<std::string> mates = enscendb->materialNames();
    std::vector<CLHEP::Hep3Vector> sitees = enscendb->centersOfBoxes();
    
    nBox = dimes.size();

    for(int i = 0; i < nBox; i++)
      {

	std::vector<double> lwhs = dimes[i];

	//  Make the name of the box
	std::ostringstream name;
	name << "ExtNeutShieldCendBox_" << i+1 ;

	// Make a non-rotating rotation
	static CLHEP::HepRotation fakeRotat(CLHEP::HepRotation::IDENTITY);

	// Build each box here

	if ( enscendb->hasHole(i) ) {

	  // This box has a window.  implemented as a
          // G4SubtractionSolid to alow for another volume placement
          // through it

 	  // Now make the window (AKA "Hole")
 	  name << "window";

          std::cout << " making " << name.str() << std::endl;

          std::ostringstream name1;
          name1 << "ExtNeutShieldCendBox_sub_" << i+1;

          G4Box* awindBox1 = new G4Box(name1.str(),lwhs[0],lwhs[1],lwhs[2]);

 	  // Find the index of this hole, for accessing info lists...
 	  int hID = enscendb->holeIndex(i);

 	  const TubsParams windparams(0.0,  //inner radius
 				      enscendb->holeRadius(hID), // outer
 				      enscendb->holeHalfLength(hID)  //obvious?
 				      );

          std::ostringstream name2;
          name2 << "ExtNeutShieldCendBox_sub_" << i+1 << "window";

          G4Tubs* awindTub = new G4Tubs( name2.str(), 
                                         windparams.data()[0], 
                                         windparams.data()[1], 
                                         windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                         windparams.data()[3], 
                                         windparams.data()[4]);

	  VolumeInfo awindBox;
          awindBox.name = name.str();
          
          // we need to put the window on the DS axis

          // fixme, this is specific to the Box #4 (but windows in any
          // other boxes will probably never be needed anyway)

          GeomHandle<DetectorSolenoid> ds;
          G4ThreeVector const & dsP ( ds->position() );

          G4ThreeVector offsetWRTDS(sitees[i].x()-dsP.x(), sitees[i].y()-dsP.y(), 0.0);

          awindBox.solid = new G4SubtractionSolid(awindBox.name,awindBox1,awindTub,0,-offsetWRTDS);

          finishNesting(awindBox,
                        findMaterialOrThrow(mates[i]),
                        0,
                        sitees[i]-parent.centerInMu2e(),
                        parent.logical,
                        0,
                        config.getBool("ExtNeutShieldCendBoxes.visible"),
                        G4Colour::Magenta(),
                        config.getBool("ExtNeutShieldCendBoxes.solid"),
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck);

	} else {

	  // Just put the box in the world, no window
	  nestBox( name.str(), lwhs, findMaterialOrThrow(mates[i]),
		   &fakeRotat, sitees[i]-parent.centerInMu2e(),
		   parent.logical,
		   0,
		   config.getBool("ExtNeutShieldCendBoxes.visible"),
		   G4Colour::Magenta(),
		   config.getBool("ExtNeutShieldCendBoxes.solid"),
		   forceAuxEdgeVisible,
		   placePV,
		   doSurfaceCheck);

	}

      }
  
  }

}
