//
// Free function to create the Vane calorimeter.
//
// $Id: constructVaneCalorimeter.cc,v 1.14 2014/08/01 23:14:36 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 23:14:36 $
//
// Original author Ivan Logashenko
// Modified by Bertrand Echenard
//
// Notes
//
//  1. a crystal has readouts at the back, both are surrounded by the wrapping, and the wrapping by a shell
//  2. by default, the wrapping surrounds the front/back face of the crystal+ro, the shell does not (shell is a casing)
//  3. The vanes are placed directly into DS3.  We did not make a mother volume for them.
//  4. The argument zOff is the zlocation of the center of the mother volume, as mesaured in the mu2e coordinate system.
//
//  5) Modified version  builds the calorimeter by making a physical mother volume.

// Important: The crystals need to be modelled as G4Polyhedra (not G4box) to be compatible with the disks

#include <iostream>

// Mu2e includes.
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4/inc/constructVaneCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructVaneCalorimeter( VolumeInfo const &  mother,SimpleConfig const& config )
  {


      int  const verbosityLevel              = config.getInt("calorimeter.verbosityLevel",0);
      bool const isCalorimeterVisible        = config.getBool("calorimeter.calorimeterVisible",false);
      bool const isCalorimeterSolid          = config.getBool("calorimeter.calorimeterSolid",false);
      bool const isVaneBoxVisible            = config.getBool("calorimeter.boxVisible",true);
      bool const isVaneBoxSolid              = config.getBool("calorimeter.boxSolid",true);
      bool const isAbsorberBoxVisible        = config.getBool("calorimeter.absorberBoxVisible",true);
      bool const isAbsorberBoxSolid          = config.getBool("calorimeter.absorberBoxSolid",true);
      bool const isCrystalVisible            = config.getBool("calorimeter.crystalVisible",false);
      bool const isCrystalSolid              = config.getBool("calorimeter.crystalSolid",true);
      bool const forceAuxEdgeVisible         = config.getBool("g4.forceAuxEdgeVisible",false);
      bool const doSurfaceCheck              = config.getBool("g4.doSurfaceCheck",false);


      //-- A helper class for parsing the config file.
      MaterialFinder materialFinder(config);
      G4Material* fillMaterial               = materialFinder.get("calorimeter.calorimeterFillMaterial");
      G4Material* crysMaterial               = materialFinder.get("calorimeter.crystalMaterial");
      G4Material* wrapMaterial               = materialFinder.get("calorimeter.crystalWrapper");
      G4Material* readMaterial               = materialFinder.get("calorimeter.crystalReadoutMaterial");
      G4Material* shieldMaterial             = materialFinder.get("calorimeter.shieldMaterial");
      G4Material* neutronAbsorberMaterial    = materialFinder.get("calorimeter.neutronAbsorberMaterial");

      G4VPhysicalVolume* pv;

      //-- Get calorimeter handle
      VaneCalorimeter const & cal = *(GeomHandle<VaneCalorimeter>());


      //calorimeter mother enveloppe
      G4double mother_inRadius               = cal.caloGeomInfo().enveloppeInRadius();
      G4double mother_outRadius              = cal.caloGeomInfo().enveloppeOutRadius();
      G4double mother_z0                     = cal.caloGeomInfo().enveloppeZ0();
      G4double mother_z1                     = cal.caloGeomInfo().enveloppeZ1();

      //crystal properties
      G4int    nRO                           = cal.caloGeomInfo().nROPerCrystal();
      G4double ROHalfThickness               = cal.caloGeomInfo().roHalfThickness();
      G4double ROHalfTrans                   = cal.caloGeomInfo().roHalfTrans();

      G4int    crystalnEdges                 = 4;
      G4double crystalPolysize               = cal.caloGeomInfo().crystalHalfTrans();
      G4double crystalDepth                  = 2.0*cal.caloGeomInfo().crystalHalfLength();


      G4double wrapThickness                 = cal.caloGeomInfo().wrapperThickness();
      G4double wrapPolysize                  = crystalPolysize + wrapThickness;
      G4double wrapDepth                     = crystalDepth + 2.0*ROHalfThickness + 2.0*wrapThickness; 

      const int nvane                        = cal.nVane();
      G4double  shieldHalfThickness          = config.getDouble("calorimeter.shieldHalfThickness");
      G4double  neutronAbsorberHalfThickness = config.getDouble("calorimeter.neutronAbsorberHalfThickness");
      G4double  absorberHalfThickness        = shieldHalfThickness + neutronAbsorberHalfThickness;


      bool       hasScint                    = config.getBool("calorimeter.hasScint", false);
      G4Material*scintMaterial               = materialFinder.get("calorimeter.scintMaterial", "G4_PLASTIC_SC_VINYLTOLUENE");
      double     scintDx                     = config.getDouble("calorimeter.scintDz", 20.);
      double     scintDy                     = config.getDouble("calorimeter.scintDz", 20.);
      double     scintDz                     = config.getDouble("calorimeter.scintDz", 10.);
      std::vector<double> scintPos;            config.getVectorDouble("calorimeter.scintPos", scintPos, 3);





      // Readout positions
      std::vector<double> XposRO, YposRO;

      double R0disp = 0.5*crystalPolysize;
      if (nRO==1) {XposRO.push_back(0);
                   YposRO.push_back(0);}
      if (nRO==2) {XposRO.push_back(0);XposRO.push_back(0); 
                   YposRO.push_back(-R0disp);YposRO.push_back(R0disp);}
      if (nRO==4) {XposRO.push_back(-R0disp);XposRO.push_back(-R0disp);XposRO.push_back(R0disp);XposRO.push_back(R0disp);
                   YposRO.push_back(-R0disp);YposRO.push_back(R0disp);YposRO.push_back(-R0disp);YposRO.push_back(R0disp);}



      // crystal z position
      G4double ZPoscrystal  = wrapThickness;
      G4double ZPosR0       = wrapThickness+crystalDepth+ROHalfThickness;



      //--------------------------------------
      // Building blocks for a crystal
      //

      //
      // define required solids


      double offsetAngle = CLHEP::pi/4.0;  //need to add an offset to the phi angle to have them rotated properly

      G4double crystalWrapZplanes[2] = {0,wrapDepth};
      G4double crystalWrapRinner[2]  = {0,0};
      G4double crystalWrapRouter[2]  = {wrapPolysize,wrapPolysize};
      G4Polyhedra* crystalWrap       = new G4Polyhedra("CrystalWrap",
                                           offsetAngle,CLHEP::twopi+offsetAngle, 
                                	   crystalnEdges,2, 
					   crystalWrapZplanes,crystalWrapRinner,crystalWrapRouter);

      G4double crystalZplanes[2] = {0,crystalDepth};
      G4double crystalRinner[2]  = {0,0};
      G4double crystalRouter[2]  = {crystalPolysize,crystalPolysize};
      G4Polyhedra* crystal       = new G4Polyhedra("Crystal",
                                       offsetAngle,CLHEP::twopi+offsetAngle, 
                                       crystalnEdges,2, 
				       crystalZplanes,crystalRinner,crystalRouter);

      G4Box *crystalRO = new G4Box("CrystalRO",ROHalfTrans,ROHalfTrans,ROHalfThickness);


      //
      // define required logical volumes

      G4LogicalVolume *CrystalLog  = new G4LogicalVolume(crystal,   crysMaterial, "CrystalLog");
      G4LogicalVolume *ROLog       = new G4LogicalVolume(crystalRO, readMaterial, "CrystalROLog" );    

      G4VisAttributes* crys_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Green());
      crys_visAtt->SetForceSolid(isCrystalSolid);
      crys_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);

      G4VisAttributes* ro_visAtt = new G4VisAttributes(isCrystalVisible, G4Color::Cyan());
      ro_visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);

      CrystalLog->SetVisAttributes(crys_visAtt);    
      ROLog->SetVisAttributes(ro_visAtt);

      //-- Sensitive detector
      G4VSensitiveDetector* ccSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal());
      G4VSensitiveDetector* crSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::CaloReadout());

      CrystalLog->SetSensitiveDetector(ccSD);
      ROLog->SetSensitiveDetector(crSD);




      //-- Construct calorrimeter mother volume

      double mother_zlength  = mother_z1-mother_z0;
      double mother_zCenter  = (mother_z1+mother_z0)/2.0;

      //  Make the mother volume for the calorimeter.
      CLHEP::Hep3Vector const& posDS3  = mother.centerInMu2e();
      G4ThreeVector posCaloMother      = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
      G4ThreeVector posCaloMotherInDS  = posCaloMother - posDS3;

      TubsParams caloParams(mother_inRadius,mother_outRadius,mother_zlength/2.0, 0., CLHEP::twopi);
      VolumeInfo calorimeterInfo = nestTubs( "CalorimeterMother",
					     caloParams,fillMaterial,0,posCaloMotherInDS,
					     mother,
					     0,isCalorimeterVisible,G4Colour::Blue(),isCalorimeterSolid,forceAuxEdgeVisible,
					     true,doSurfaceCheck);

      if ( verbosityLevel > 0) 
      {
	double zhl         = static_cast<G4Tubs*>(calorimeterInfo.solid)->GetZHalfLength();
	CLHEP::Hep3Vector const & CalorimeterOffsetInMu2e = calorimeterInfo.centerInMu2e();
	double CalorimeterOffsetInMu2eZ = CalorimeterOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " Calorimeter mother center in Mu2e   : " << CalorimeterOffsetInMu2e << endl;
	cout << __func__ << " Calorimeter mother Z extent in Mu2e    : " <<CalorimeterOffsetInMu2eZ - zhl << ", " << CalorimeterOffsetInMu2eZ + zhl << endl;
      }






      //--------------------------------------
      // Construct vznes: vaneInInfo contains the crystals. vaneOutInfo contains the absorbers/shields and the crystals
      //
      VolumeInfo vaneOutInfo[nvane];
      VolumeInfo vaneInInfo[nvane];
      VolumeInfo shieldInfo[nvane];
      VolumeInfo neutronAbsorberInfo[nvane];
      VolumeInfo scintInfo;


      for( int iv=0; iv<nvane; ++iv ) 
      {

	  ostringstream nameOutVane;           nameOutVane          << "CalorimeterOutVane_"            << iv;
	  ostringstream nameInVane;            nameInVane           << "CalorimeterInVane_"             << iv;
	  ostringstream nameShield;            nameShield           << "CalorimeterShield_"             << iv;
	  ostringstream nameNeutronAbsorber;   nameNeutronAbsorber  << "CalorimeterNeutronAbsorber_"    << iv;

	  double caseThickness = cal.caloGeomInfo().caseThickness();

	  const CLHEP::Hep3Vector & sizeOut = cal.vane(iv).size();
	  const CLHEP::Hep3Vector & sizeIn  = cal.vane(iv).size() - CLHEP::Hep3Vector(caseThickness,caseThickness,caseThickness);

	  double dimOutVane[3]             = {sizeOut.x(),                       sizeOut.y(), sizeOut.z()};
	  double dimInVane[3]              = {sizeIn.x()- absorberHalfThickness, sizeIn.y(),  sizeIn.z()};      
	  double dimShield[3]              = {shieldHalfThickness,               sizeIn.y(),  sizeIn.z()};
	  double dimNeutronAbsorber[3]     = {neutronAbsorberHalfThickness,      sizeIn.y(),  sizeIn.z()};

	  G4ThreeVector posVane            = cal.vane(iv).origin() - posCaloMother;
	  G4ThreeVector posInVane          = G4ThreeVector(absorberHalfThickness, 0. , 0.);
	  G4ThreeVector posShield          = G4ThreeVector(-sizeIn.x() + shieldHalfThickness, 0. , 0.);
	  G4ThreeVector posNeutronAbsorber = G4ThreeVector(-sizeIn.x() + 2*shieldHalfThickness + neutronAbsorberHalfThickness, 0. , 0.);

	  vaneOutInfo[iv]  = nestBox(nameOutVane.str(),
				     dimOutVane,fillMaterial,&cal.vane(iv).rotation(),posVane,
				     calorimeterInfo,
				     iv,
				     isVaneBoxVisible,
				     G4Colour::Yellow(),isVaneBoxSolid,forceAuxEdgeVisible,
				     true,doSurfaceCheck );


	  vaneInInfo[iv]   = nestBox(nameInVane.str(),
				     dimInVane,fillMaterial,0,posInVane,
				     vaneOutInfo[iv],
				     iv,
				     isVaneBoxVisible,G4Colour::Yellow(),isVaneBoxSolid,
				     forceAuxEdgeVisible,
				     true,doSurfaceCheck );


	  if( shieldHalfThickness > 0.0)     
	      shieldInfo[iv]  = nestBox(nameShield.str(),
					dimShield,shieldMaterial,0,posShield,
					vaneOutInfo[iv] ,
					iv,
					isAbsorberBoxVisible,G4Colour::Blue(),isAbsorberBoxSolid,forceAuxEdgeVisible,
					true,doSurfaceCheck );


	  if( neutronAbsorberHalfThickness > 0) 
               neutronAbsorberInfo[iv]  = nestBox(nameNeutronAbsorber.str(),
						  dimNeutronAbsorber,neutronAbsorberMaterial,0,posNeutronAbsorber,
						  vaneOutInfo[iv] ,
						  iv,
						  isAbsorberBoxVisible,G4Colour::Cyan(),isAbsorberBoxSolid,forceAuxEdgeVisible,
						  true,doSurfaceCheck );


	  if ( verbosityLevel > 0) 
	  {
	      double xhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetXHalfLength();
	      double yhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetYHalfLength();
	      double zhl  = static_cast<G4Box*>(vaneOutInfo[iv].solid)->GetZHalfLength();
	      cout << __func__ << " center in Mu2e    : " <<vaneOutInfo[iv].centerInMu2e() << endl;
	      cout << __func__ << " X extent in Mu2e (unrotated) : " <<vaneOutInfo[iv].centerInMu2e().x() - yhl << ", " << vaneOutInfo[iv].centerInMu2e().x() + yhl << endl;
	      cout << __func__ << " Y extent in Mu2e (unrotated) : " <<vaneOutInfo[iv].centerInMu2e().y() - zhl << ", " << vaneOutInfo[iv].centerInMu2e().y() + zhl << endl;
	      cout << __func__ << " Z extent in Mu2e (unrotated) : " <<vaneOutInfo[iv].centerInMu2e().z() - xhl << ", " << vaneOutInfo[iv].centerInMu2e().z() + xhl << endl;
	  }



	  //-- place crystals inside vanes  ---  see Note about rotation of polyhedra ---
	  G4int ncrys = cal.vane(iv).nCrystals();
	  for( int ic=0; ic<ncrys; ++ic ) 
	  {

	      // IDs
	      G4int id       = iv*ncrys + ic;       // Crystal ID
	      G4int roidBase = cal.ROBaseByCrystal(id);

	      //this is the position of the wrapper
	      CLHEP::Hep3Vector unitPosition = cal.vane(iv).crystal(ic).localPosition();
	      double x = unitPosition.x() - posInVane.x();  //local Position is in vaneOut, but the crystals are placed in VaneIn, so need this offset
	      double y = unitPosition.y();  
	      double z = -wrapDepth/2.0;

	      G4LogicalVolume *thisWrapLog = new G4LogicalVolume(crystalWrap, wrapMaterial, "WrapLog");
	      //thisWrapLog->SetVisAttributes(G4Colour::Red());
              thisWrapLog->SetVisAttributes(G4VisAttributes::Invisible);
	      pv = new G4PVPlacement(0,G4ThreeVector(x,y,z),thisWrapLog,"CrysWrapPV",vaneInInfo[iv].logical,false,id,false);   	      
              doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);


	      // -- place crystal inside warp -- 
	      pv = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,ZPoscrystal),CrystalLog,"CrysPV",thisWrapLog,false,id,false);
              doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);

              // add the readout
	      for (unsigned int iro=0;iro < XposRO.size();++iro)
	      {
		 pv = new G4PVPlacement(0,G4ThreeVector(XposRO[iro],YposRO[iro],ZPosR0),ROLog,"CrysROPV_0",thisWrapLog,true,roidBase+iro,false);
		 doSurfaceCheck && checkForOverlaps( pv, config, verbosityLevel>0);
	      }


	  }
      }
 
      //adding scintillators if needed for test beam
      if (hasScint)
      {
	  double posx = scintPos[0] - posCaloMother.x();
	  double posy = scintPos[1] - posCaloMother.y();
	  double posz = scintPos[2] - posCaloMother.z();
	  G4ThreeVector pos(posx, posy, posz);
	  double dimScint[3] = {scintDx, scintDy, scintDz};

	  scintInfo =  nestBox("ScintillatorFinger",
			       dimScint,scintMaterial,0,pos,
			       calorimeterInfo,
			       0,
			       true,G4Colour::Red(),isVaneBoxSolid,forceAuxEdgeVisible,
			       true,doSurfaceCheck );
      }

      return calorimeterInfo;

  }


} 
