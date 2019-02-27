#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"

#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSMaterials.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

WLSDetectorConstruction* WLSDetectorConstruction::_fgInstance = NULL;

WLSDetectorConstruction::WLSDetectorConstruction(double lengthOption, int reflectorOption)
{
  _fgInstance = this;

  _checkOverlaps = true;

  _materials = NULL;
  _physiWorld = NULL;

  _mppcPolish = 1.;
  _mppcReflectivity = 0.30;  

  _reflectorPolish = 1.;
  _reflectorReflectivity = 0.90;  //reflector at test beam
//  _reflectorReflectivity = 0.10;  //black tape at test beam

  _extrusionPolish = 0.3;

  _fiberGuideBarLength = 1.0*cm;
  _fiberGuideBarPolish = 1.;
  _fiberGuideBarReflectivity = 0.;

  _holePolish = 0.;
 
  _barLength        = lengthOption*mm;
//#define TESTBEAM
#ifndef TESTBEAM
  _barThickness     = 19.78*mm;
  _barWidth         = 51.34*mm;
#else
#pragma message "USING TESTBEAM"
  _barThickness     = 19.8*mm;  //FIXME testbeam
  _barWidth         = 49.4*mm;  //FIXME testbeam
#endif
  _fiberSeparation  = 2.6*cm;
  _holeRadiusX      = 2.00*mm;
  _holeRadiusY      = 1.00*mm;
  _coatingThickness = 0.25*mm;
  _extrusionCornerRadius = 2.00*mm;
  _fiberRadius      = 0.70*mm - 0.042*mm - 0.042*mm;
  _clad1Radius      = 0.70*mm - 0.042*mm;
  _clad2Radius      = 0.70*mm;
  _sipmLength       = 1.*mm;
  _sipmRadius       = 0.70*mm;
  _airGap           = 0.0*mm;

#ifdef FIBERTEST
#pragma message "USING FIBERTEST"
  _mppcReflectivity       = 0.0;
  _reflectorReflectivity  = 0.0;
  _fiberGuideBarLength = 0.0*cm;
  _airGap              = 0.0*mm;
#endif

  _reflectorOption = reflectorOption;
  _reflectorAtPositiveSide = (reflectorOption==1?true:false);
  _reflectorAtNegativeSide = (reflectorOption==-1?true:false);

//use only the positive values of x (assuming symmetry) (0 is at the center)
#ifndef TESTBEAM
  double xbinsTmp[11] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 6.0, 9.0, 9.5, 9.89};
  double ybinsTmp[42] = {-25.67, -25.2, -24.7, -23.5, -20.5, -17.5, -15.5, -15.0, -14.5, -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -8.5, -6.5, -4.5, -1.5, 1.5, 4.5, 6.5, 8.5, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 17.5, 20.5, 23.5, 24.7, 25.2, 25.67}; 
#else
#pragma message "USING TESTBEAM"
  double xbinsTmp[11] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 6.0, 9.0, 9.5, 9.9}; //FIXME testbeam
  double ybinsTmp[42] = {-24.7, -24.25, -23.75, -22.5, -20.5, -17.5, -15.5, -15.0, -14.5, -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -8.5, -6.5, -4.5, -1.5, 1.5, 4.5, 6.5, 8.5, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 17.5, 20.5, 22.5, 23.75, 24.25, 24.7};  //FIXME testbeam
#endif

  for(int i=0; i<11; i++) _xbins.push_back(xbinsTmp[i]*mm); //10 bins
  for(int i=0; i<42; i++) _ybins.push_back(ybinsTmp[i]*mm); //41 bins

  double halfLength = _barLength/2.0;
  int    nBinsMainSection = lrint((_barLength-600.0)/100.0);
  double binWidthMainSection = 0;
  if(nBinsMainSection>0) binWidthMainSection = (_barLength-600.0)/nBinsMainSection;
  for(int i=0; i<=5;  i++)               _zbins.push_back(-halfLength*mm+                        10.0*mm*i);
  for(int i=1; i<=10; i++)               _zbins.push_back(-halfLength*mm+ 50.0*mm+               25.0*mm*i);
  for(int i=1; i<=nBinsMainSection; i++) _zbins.push_back(-halfLength*mm+300.0*mm+binWidthMainSection*mm*i);
  for(int i=1; i<=10; i++)               _zbins.push_back( halfLength*mm-300.0*mm+               25.0*mm*i);
  for(int i=1; i<=5;  i++)               _zbins.push_back( halfLength*mm- 50.0*mm+               10.0*mm*i);

  double inverseBetaBins[6] = {2.20, 1.80, 1.64, 1.58, 1.24, 1.0};     //values of interest from the index of refraction of polystyrene
  for(int i=0; i<6; i++)  _betabins.push_back(1.0/inverseBetaBins[i]); //5 bins

  _thetabins.push_back(0);                                                        //0
  for(int i=0; i<10; i++) _thetabins.push_back(CLHEP::pi/20.0+CLHEP::pi*i/10.0);  //1/20*pi ... 19/20*pi
  _thetabins.push_back(CLHEP::pi);                                                //pi    --> 11 bins

  _phibins.push_back(0);                                                      //0
  for(int i=0; i<4; i++)  _phibins.push_back(CLHEP::pi/8.0+CLHEP::pi*i/4.0);  //1/8*pi ... 7/8*pi
  _phibins.push_back(CLHEP::pi);                                              //pi    --> 5 bins
                                                                              //don't need pi...2*pi due to symmetry

  for(int i=0; i<6; i++)  _rbins.push_back(_fiberRadius*i/5.0); //5 bins

  _worldSizeX = _barThickness + 1.*cm;
  _worldSizeY = _barWidth + 1.*cm;
  _worldSizeZ = _barLength + _fiberGuideBarLength+_airGap+_sipmLength + 1.*cm;
}

WLSDetectorConstruction::~WLSDetectorConstruction()
{
  if(_materials) delete _materials;
}

G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
  _materials = WLSMaterials::GetInstance();

  return ConstructDetector();
}

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
  //--------------------------------------------------
  // World
  //--------------------------------------------------

  G4VSolid* solidWorld = new G4Box("World", _worldSizeX, _worldSizeY, _worldSizeZ);

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    FindMaterial("G4_AIR"),
                                                    "World");

  _physiWorld = new G4PVPlacement(0,
                                  G4ThreeVector(),
                                  logicWorld,
                                  "World",
                                  0,
                                  _checkOverlaps,
                                  0);

#ifndef FIBERTEST
  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  G4VSolid* solidScintillatorBase = new G4Box("ScintillatorBase",
                                          _barThickness/2.0-_coatingThickness,
                                          _barWidth/2.0-_coatingThickness,
                                          _barLength/2.0);

  double scintillatorCornerRadius = _extrusionCornerRadius-_coatingThickness;
  G4VSolid* solidScintillatorCorner1 = new G4Tubs("ScintillatorCorner1", scintillatorCornerRadius, 1.5*scintillatorCornerRadius, _barLength/2.0, 0.0, halfpi);
  G4VSolid* solidScintillatorCorner2 = new G4Tubs("ScintillatorCorner2", scintillatorCornerRadius, 1.5*scintillatorCornerRadius, _barLength/2.0, halfpi, halfpi);
  G4VSolid* solidScintillatorCorner3 = new G4Tubs("ScintillatorCorner3", scintillatorCornerRadius, 1.5*scintillatorCornerRadius, _barLength/2.0, pi, halfpi);
  G4VSolid* solidScintillatorCorner4 = new G4Tubs("ScintillatorCorner4", scintillatorCornerRadius, 1.5*scintillatorCornerRadius, _barLength/2.0, 3.0*halfpi, halfpi);

  double absScintCornerPosX = _barThickness/2.0-_coatingThickness-scintillatorCornerRadius;
  double absScintCornerPosY = _barWidth/2.0-_coatingThickness-scintillatorCornerRadius;
  G4SubtractionSolid *solidScintillator1Corner = new G4SubtractionSolid("Scintillator1Corner",solidScintillatorBase, solidScintillatorCorner1, NULL,
                                                                        G4ThreeVector(absScintCornerPosX, absScintCornerPosY, 0.0));
  G4SubtractionSolid *solidScintillator2Corner = new G4SubtractionSolid("Scintillator2Corner",solidScintillator1Corner, solidScintillatorCorner2, NULL,
                                                                        G4ThreeVector(-absScintCornerPosX, absScintCornerPosY, 0.0));
  G4SubtractionSolid *solidScintillator3Corner = new G4SubtractionSolid("Scintillator3Corner",solidScintillator2Corner, solidScintillatorCorner3, NULL,
                                                                        G4ThreeVector(-absScintCornerPosX, -absScintCornerPosY, 0.0));
  G4SubtractionSolid *solidScintillator4Corner = new G4SubtractionSolid("Scintillator4Corner",solidScintillator3Corner, solidScintillatorCorner4, NULL,
                                                                        G4ThreeVector(absScintCornerPosX, -absScintCornerPosY, 0.0));

  G4VSolid* solidScintillatorHole = new G4EllipticalTube("ScintillatorHole",
                                                          _holeRadiusX,
                                                          _holeRadiusY,
                                                          _barLength/2.0);

  G4SubtractionSolid * solidScintillator1Hole = new G4SubtractionSolid("Scintillator1Hole", 
                                                                        solidScintillator4Corner, solidScintillatorHole, NULL, 
                                                                        G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0));
  G4SubtractionSolid * solidScintillator2Hole = new G4SubtractionSolid("Scintillator2Hole", 
                                                                        solidScintillator1Hole, solidScintillatorHole, NULL, 
                                                                        G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0));

  G4LogicalVolume* logicScintillator = new G4LogicalVolume(solidScintillator2Hole,
                                                           FindMaterial("PolystyreneScint"),
                                                           "Scintillator");

  _physiScintillator = new G4PVPlacement(0,
                                         G4ThreeVector(),
                                         logicScintillator,
                                         "Scintillator",
                                         logicWorld,
                                         _checkOverlaps,
                                         0);
  _physiScintillator->CheckOverlaps();

  //--------------------------------------------------
  // Extrusion (TiO2 Coating)
  //--------------------------------------------------

  G4VSolid* solidExtrusionBase = new G4Box("ExtrusionBase",_barThickness/2.0,_barWidth/2.0,_barLength/2.0);

  G4VSolid* solidExtrusionCorner1 = new G4Tubs("ExtrusionCorner1", _extrusionCornerRadius, 1.5*_extrusionCornerRadius, _barLength/2.0, 0.0, halfpi);
  G4VSolid* solidExtrusionCorner2 = new G4Tubs("ExtrusionCorner2", _extrusionCornerRadius, 1.5*_extrusionCornerRadius, _barLength/2.0, halfpi, halfpi);
  G4VSolid* solidExtrusionCorner3 = new G4Tubs("ExtrusionCorner3", _extrusionCornerRadius, 1.5*_extrusionCornerRadius, _barLength/2.0, pi, halfpi);
  G4VSolid* solidExtrusionCorner4 = new G4Tubs("ExtrusionCorner4", _extrusionCornerRadius, 1.5*_extrusionCornerRadius, _barLength/2.0, 3.0*halfpi, halfpi);

  double absExtrCornerPosX = _barThickness/2.0-_extrusionCornerRadius;
  double absExtrCornerPosY = _barWidth/2.0-_extrusionCornerRadius;
  G4SubtractionSolid *solidExtrusion1Corner = new G4SubtractionSolid("Extrusion1Corner",solidExtrusionBase, solidExtrusionCorner1, NULL,
                                                                        G4ThreeVector(absExtrCornerPosX, absExtrCornerPosY, 0.0));
  G4SubtractionSolid *solidExtrusion2Corner = new G4SubtractionSolid("Extrusion2Corner",solidExtrusion1Corner, solidExtrusionCorner2, NULL,
                                                                        G4ThreeVector(-absExtrCornerPosX, absExtrCornerPosY, 0.0));
  G4SubtractionSolid *solidExtrusion3Corner = new G4SubtractionSolid("Extrusion3Corner",solidExtrusion2Corner, solidExtrusionCorner3, NULL,
                                                                        G4ThreeVector(-absExtrCornerPosX, -absExtrCornerPosY, 0.0));
  G4SubtractionSolid *solidExtrusion4Corner = new G4SubtractionSolid("Extrusion4Corner",solidExtrusion3Corner, solidExtrusionCorner4, NULL,
                                                                        G4ThreeVector(absExtrCornerPosX, -absExtrCornerPosY, 0.0));

  G4SubtractionSolid * solidExtrusion = new G4SubtractionSolid("Extrusion", 
                                                               solidExtrusion4Corner, solidScintillator4Corner, NULL, 
                                                               G4ThreeVector(0.0, 0.0, 0.0));

  G4LogicalVolume* logicExtrusion = new G4LogicalVolume(solidExtrusion,
                                                        FindMaterial("Coating"),
                                                        "Extrusion");

  G4VPhysicalVolume *physiExtrusion = new G4PVPlacement(0,
                                                        G4ThreeVector(),
                                                        logicExtrusion,
                                                        "Extrusion",
                                                        logicWorld,
                                                        _checkOverlaps,
                                                        0);
  physiExtrusion->CheckOverlaps();

  //--------------------------------------------------
  // Extrusion (TiO2 Coating) Surface
  //--------------------------------------------------

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       unified,
                                                       ground,
                                                       dielectric_metal,
                                                       1.5);

  G4MaterialPropertiesTable* TiO2SurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_TiO2[11] =    {2.00*eV, 2.75*eV, 2.88*eV, 2.95*eV, 3.02*eV, 3.10*eV, 3.18*eV, 3.26*eV, 3.35*eV, 3.44*eV, 15.75*eV};
  G4double refl_TiO2[11] = {0.91,    0.91,    0.90,    0.85,    0.69,    0.44,    0.27,    0.13,    0.08,    0.07,    0.07}; //assume a constant value for energies > 3.44eV (most photons with these energies get absorbed and wave length shifted)
  G4double effi_TiO2[11] = {0,       0 ,      0,       0,       0,       0,       0,       0,       0,       0,       0};
  for(int i=0; i<11; i++) refl_TiO2[i]=1.0-0.9*(1.0-refl_TiO2[i]);  //a higher reflectivities comparared to the numbers given by Anna 
                                                                    //improves the match between testbeam data and MC

  TiO2SurfaceProperty->AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,11);
  TiO2SurfaceProperty->AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,11);

  G4double pp[2] = {2.00*eV, 15.75*eV};
  G4double specularlobe[2] = {1.0, 1.0};
  G4double specularspike[2] = {0.0, 0.0};
  G4double backscatter[2] = {0.0, 0.0};
  TiO2SurfaceProperty->AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,2);
  TiO2SurfaceProperty->AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,2);
  TiO2SurfaceProperty->AddProperty("BACKSCATTERCONSTANT",pp,backscatter,2);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4LogicalBorderSurface("TiO2Surface",  
                             _physiScintillator,
                             physiExtrusion,
                             TiO2Surface);

  //--------------------------------------------------
  // Plastic Fiber Guide Bar
  //--------------------------------------------------

  G4VSolid* solidFiberGuideBarBase = new G4Box("FiberGuideBarBase",
                                      _barThickness/2.0,
                                      _barWidth/2.0,
                                      _fiberGuideBarLength/2.0);

  G4VSolid* solidFiberGuideBarHole = new G4Tubs("FiberGuideBarHole",0.,_clad2Radius,_fiberGuideBarLength/2.0,0.0*rad,twopi*rad);

  G4SubtractionSolid * solidFiberGuideBar1Hole = new G4SubtractionSolid("FiberGuideBar1Hole", 
                                                                        solidFiberGuideBarBase, solidFiberGuideBarHole, NULL, 
                                                                        G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0));
  G4SubtractionSolid * solidFiberGuideBar2Hole = new G4SubtractionSolid("FiberGuideBar2Hole", 
                                                                        solidFiberGuideBar1Hole, solidFiberGuideBarHole, NULL, 
                                                                        G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0));

  G4LogicalVolume* logicFiberGuideBar =  new G4LogicalVolume(solidFiberGuideBar2Hole,
                                                             FindMaterial("G4_POLYVINYL_CHLORIDE"),
                                                             "FiberGuideBar");

  G4VPhysicalVolume *physiFiberGuideBar0=
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0,-_barLength/2.0-_fiberGuideBarLength/2.0),
                                                        logicFiberGuideBar,
                                                        "FiberGuideBar",
                                                        logicWorld,
                                                        _checkOverlaps,
                                                        0);

  G4VPhysicalVolume *physiFiberGuideBar1=
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, _barLength/2.0+_fiberGuideBarLength/2.0),
                                                        logicFiberGuideBar,
                                                        "FiberGuideBar",
                                                        logicWorld,
                                                        _checkOverlaps,
                                                        1);
  physiFiberGuideBar0->CheckOverlaps();
  physiFiberGuideBar1->CheckOverlaps();


  //surface

  G4OpticalSurface* FiberGuideBarSurface = new G4OpticalSurface("FiberGuideBarSurface",
                                                           glisur,
                                                           ground,
                                                           dielectric_metal,
                                                           _fiberGuideBarPolish);

  G4MaterialPropertiesTable* FiberGuideBarSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_FiberGuideBar[2] = {2.00*eV, 15.75*eV};
  G4double refl_FiberGuideBar[2] = {_fiberGuideBarReflectivity,_fiberGuideBarReflectivity};
  G4double effi_FiberGuideBar[2] = {0, 0};

  FiberGuideBarSurfaceProperty -> AddProperty("REFLECTIVITY",p_FiberGuideBar,refl_FiberGuideBar,2);
  FiberGuideBarSurfaceProperty -> AddProperty("EFFICIENCY",p_FiberGuideBar,effi_FiberGuideBar,2);

  FiberGuideBarSurface -> SetMaterialPropertiesTable(FiberGuideBarSurfaceProperty);

  new G4LogicalSkinSurface("FiberGuideBarSurface",logicFiberGuideBar,FiberGuideBarSurface); 

#endif

  //--------------------------------------------------
  // Cladding 2
  //--------------------------------------------------

  G4VSolid* solidClad2 = new G4Tubs("Clad2",0.,_clad2Radius,_barLength/2.0+_fiberGuideBarLength,0.0*rad,twopi*rad);

  G4LogicalVolume* logicClad2  = new G4LogicalVolume(solidClad2,
                                                     FindMaterial("FPethylene"),
                                                     "Clad2");

  G4VPhysicalVolume *physiClad2Hole0=
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0),
                    logicClad2,
                    "Clad2Hole0",
                    logicWorld,
                    _checkOverlaps,
                    0);
  G4VPhysicalVolume *physiClad2Hole1=
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0),
                    logicClad2,
                    "Clad2Hole1",
                    logicWorld,
                    _checkOverlaps,
                    1);

  physiClad2Hole0->CheckOverlaps();
  physiClad2Hole1->CheckOverlaps();

  //--------------------------------------------------
  // Cladding 1
  //--------------------------------------------------

  G4VSolid* solidClad1 = new G4Tubs("Clad1",0.,_clad1Radius,_barLength/2.0+_fiberGuideBarLength,0.0*rad,twopi*rad);

  G4LogicalVolume* logicClad1 = new G4LogicalVolume(solidClad1,
                                                    FindMaterial("PMMA"),
                                                    "Clad1");

  G4VPhysicalVolume *physiClad1=
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicClad1,
                    "Clad1",
                    logicClad2,
                    _checkOverlaps,
                    0);
  physiClad1->CheckOverlaps();

  //--------------------------------------------------
  // WLS Fiber
  //--------------------------------------------------

  G4VSolid* solidWLSfiber = new G4Tubs("WLSFiber",0.,_fiberRadius,_barLength/2.0+_fiberGuideBarLength,0.0*rad,twopi*rad);

  G4LogicalVolume* logicWLSfiber = new G4LogicalVolume(solidWLSfiber,
                                                       FindMaterial("PolystyreneFiber"),
                                                       "WLSFiber");

  G4VPhysicalVolume *physiWLSfiber=
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicWLSfiber,
                    "WLSFiber",
                    logicClad1,
                    _checkOverlaps,
                    0);
  physiWLSfiber->CheckOverlaps();

  //--------------------------------------------------
  // PhotonDet (Sensitive Detector) Or Reflector
  //--------------------------------------------------  

  // Physical Construction
  G4VSolid* solidPhotonDet = new G4Tubs("PhotonDet",0.,_sipmRadius,_sipmLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume* logicPhotonDet = new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("G4_Si"),
                                                        "PhotonDet");
  G4LogicalVolume* logicReflector = new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("G4_Al"),
                                                        "Reflector");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, -_barLength/2.0-_fiberGuideBarLength-_airGap-_sipmLength/2.0),
                    _reflectorAtNegativeSide?logicReflector:logicPhotonDet,
                    _reflectorAtNegativeSide?"Reflector":"PhotonDet",
                    logicWorld,
                    _checkOverlaps,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, _barLength/2.0+_fiberGuideBarLength+_airGap+_sipmLength/2.0),
                    _reflectorAtPositiveSide?logicReflector:logicPhotonDet,
                    _reflectorAtPositiveSide?"Reflector":"PhotonDet",
                    logicWorld,
                    _checkOverlaps,
                    1);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, -_barLength/2.0-_fiberGuideBarLength-_airGap-_sipmLength/2.0),
                    _reflectorAtNegativeSide?logicReflector:logicPhotonDet,
                    _reflectorAtNegativeSide?"Reflector":"PhotonDet",
                    logicWorld,
                    _checkOverlaps,
                    2);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, _barLength/2.0+_fiberGuideBarLength+_airGap+_sipmLength/2.0),
                    _reflectorAtPositiveSide?logicReflector:logicPhotonDet,
                    _reflectorAtPositiveSide?"Reflector":"PhotonDet",
                    logicWorld,
                    _checkOverlaps,
                    3);

  // PhotonDet Surface Properties
  G4OpticalSurface* PhotonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                            glisur,
                                                            ground,
                                                            dielectric_metal,
                                                            _mppcPolish);

  G4MaterialPropertiesTable* PhotonDetSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_mppc[21] = {2.0*eV, 2.1*eV, 2.2*eV, 2.3*eV, 2.4*eV, 
                         2.5*eV, 2.6*eV, 2.7*eV, 2.8*eV, 2.9*eV, 
                         3.0*eV, 3.1*eV, 3.2*eV, 3.3*eV, 3.4*eV, 
                         3.5*eV, 3.6*eV, 3.7*eV, 3.8*eV, 3.9*eV,
                       15.75*eV};
  G4double refl_mppc[21];
  for(int i=0; i<21; i++) refl_mppc[i]=_mppcReflectivity;

  //if a photon gets absorbed at a surface, it gets "detected" with the probability EFFICIENCY
  //(the status will be Detection) - see G4OpBoundaryProcess::DoAbsorption()
  //the probability to produce a PE is done by the SiPM simulation
  G4double effi_mppc[21] = {0.602, 0.697, 0.784, 0.858, 0.933,
                            0.970, 0.990, 1.000, 0.978, 0.958,
                            0.920, 0.871, 0.803, 0.709, 0.622,
                            0.547, 0.398, 0.249, 0.100, 0.000,
                            0.000};

#ifdef FIBERTEST
  for(int i=0; i<21; i++) effi_mppc[i]=1.0;
#endif

  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,21);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,21);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  // Reflector Surface Properties
  G4OpticalSurface* ReflectorSurface = new G4OpticalSurface("ReflectorSurface",
                                                         glisur,
                                                         ground,
                                                         dielectric_metal,
                                                         _reflectorPolish);

  G4MaterialPropertiesTable* ReflectorSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_reflector[2] = {2.00*eV, 15.75*eV};
  G4double refl_reflector[2] = {_reflectorReflectivity,_reflectorReflectivity};
  G4double effi_reflector[2] = {0.0, 0.0};  
  
  ReflectorSurfaceProperty -> AddProperty("REFLECTIVITY",p_reflector,refl_reflector,2);
  ReflectorSurfaceProperty -> AddProperty("EFFICIENCY",p_reflector,effi_reflector,2);

  ReflectorSurface -> SetMaterialPropertiesTable(ReflectorSurfaceProperty);
 
  new G4LogicalSkinSurface("ReflectorSurface",logicReflector,ReflectorSurface); 

  //--------------------------------------------------
  // End of Construction
  //--------------------------------------------------

  return _physiWorld;
}

void WLSDetectorConstruction::UpdateGeometry()
{
  if (!_physiWorld) return;

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  G4RegionStore::GetInstance()->UpdateMaterialList(_physiWorld);
}

void WLSDetectorConstruction::UpdateGeometryParameters()
{
  std::cout<<"DOES ANYONE CALL THIS???"<<std::endl;
}

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) 
{
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
