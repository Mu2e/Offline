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

WLSDetectorConstruction::WLSDetectorConstruction(int lengthOption)
{
  _fgInstance = this;

  _lengthOption = lengthOption;
  _checkOverlaps = true;

  _materials = NULL;
  _physiWorld = NULL;

  _mppcPolish = 1.;
  _mppcReflectivity = 0.30;  

  _reflectorPolish = 1.;
  _reflectorReflectivity = 0.90;  //reflector
//  _reflectorReflectivity = 0.10;  //black tape

  _extrusionPolish = 0.3;

  _fiberGuideBarLength = 1.0*cm;
  _fiberGuideBarPolish = 1.;
  _fiberGuideBarReflectivity = 0.;

  _holePolish = 0.;
 
  _barWidth         = 5.*cm;
  _barThickness     = 2.*cm;
  _fiberSeparation  = 2.6*cm;
  _holeRadiusX      = 2.00*mm;
  _holeRadiusY      = 1.00*mm;
  _coatingThickness = 0.25*mm;
  _extrusionCornerRadius = 2.00*mm;
  _fiberRadius      = 0.70*mm - 0.021*mm - 0.021*mm;
  _clad1Radius      = 0.70*mm - 0.021*mm;
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

  _reflectorAtPositiveSide = false;
  _reflectorAtNegativeSide = false;

  double xbinsTmp[17] = {-10.0, -8.5, -5.5, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.5, 8.5, 10.0};
  double ybinsTmp[40] = {-25.0, -24.5, -24.0, -23.5, -20.5, -17.5, -15.5, -15.0, -14.5, -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -8.5, -5.5, -2.0, 2.0, 5.5, 8.5, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 17.5, 20.5, 23.5, 24.0, 24.5, 25.0};

  for(int i=0; i<17; i++) _xbins.push_back(xbinsTmp[i]*mm); //16 bins
  for(int i=0; i<40; i++) _ybins.push_back(ybinsTmp[i]*mm); //39 bins

  switch(_lengthOption)
  {
    case 7600: _barLength        = 760.*cm;
            //100 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3800.0*mm+10.0*mm*i);       // -3800 ... -3750
            for(int i=6; i<16; i++)  _zbins.push_back(-3750.0*mm+25.0*mm*(i-5));   // -3725 ... -3500
            for(int i=16; i<85; i++) _zbins.push_back(-3500.0*mm+100.0*mm*(i-15)); // -3400 ...  3400
            for(int i=85; i<95; i++) _zbins.push_back(3475.0*mm+25.0*mm*(i-84));   //  3500 ...  3725
            for(int i=95; i<101; i++) _zbins.push_back(3740.0*mm+10.0*mm*(i-94));  //  3750 ...  3800
            _reflectorAtPositiveSide = true;
            break;
    case 7100: _barLength        = 710.*cm;
            //95 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3550.0*mm+10.0*mm*i);       // -3550 ... -3500
            for(int i=6; i<16; i++)  _zbins.push_back(-3500.0*mm+25.0*mm*(i-5));   // -3475 ... -3250
            for(int i=16; i<80; i++) _zbins.push_back(-3250.0*mm+100.0*mm*(i-15)); // -3150 ...  3150
            for(int i=80; i<90; i++) _zbins.push_back(3225.0*mm+25.0*mm*(i-79));   //  3250 ...  3475
            for(int i=90; i<96; i++) _zbins.push_back(3490.0*mm+10.0*mm*(i-89));   //  3500 ...  3550
            _reflectorAtPositiveSide = true;
            break;
    case 6900: _barLength        = 690.*cm;
            //93 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3450.0*mm+10.0*mm*i);       // -3450 ... -3400
            for(int i=6; i<16; i++)  _zbins.push_back(-3400.0*mm+25.0*mm*(i-5));   // -3375 ... -3150
            for(int i=16; i<78; i++) _zbins.push_back(-3150.0*mm+100.0*mm*(i-15)); // -3050 ...  3050
            for(int i=78; i<88; i++) _zbins.push_back(3125.0*mm+25.0*mm*(i-77));   //  3150 ...  3375
            for(int i=88; i<94; i++) _zbins.push_back(3390.0*mm+10.0*mm*(i-87));   //  3400 ...  3450
            _reflectorAtPositiveSide = true;
            break;
    case 6600: _barLength        = 660.*cm;
            //90 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3300.0*mm+10.0*mm*i);       // -3300 ... -3250
            for(int i=6; i<16; i++)  _zbins.push_back(-3250.0*mm+25.0*mm*(i-5));   // -3225 ... -3000
            for(int i=16; i<75; i++) _zbins.push_back(-3000.0*mm+100.0*mm*(i-15)); // -2900 ...  2900
            for(int i=75; i<85; i++) _zbins.push_back(2975.0*mm+25.0*mm*(i-74));   //  3000 ...  3225
            for(int i=85; i<91; i++) _zbins.push_back(3240.0*mm+10.0*mm*(i-84));   //  3250 ...  3300
            _reflectorAtPositiveSide = true;
            break;
    case 6001:
            _reflectorAtPositiveSide = true;
    case 6000: _barLength        = 600.*cm;
            //84 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3000.0*mm+10.0*mm*i);       // -3000 ... -2950
            for(int i=6; i<16; i++)  _zbins.push_back(-2950.0*mm+25.0*mm*(i-5));   // -2925 ... -2700
            for(int i=16; i<69; i++) _zbins.push_back(-2700.0*mm+100.0*mm*(i-15)); // -2600 ...  2600
            for(int i=69; i<79; i++) _zbins.push_back(2675.0*mm+25.0*mm*(i-68));   //  2700 ...  2925
            for(int i=79; i<85; i++) _zbins.push_back(2940.0*mm+10.0*mm*(i-78));   //  2950 ...  3000
            break;
    case 5600: _barLength        = 560.*cm;
            //80 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2800.0*mm+10.0*mm*i);       // -2800 ... -2750
            for(int i=6; i<16; i++)  _zbins.push_back(-2750.0*mm+25.0*mm*(i-5));   // -2725 ... -2500
            for(int i=16; i<65; i++) _zbins.push_back(-2500.0*mm+100.0*mm*(i-15)); // -2400 ...  2400
            for(int i=65; i<75; i++) _zbins.push_back(2475.0*mm+25.0*mm*(i-64));   //  2500 ...  2725
            for(int i=75; i<81; i++) _zbins.push_back(2740.0*mm+10.0*mm*(i-74));   //  2750 ...  2800
            break;
    case 5000: _barLength        = 500.*cm;
            //74 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2500.0*mm+10.0*mm*i);       // -2500 ... -2450
            for(int i=6; i<16; i++)  _zbins.push_back(-2450.0*mm+25.0*mm*(i-5));   // -2425 ... -2200
            for(int i=16; i<59; i++) _zbins.push_back(-2200.0*mm+100.0*mm*(i-15)); // -2100 ...  2100
            for(int i=59; i<69; i++) _zbins.push_back(2175.0*mm+25.0*mm*(i-58));   //  2200 ...  2425
            for(int i=69; i<75; i++) _zbins.push_back(2440.0*mm+10.0*mm*(i-68));   //  2450 ...  2500
            _reflectorAtNegativeSide = true;
            break;
    case 4500: _barLength        = 450.*cm;
            //69 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2250.0*mm+10.0*mm*i);       // -2250 ... -2200
            for(int i=6; i<16; i++)  _zbins.push_back(-2200.0*mm+25.0*mm*(i-5));   // -2175 ... -1950
            for(int i=16; i<54; i++) _zbins.push_back(-1950.0*mm+100.0*mm*(i-15)); // -1850 ...  1850
            for(int i=54; i<64; i++) _zbins.push_back(1925.0*mm+25.0*mm*(i-53));   //  1950 ...  2175
            for(int i=64; i<70; i++) _zbins.push_back(2190.0*mm+10.0*mm*(i-63));   //  2200 ...  2250
            break;
    case 3001:
            _reflectorAtNegativeSide = true;
    case 3000: _barLength        = 300.*cm;  //test beam
            //54 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-1500.0*mm+10.0*mm*i);       // -1500 ... -1450
            for(int i=6; i<16; i++)  _zbins.push_back(-1450.0*mm+25.0*mm*(i-5));   // -1425 ... -1200
            for(int i=16; i<39; i++) _zbins.push_back(-1200.0*mm+100.0*mm*(i-15)); // -1100 ...  1100
            for(int i=39; i<49; i++) _zbins.push_back(1175.0*mm+25.0*mm*(i-38));   //  1200 ...  1425
            for(int i=49; i<55; i++) _zbins.push_back(1440.0*mm+10.0*mm*(i-48));   //  1450 ...  1500
            break;
    case 2300: _barLength        = 230.*cm;
            //47 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-1150.0*mm+10.0*mm*i);       // -1150 ... -1100
            for(int i=6; i<16; i++)  _zbins.push_back(-1100.0*mm+25.0*mm*(i-5));   // -1075 ...  -850
            for(int i=16; i<32; i++) _zbins.push_back(-850.0*mm+100.0*mm*(i-15));  //  -750 ...   750
            for(int i=32; i<42; i++) _zbins.push_back(825.0*mm+25.0*mm*(i-31));    //   850 ...  1075
            for(int i=42; i<48; i++) _zbins.push_back(1090.0*mm+10.0*mm*(i-41));   //  1100 ...  1150
            break;
    case 900: _barLength        = 90.*cm;
            //33 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-450.0*mm+10.0*mm*i);        //  -450 ...  -400
            for(int i=6; i<16; i++)  _zbins.push_back(-400.0*mm+25.0*mm*(i-5));    //  -375 ...  -150
            for(int i=16; i<18; i++) _zbins.push_back(-150.0*mm+100.0*mm*(i-15));  //   -50 ...    50
            for(int i=18; i<28; i++) _zbins.push_back(125.0*mm+25.0*mm*(i-17));    //   150 ...   375
            for(int i=28; i<34; i++) _zbins.push_back(390.0*mm+10.0*mm*(i-27));    //   400 ...   450
            break;
  }

  for(int i=0; i<5; i++)  _betabins.push_back(0.62+0.38*i/4.0); //4 bins

  _thetabins.push_back(0);                                                        //0
  for(int i=0; i<10; i++) _thetabins.push_back(CLHEP::pi/20.0+CLHEP::pi*i/10.0);  //1/20*pi ... 19/20*pi
  _thetabins.push_back(CLHEP::pi);                                                //pi    --> 11 bins

  for(int i=0; i<7; i++)  _phibins.push_back(-CLHEP::pi+i*CLHEP::pi/3.0); //6 bins
  for(int i=0; i<5; i++)  _rbins.push_back(_fiberRadius*i/4.0); //4 bins

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
                                                           FindMaterial("Polystyrene"),
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
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       _extrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_TiO2[11] =    {2.00*eV, 2.75*eV, 2.88*eV, 2.95*eV, 3.02*eV, 3.10*eV, 3.18*eV, 3.26*eV, 3.35*eV, 3.44*eV, 3.80*eV};
  G4double refl_TiO2[11] = {0.91,    0.91,    0.90,    0.85,    0.69,    0.44,    0.27,    0.13,    0.08,    0.07,    0.07}; 
  G4double effi_TiO2[11] = {0,       0 ,      0,       0,       0,       0,       0,       0,       0,       0,       0};
  for(int i=0; i<11; i++) refl_TiO2[i]=1.0-0.5*(1.0-refl_TiO2[i]);  //a higher reflectivities comparared to the numbers given by Anna 
                                                                    //improves the match between testbeam data and MC

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,11);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,11);

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

  G4double p_FiberGuideBar[2] = {2.00*eV, 3.80*eV};
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
                                                    FindMaterial("Pethylene"),
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
                                                       FindMaterial("PMMA"),
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

  G4double p_mppc[19] = {2.0*eV, 2.1*eV, 2.2*eV, 2.3*eV, 2.4*eV, 
                         2.5*eV, 2.6*eV, 2.7*eV, 2.8*eV, 2.9*eV, 
                         3.0*eV, 3.1*eV, 3.2*eV, 3.3*eV, 3.4*eV, 
                         3.5*eV, 3.6*eV, 3.7*eV, 3.8*eV};
  G4double refl_mppc[19];
  for(int i=0; i<19; i++) refl_mppc[i]=_mppcReflectivity;

  //if a photon gets absorbed at a surface, it gets "detected" with the probability EFFICIENCY
  //(the status will be Detection) - see G4OpBoundaryProcess::DoAbsorption()
  //the probability to produce a PE is done by the SiPM simulation
  G4double effi_mppc[19] = {0.602, 0.697, 0.784, 0.858, 0.933,
                            0.970, 0.990, 1.000, 0.978, 0.958,
                            0.920, 0.871, 0.803, 0.709, 0.622,
                            0.547, 0.398, 0.249, 0.100};

#ifdef FIBERTEST
  for(int i=0; i<19; i++) effi_mppc[i]=1.0;
#endif

  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,19);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,19);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  // Reflector Surface Properties
  G4OpticalSurface* ReflectorSurface = new G4OpticalSurface("ReflectorSurface",
                                                         glisur,
                                                         ground,
                                                         dielectric_metal,
                                                         _reflectorPolish);

  G4MaterialPropertiesTable* ReflectorSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_reflector[2] = {2.00*eV, 3.8*eV};
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
