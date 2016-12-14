#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

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
  _mppcReflectivity = 0.;

  _mirrorPolish = 1.;
  _mirrorReflectivity = 0.80;

  _extrusionPolish = 1.;
  _extrusionReflectivity = 0.90;

  _manifoldPolish = 1.;
  _manifoldReflectivity = 0.;

  _surfaceRoughness = 1;
 
  _barWidth         = 5.*cm;
  _barThickness     = 2.*cm;
  _fiberSeparation  = 2.6*cm;
  _holeRadius       = 1.30*mm;
  _coatingThickness = 0.25*mm;
  _fiberRadius      = 0.70*mm - 0.021*mm - 0.021*mm;
  _clad1Radius      = 0.70*mm - 0.021*mm;
  _clad2Radius      = 0.70*mm;
  _sipmLength       = 1.*mm;
  _sipmRadius       = 0.70*mm;

  double xbinsTmp[17] = {-10.0, -8.5, -5.5, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.5, 8.5, 10.0};
  double ybinsTmp[36] = {-25.0, -23.5, -22.0, -17.0, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -6.0, -4.5, -1.5, 1.5, 4.5, 6.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 17.0, 22.0, 23.5, 25.0};

  for(int i=0; i<17; i++) _xbins.push_back(xbinsTmp[i]*mm); //16 bins
  for(int i=0; i<36; i++) _ybins.push_back(ybinsTmp[i]*mm); //35 bins

  switch(_lengthOption)
  {
    case 7600: _barLength        = 760.*cm;
            //100 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3800.0*mm+10.0*mm*i);       // -3800 ... -3750
            for(int i=6; i<16; i++)  _zbins.push_back(-3750.0*mm+25.0*mm*(i-5));   // -3725 ... -3500
            for(int i=16; i<85; i++) _zbins.push_back(-3500.0*mm+100.0*mm*(i-15)); // -3400 ...  3400
            for(int i=85; i<95; i++) _zbins.push_back(3475.0*mm+25.0*mm*(i-84));   //  3500 ...  3725
            for(int i=95; i<101; i++) _zbins.push_back(3740.0*mm+10.0*mm*(i-94));  //  3750 ...  3800
            break;
    case 7100: _barLength        = 710.*cm;
            //95 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3550.0*mm+10.0*mm*i);       // -3550 ... -3500
            for(int i=6; i<16; i++)  _zbins.push_back(-3500.0*mm+25.0*mm*(i-5));   // -3475 ... -3250
            for(int i=16; i<80; i++) _zbins.push_back(-3250.0*mm+100.0*mm*(i-15)); // -3150 ...  3150
            for(int i=80; i<90; i++) _zbins.push_back(3225.0*mm+25.0*mm*(i-79));   //  3250 ...  3475
            for(int i=90; i<96; i++) _zbins.push_back(3490.0*mm+10.0*mm*(i-89));   //  3500 ...  3550
            break;
    case 6600: _barLength        = 660.*cm;
            //90 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3300.0*mm+10.0*mm*i);       // -3300 ... -3250
            for(int i=6; i<16; i++)  _zbins.push_back(-3250.0*mm+25.0*mm*(i-5));   // -3225 ... -3000
            for(int i=16; i<75; i++) _zbins.push_back(-3000.0*mm+100.0*mm*(i-15)); // -2900 ...  2900
            for(int i=75; i<85; i++) _zbins.push_back(2975.0*mm+25.0*mm*(i-74));   //  3000 ...  3225
            for(int i=85; i<91; i++) _zbins.push_back(3240.0*mm+10.0*mm*(i-84));   //  3250 ...  3300
            break;
    case 5600: _barLength        = 560.*cm;
            //80 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2800.0*mm+10.0*mm*i);       // -2800 ... -2750
            for(int i=6; i<16; i++)  _zbins.push_back(-2750.0*mm+25.0*mm*(i-5));   // -2725 ... -2500
            for(int i=16; i<65; i++) _zbins.push_back(-2500.0*mm+100.0*mm*(i-15)); // -2400 ...  2400
            for(int i=65; i<75; i++) _zbins.push_back(2475.0*mm+25.0*mm*(i-64));   //  2500 ...  2725
            for(int i=75; i<81; i++) _zbins.push_back(2740.0*mm+10.0*mm*(i-74));   //  2750 ...  2800
            break;
    case 4500: _barLength        = 450.*cm;
            //69 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2250.0*mm+10.0*mm*i);       // -2250 ... -2200
            for(int i=6; i<16; i++)  _zbins.push_back(-2200.0*mm+25.0*mm*(i-5));   // -2175 ... -1950
            for(int i=16; i<54; i++) _zbins.push_back(-1950.0*mm+100.0*mm*(i-15)); // -1850 ...  1850
            for(int i=54; i<64; i++) _zbins.push_back(1925.0*mm+25.0*mm*(i-53));   //  1950 ...  2175
            for(int i=64; i<70; i++) _zbins.push_back(2190.0*mm+10.0*mm*(i-63));   //  2200 ...  2250
            break;
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

  UpdateGeometryParameters();
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

  //--------------------------------------------------
  // Extrusion (TiO2 Coating)
  //--------------------------------------------------

  G4VSolid* solidExtrusion = new G4Box("Extrusion",_barThickness/2.0,_barWidth/2.0,_barLength/2.0);

  G4LogicalVolume* logicExtrusion = new G4LogicalVolume(solidExtrusion,
                                                        FindMaterial("Coating"),
                                                        "Extrusion");

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       _extrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_TiO2[11] =    {2.00*eV, 2.75*eV, 2.88*eV, 2.95*eV, 3.02*eV, 3.10*eV, 3.18*eV, 3.26*eV, 3.35*eV, 3.44*eV, 3.47*eV};
  G4double refl_TiO2[11] = {0.91,    0.91,    0.90,    0.85,    0.69,    0.44,    0.27,    0.13,    0.08,    0.07,    0.07};
  G4double effi_TiO2[11] = {0,       0 ,      0,       0,       0,       0,       0,       0,       0,       0,       0};

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,11);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,11);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicExtrusion,
                    "Extrusion",
                    logicWorld,
                    _checkOverlaps,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",logicExtrusion,TiO2Surface);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  G4VSolid* solidScintillator = new G4Box("Scintillator",
                                          _scintillatorHalfThickness,
                                          _scintillatorHalfWidth,
                                          _scintillatorHalfLength);

  G4LogicalVolume* logicScintillator = new G4LogicalVolume(solidScintillator,
                                                           FindMaterial("Polystyrene"),
                                                           "Scintillator");

  _physiScintillator = new G4PVPlacement(0,
                                         G4ThreeVector(),
                                         logicScintillator,
                                         "Scintillator",
                                         logicExtrusion,
                                         _checkOverlaps,
                                         0);


  //--------------------------------------------------
  // Plastic Manifold
  //--------------------------------------------------

  G4VSolid* solidManifold = new G4Box("Manifold",
                                      _scintillatorHalfThickness,
                                      _scintillatorHalfWidth,
                                      _sipmLength/2.0);

  G4LogicalVolume* logicManifold0 =  new G4LogicalVolume(solidManifold,
                                                         FindMaterial("G4_POLYVINYL_CHLORIDE"),
                                                         "Manifold0");
  G4LogicalVolume* logicManifold1 =  new G4LogicalVolume(solidManifold,                          //the nested volumes may be different (sipm vs. mirror)
                                                         FindMaterial("G4_POLYVINYL_CHLORIDE"),  //that's why we need to different logically volumes
                                                         "Manifold1");

  G4VPhysicalVolume *physiManifold0 = new G4PVPlacement(0,
                                                        G4ThreeVector(0.0, 0.0, -_barLength/2.0-_sipmLength/2.0),
                                                        logicManifold0,
                                                        "Manifold0",
                                                        logicWorld,
                                                        _checkOverlaps,
                                                        0);

  G4VPhysicalVolume *physiManifold1 = new G4PVPlacement(0,
                                                        G4ThreeVector(0.0, 0.0, _barLength/2.0+_sipmLength/2.0),
                                                        logicManifold1,
                                                        "Manifold1",
                                                        logicWorld,
                                                        _checkOverlaps,
                                                        0);


  //surface

  G4OpticalSurface* ManifoldSurface = new G4OpticalSurface("ManifoldSurface",
                                                           glisur,
                                                           ground,
                                                           dielectric_metal,
                                                           _manifoldPolish);

  G4MaterialPropertiesTable* ManifoldSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_Manifold[2] = {2.00*eV, 3.47*eV};
  G4double refl_Manifold[2] = {_manifoldReflectivity,_manifoldReflectivity};
  G4double effi_Manifold[2] = {0, 0};

  ManifoldSurfaceProperty -> AddProperty("REFLECTIVITY",p_Manifold,refl_Manifold,2);
  ManifoldSurfaceProperty -> AddProperty("EFFICIENCY",p_Manifold,effi_Manifold,2);

  ManifoldSurface -> SetMaterialPropertiesTable(ManifoldSurfaceProperty);

  new G4LogicalBorderSurface("surfaceScintillatorManifold0",
                             _physiScintillator,
                             physiManifold0,
                             ManifoldSurface);
  new G4LogicalBorderSurface("surfaceScintillatorManifold1",
                             _physiScintillator,
                             physiManifold1,
                             ManifoldSurface);


  //--------------------------------------------------
  // Fiber Holes
  //--------------------------------------------------

  G4VSolid* solidHole = new G4Tubs("Hole",
                                   0,
                                   _holeRadius,
                                   _barLength/2.0,
                                   0.*deg,
                                   360.*deg);
  G4LogicalVolume *logicHole = new G4LogicalVolume(solidHole,
                                                   FindMaterial("G4_AIR"),
                                                   "Hole");

  G4VPhysicalVolume *physiHole1 = new G4PVPlacement(0,
                                                    G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0),
                                                    logicHole,
                                                    "Hole",
                                                    logicScintillator,
                                                    _checkOverlaps,
                                                    0);
  G4VPhysicalVolume *physiHole2 = new G4PVPlacement(0,
                                                    G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0),
                                                    logicHole,
                                                    "Hole",
                                                    logicScintillator,
                                                    _checkOverlaps,
                                                    1);
  //--------------------------------------------------
  // Fiber Construction
  //-------------------------------------------------- 

  // Boundary Surface Properties
  G4OpticalSurface* OpSurface = NULL;
 
  if(_surfaceRoughness < 1.)
     OpSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      _surfaceRoughness);      // SetPolish

  //--------------------------------------------------
  // Cladding 2
  //--------------------------------------------------

  G4VSolid* solidClad2 = new G4Tubs("Clad2",0.,_clad2Radius,_barLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume* logicClad2  = new G4LogicalVolume(solidClad2,
                                                     FindMaterial("FPethylene"),
                                                     "Clad2");

  G4VPhysicalVolume* physiClad2 = new G4PVPlacement(0,
                                                    G4ThreeVector(),
                                                    logicClad2,
                                                    "Clad2",
                                                    logicHole,
                                                    _checkOverlaps,
                                                    0);

  // Place the rough surface only if needed
  if(OpSurface) 
  {
       new G4LogicalBorderSurface("surfaceFiber1Clad2Out",
                                  physiClad2,
                                  physiHole1,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceFiber1Clad2In",
                                  physiHole1,
                                  physiClad2,
                                  OpSurface);
  }
  if(OpSurface) 
  {
       new G4LogicalBorderSurface("surfaceFiber2Clad2Out",
                                  physiClad2,
                                  physiHole2,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceFiber2Clad2In",
                                  physiHole2,
                                  physiClad2,
                                  OpSurface);
  }

  //--------------------------------------------------
  // Cladding 1
  //--------------------------------------------------

  G4VSolid* solidClad1 = new G4Tubs("Clad1",0.,_clad1Radius,_barLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume* logicClad1 = new G4LogicalVolume(solidClad1,
                                                    FindMaterial("Pethylene"),
                                                    "Clad1");

  G4VPhysicalVolume* physiClad1 = new G4PVPlacement(0,
                                                   G4ThreeVector(),
                                                   logicClad1,
                                                   "Clad1",
                                                   logicClad2,
                                                   _checkOverlaps,
                                                   0);

  // Place the rough surface only if needed
  if(OpSurface) 
  {
       new G4LogicalBorderSurface("surfaceClad1Out",
                                  physiClad1,
                                  physiClad2,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceClad1In",
                                  physiClad2,
                                  physiClad1,
                                  OpSurface);
  }

  //--------------------------------------------------
  // WLS Fiber
  //--------------------------------------------------

  G4VSolid* solidWLSfiber = new G4Tubs("WLSFiber",0.,_fiberRadius,_barLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume* logicWLSfiber = new G4LogicalVolume(solidWLSfiber,
                                                       FindMaterial("PMMA"),
                                                       "WLSFiber");

  G4VPhysicalVolume* physiWLSfiber = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicWLSfiber,
                                                      "WLSFiber",
                                                      logicClad1,
                                                      _checkOverlaps,
                                                      0);

  // Place the rough surface only if needed
  if (OpSurface) {
       new G4LogicalBorderSurface("surfaceWLSOut",
                                  physiWLSfiber,
                                  physiClad1,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceWLSIn",
                                  physiClad1,
                                  physiWLSfiber,
                                  OpSurface);
  }

  //--------------------------------------------------
  // PhotonDet (Sensitive Detector) Or Mirror
  //--------------------------------------------------  

  // Physical Construction
  G4VSolid* solidPhotonDet = new G4Tubs("PhotonDet",0.,_sipmRadius,_sipmLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume* logicPhotonDet = new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("G4_Al"),
                                                        "PhotonDet");
  G4LogicalVolume* logicMirror = new G4LogicalVolume(solidPhotonDet,
                                                     FindMaterial("G4_Al"),
                                                     "Mirror");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicManifold0,
                    _checkOverlaps,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0),
                    _lengthOption>=6600?logicMirror:logicPhotonDet,
                    _lengthOption>=6600?"Mirror":"PhotonDet",
                    logicManifold1,
                    _checkOverlaps,
                    1);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicManifold0,
                    _checkOverlaps,
                    2);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0),
                    _lengthOption>=6600?logicMirror:logicPhotonDet,
                    _lengthOption>=6600?"Mirror":"PhotonDet",
                    logicManifold1,
                    _checkOverlaps,
                    3);

  // PhotonDet Surface Properties
  G4OpticalSurface* PhotonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                            glisur,
                                                            ground,
                                                            dielectric_metal,
                                                            _mppcPolish);

  G4MaterialPropertiesTable* PhotonDetSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_mppc[2] = {2.00*eV, 3.47*eV};
  G4double refl_mppc[2] = {_mppcReflectivity,_mppcReflectivity};
  G4double effi_mppc[2] = {1.0, 1.0};   //if a photon gets absorbed at a surface, 
                                        //it gets "detected" with the probability EFFICIENCY
                                        //(the status will be Detection) - see G4OpBoundaryProcess::DoAbsorption()

  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  // Mirror Surface Properties
  G4OpticalSurface* MirrorSurface = new G4OpticalSurface("PhotonDetSurface",
                                                         glisur,
                                                         ground,
                                                         dielectric_metal,
                                                         _mirrorPolish);

  G4MaterialPropertiesTable* MirrorSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_mirror[2] = {2.00*eV, 3.47*eV};
  G4double refl_mirror[2] = {_mirrorReflectivity,_mirrorReflectivity};
  G4double effi_mirror[2] = {0.0, 0.0};  
  
  MirrorSurfaceProperty -> AddProperty("REFLECTIVITY",p_mirror,refl_mirror,2);
  MirrorSurfaceProperty -> AddProperty("EFFICIENCY",p_mirror,effi_mirror,2);

  MirrorSurface -> SetMaterialPropertiesTable(MirrorSurfaceProperty);
 
  new G4LogicalSkinSurface("MirrorSurface",logicMirror,MirrorSurface); 

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
  UpdateGeometryParameters();
 
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  G4RegionStore::GetInstance()->UpdateMaterialList(_physiWorld);
}

void WLSDetectorConstruction::UpdateGeometryParameters()
{
  _worldSizeX = _barThickness + 1.*cm;
  _worldSizeY = _barWidth + 1.*cm;
  _worldSizeZ = _barLength + _sipmLength + 1.*cm;

  _scintillatorHalfThickness = _barThickness/2.0 - _coatingThickness;
  _scintillatorHalfWidth = _barWidth/2.0 - _coatingThickness;
  _scintillatorHalfLength = _barLength/2.0;
}

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) 
{
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
