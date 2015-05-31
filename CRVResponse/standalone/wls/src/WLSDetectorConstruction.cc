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

WLSDetectorConstruction::WLSDetectorConstruction(int lengthOption) : physiWorld(NULL)
{
  _fgInstance = this;

  materials = NULL;

  surfaceRoughness = 1;
 
  mppcPolish = 1.;
  mppcReflectivity = 0.;

  extrusionPolish = 1.;
  extrusionReflectivity = 0.95;

  _barWidth         = 5.*cm;
  _barThickness     = 2.*cm;
  _fiberSeparation  = 2.*cm;
  _holeRadius       = 1.30*mm;
  _coatingThickness = 0.25*mm;
  _fiberRadius      = 0.70*mm - 0.021*mm - 0.021*mm;
  _clad1Radius      = 0.70*mm - 0.021*mm;
  _clad2Radius      = 0.70*mm;
  _sipmLength       = 1.*mm;
  _sipmRadius       = 0.70*mm;

  double xbinsTmp[17] = {-10.0, -7.5, -5.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5, 10.0};
  double ybinsTmp[36] = {-25.0, -20.0, -15.0, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -6.5, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 6.5, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 15.0, 20.0, 25.0};

  for(int i=0; i<17; i++) _xbins.push_back(xbinsTmp[i]*mm); //16 bins
  for(int i=0; i<36; i++) _ybins.push_back(ybinsTmp[i]*mm); //35 bins

  switch(lengthOption)
  {
    case 0: _barLength        = 660.*cm;
            //90 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-3300.0*mm+10.0*mm*i);       // -3300 ... -3250
            for(int i=6; i<16; i++)  _zbins.push_back(-3250.0*mm+25.0*mm*(i-5));   // -3225 ... -3000
            for(int i=16; i<75; i++) _zbins.push_back(-3000.0*mm+100.0*mm*(i-15)); // -2900 ...  2900
            for(int i=75; i<85; i++) _zbins.push_back(2975.0*mm+25.0*mm*(i-74));   //  3000 ...  3225
            for(int i=85; i<91; i++) _zbins.push_back(3240.0*mm+10.0*mm*(i-84));   //  3250 ...  3300
            break;
    case 1: _barLength        = 560.*cm;
            //80 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2800.0*mm+10.0*mm*i);       // -2800 ... -2750
            for(int i=6; i<16; i++)  _zbins.push_back(-2750.0*mm+25.0*mm*(i-5));   // -2725 ... -2500
            for(int i=16; i<65; i++) _zbins.push_back(-2500.0*mm+100.0*mm*(i-15)); // -2400 ...  2400
            for(int i=65; i<75; i++) _zbins.push_back(2475.0*mm+25.0*mm*(i-64));   //  2500 ...  2725
            for(int i=75; i<81; i++) _zbins.push_back(2740.0*mm+10.0*mm*(i-74));   //  2750 ...  2800
            break;
    case 2: _barLength        = 450.*cm;
            //69 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-2250.0*mm+10.0*mm*i);       // -2250 ... -2200
            for(int i=6; i<16; i++)  _zbins.push_back(-2200.0*mm+25.0*mm*(i-5));   // -2175 ... -1950
            for(int i=16; i<54; i++) _zbins.push_back(-1950.0*mm+100.0*mm*(i-15)); // -1850 ...  1850
            for(int i=54; i<64; i++) _zbins.push_back(1925.0*mm+25.0*mm*(i-53));   //  1950 ...  2175
            for(int i=64; i<70; i++) _zbins.push_back(2190.0*mm+10.0*mm*(i-63));   //  2200 ...  2250
            break;
    case 3: _barLength        = 300.*cm;  //test beam
            //54 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-1500.0*mm+10.0*mm*i);       // -1500 ... -1450
            for(int i=6; i<16; i++)  _zbins.push_back(-1450.0*mm+25.0*mm*(i-5));   // -1425 ... -1200
            for(int i=16; i<39; i++) _zbins.push_back(-1200.0*mm+100.0*mm*(i-15)); // -1100 ...  1100
            for(int i=39; i<49; i++) _zbins.push_back(1175.0*mm+25.0*mm*(i-38));   //  1200 ...  1425
            for(int i=49; i<55; i++) _zbins.push_back(1440.0*mm+10.0*mm*(i-48));   //  1450 ...  1500
            break;
    case 4: _barLength        = 230.*cm;
            //47 bins
            for(int i=0; i<6; i++)   _zbins.push_back(-1150.0*mm+10.0*mm*i);       // -1150 ... -1100
            for(int i=6; i<16; i++)  _zbins.push_back(-1100.0*mm+25.0*mm*(i-5));   // -1075 ...  -850
            for(int i=16; i<32; i++) _zbins.push_back(-850.0*mm+100.0*mm*(i-15));  //  -750 ...   750
            for(int i=32; i<42; i++) _zbins.push_back(825.0*mm+25.0*mm*(i-31));    //   850 ...  1075
            for(int i=42; i<48; i++) _zbins.push_back(1090.0*mm+10.0*mm*(i-41));   //  1100 ...  1150
            break;
    case 5: _barLength        = 90.*cm;
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
  if (materials)         delete materials;
}

G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
  materials = WLSMaterials::GetInstance();

  return ConstructDetector();
}

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
  //--------------------------------------------------
  // World
  //--------------------------------------------------

  G4VSolid* solidWorld =
                       new G4Box("World", _worldSizeX, _worldSizeY, _worldSizeZ);

  logicWorld = new G4LogicalVolume(solidWorld,
//                                   FindMaterial("G4_AIR"),
                                   FindMaterial("G4_POLYVINYL_CHLORIDE"),   //so that all photons at the front/back sides get absorbed
                                   "World");

  physiWorld = new G4PVPlacement(0,
                                 G4ThreeVector(),
                                 logicWorld,
                                 "World",
                                 0,
                                 false,
                                 0);

  //--------------------------------------------------
  // Extrusion
  //--------------------------------------------------

  G4VSolid* solidExtrusion =
        new G4Box("Extrusion",_barThickness/2.0,_barWidth/2.0,_barLength/2.0);

  G4LogicalVolume* logicExtrusion =
                      new G4LogicalVolume(solidExtrusion,
                                          FindMaterial("Coating"),
                                          "Extrusion");

  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       extrusionPolish);

  G4MaterialPropertiesTable* TiO2SurfaceProperty =
                                             new G4MaterialPropertiesTable();

  G4double p_TiO2[2] = {2.00*eV, 3.47*eV};
  G4double refl_TiO2[2] = {extrusionReflectivity,extrusionReflectivity};
  G4double effi_TiO2[2] = {0, 0};

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,2);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,2);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicExtrusion,
                    "Extrusion",
                    logicWorld,
                    false,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",logicExtrusion,TiO2Surface);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  G4VSolid* solidScintillator = new G4Box("Scintillator",
                                _scintillatorHalfThickness,
                                _scintillatorHalfWidth,
                                _scintillatorHalfLength);

  G4Material* polystyrene = FindMaterial("Polystyrene");
  G4LogicalVolume* logicScintillator =
                             new G4LogicalVolume(solidScintillator,
                                                 polystyrene,
                                                 "Scintillator");

  physiScintillator = new G4PVPlacement(0,
                                        G4ThreeVector(),
                                        logicScintillator,
                                        "Scintillator",
                                        logicExtrusion,
                                        false,
                                        0);

  G4MaterialPropertiesTable* materialPropertiesTable = polystyrene->GetMaterialPropertiesTable();
  double scintillationYield = materialPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");
  WLSSteppingAction::Instance()->SetScintillationYield(scintillationYield);
  double scintillatorDecayTimeFast = materialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
  double scintillatorDecayTimeSlow = materialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");
  WLSSteppingAction::Instance()->SetScintillatorDecayTimeFast(scintillatorDecayTimeFast);
  WLSSteppingAction::Instance()->SetScintillatorDecayTimeSlow(scintillatorDecayTimeSlow);

  //--------------------------------------------------
  // Fiber Holes
  //--------------------------------------------------

  G4VSolid* solidHole = new G4Tubs("Hole",
                                   0,
                                   _holeRadius,
                                   _barLength/2.0,
                                   0.*deg,
                                   360.*deg);
  logicHole = new G4LogicalVolume(solidHole,
                                  FindMaterial("G4_AIR"),
                                  "Hole");

  physiHole1 = new G4PVPlacement(0,
                                 G4ThreeVector(0.0, -_fiberSeparation/2.0, 0.0),
                                 logicHole,
                                 "Hole",
                                 logicScintillator,
                                 false,
                                 0);
  physiHole2 = new G4PVPlacement(0,
                                 G4ThreeVector(0.0, _fiberSeparation/2.0, 0.0),
                                 logicHole,
                                 "Hole",
                                 logicScintillator,
                                 false,
                                 1);
  //--------------------------------------------------
  // Fiber Construction
  //-------------------------------------------------- 

  // Boundary Surface Properties
  G4OpticalSurface* OpSurface = NULL;
 
  if (surfaceRoughness < 1.)
     OpSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      surfaceRoughness);       // SetPolish

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
                                                    false,
                                                    0);

  // Place the rough surface only if needed
  if (OpSurface) {
       new G4LogicalBorderSurface("surfaceFiber1Clad2Out",
                                  physiClad2,
                                  physiHole1,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceFiber1Clad2In",
                                  physiHole1,
                                  physiClad2,
                                  OpSurface);
  }
  if (OpSurface) {
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
                                                   false,
                                                   0);

  // Place the rough surface only if needed
  if (OpSurface) {
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

  G4Material* pmma = FindMaterial("PMMA");
  G4LogicalVolume* logicWLSfiber = new G4LogicalVolume(solidWLSfiber,
                                                       FindMaterial("PMMA"),
                                                       "WLSFiber");

  G4VPhysicalVolume* physiWLSfiber = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicWLSfiber,
                                                      "WLSFiber",
                                                      logicClad1,
                                                      false,
                                                      0);

  materialPropertiesTable = pmma->GetMaterialPropertiesTable();
  double fiberDecayTime = materialPropertiesTable->GetConstProperty("WLSTIMECONSTANT");
  WLSSteppingAction::Instance()->SetFiberDecayTime(fiberDecayTime);

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
  // PhotonDet (Sensitive Detector)
  //--------------------------------------------------  

  // Physical Construction
  G4VSolid* solidPhotonDet = new G4Tubs("PhotonDet",0.,_sipmRadius,_sipmLength/2.0,0.0*rad,twopi*rad);

  G4LogicalVolume*   logicPhotonDet =
                                    new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("PMMA"),  //TODO ???
                                                        "PhotonDet");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, -_barLength/2.0-_sipmLength/2.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicWorld,
                    false,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, -_fiberSeparation/2.0, _barLength/2.0+_sipmLength/2.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicWorld,
                    false,
                    1);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, -_barLength/2.0-_sipmLength/2.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicWorld,
                    false,
                    2);
  new G4PVPlacement(0,
                    G4ThreeVector(0.0, _fiberSeparation/2.0, _barLength/2.0+_sipmLength/2.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicWorld,
                    false,
                    3);

  // PhotonDet Surface Properties
  G4OpticalSurface* PhotonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       mppcPolish);

  G4MaterialPropertiesTable* PhotonDetSurfaceProperty =
                                               new G4MaterialPropertiesTable();

  G4double p_mppc[2] = {2.00*eV, 3.47*eV};
  G4double refl_mppc[2] = {mppcReflectivity,mppcReflectivity};
  G4double effi_mppc[2] = {1, 1};   //if a photon gets absorbed at a surface, 
                                    //it gets "detected" with the probability EFFICIENCY
                                    //(the status will be Detection) - see G4OpBoundaryProcess::DoAbsorption()
 
  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  //--------------------------------------------------
  // End of Construction
  //--------------------------------------------------

  return physiWorld;
}

void WLSDetectorConstruction::UpdateGeometry()
{
  if (!physiWorld) return;

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

  G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);
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
