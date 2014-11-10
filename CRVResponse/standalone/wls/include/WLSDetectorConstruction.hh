#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;

class G4LogicalVolume;
class G4VPhysicalVolume;

class WLSMaterials;
class WLSEventAction;
class G4Material;

#include "G4VUserDetectorConstruction.hh"

class WLSDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    WLSDetectorConstruction();
    ~WLSDetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    static WLSDetectorConstruction* Instance() {return _fgInstance;}

    double GetScintillatorHalfThickness() {return _scintillatorHalfThickness;}
    double GetScintillatorHalfWidth()     {return _scintillatorHalfWidth;}
    double GetScintillatorHalfLength()    {return _scintillatorHalfLength;}

    void UpdateGeometry();
 
    // Set Material Commands for World and WLSfiber
    void SetWorldMaterial         (G4String);
    void SetWLSFiberMaterial      (G4String);
    void SetCoupleMaterial        (G4String);

    G4double GetSurfaceRoughness();
    // Set the Roughness in between each layer
    void SetSurfaceRoughness      (G4double);
    // Set the reflectivity of the mirror
    void SetPhotonDetReflectivity (G4double);
    // Set the polish of the mirror
    void SetPhotonDetPolish       (G4double);

    void SetBarLength             (double barLength) {_barLength=barLength;}
    void SetBarWidth              (double barWidth) {_barWidth=barWidth;}
    void SetBarThickness          (double barThickness) {_barThickness=barThickness;}
    void SetFiberSeparation       (double fiberSeparation) {_fiberSeparation=fiberSeparation;}
    void SetHoleRadius            (double holeRadius) {_holeRadius=holeRadius;}
    void SetCoatingThickness      (double coatingThickness) {_coatingThickness=coatingThickness;}
    void SetFiberRadius           (double fiberRadius) {_fiberRadius=fiberRadius;}
    void SetClad1Radius           (double clad1Radius) {_clad1Radius=clad1Radius;}
    void SetClad2Radius           (double clad2Radius) {_clad2Radius=clad2Radius;}
    void SetSipmLength            (double sipmLength) {_sipmLength=sipmLength;}
    void SetSipmRadius            (double sipmRadius) {_sipmRadius=sipmRadius;}

    double GetBarLength()        {return _barLength;}
    double GetBarWidth()         {return _barWidth;}
    double GetBarThickness()     {return _barThickness;}
    double GetFiberSeparation()  {return _fiberSeparation;}
    double GetHoleRadius()       {return _holeRadius;}
    double GetCoatingThickness() {return _coatingThickness;}
    double GetFiberRadius()      {return _fiberRadius;}
    double GetClad1Radius()      {return _clad1Radius;}
    double GetClad2Radius()      {return _clad2Radius;}
    double GetSipmLength()       {return _sipmLength;}
    double GetSipmRadius()       {return _sipmRadius;}
 
    G4Material*         FindMaterial(G4String);
    G4VPhysicalVolume*  GetScintillatorVolume() {return physiScintillator;}

    std::vector<double> GetXBins() {return _xbins;}
    std::vector<double> GetYBins() {return _ybins;}
    std::vector<double> GetZBins() {return _zbins;}

  private:

    static WLSDetectorConstruction*  _fgInstance;  

    WLSMaterials* materials;

    G4LogicalVolume   *logicWorld, *logicHole;
    G4VPhysicalVolume *physiWorld, *physiHole1, *physiHole2;
    G4VPhysicalVolume *physiScintillator;
 
    double           _worldSizeX;
    double           _worldSizeY;
    double           _worldSizeZ;

    G4double mppcPolish;
    G4double mppcReflectivity;
    G4double extrusionPolish;
    G4double extrusionReflectivity;
    G4double surfaceRoughness;

    double _barLength, _barWidth, _barThickness;
    double _fiberSeparation;
    double _coatingThickness;
    double _holeRadius;
    double _fiberRadius, _clad1Radius, _clad2Radius;
    double _sipmRadius, _sipmLength;

    double _scintillatorHalfThickness;
    double _scintillatorHalfWidth;
    double _scintillatorHalfLength;

    std::vector<double> _xbins, _ybins, _zbins;

    void UpdateGeometryParameters();
};

#endif
