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
    WLSDetectorConstruction();

  public:

    WLSDetectorConstruction(int lengthOption);
    ~WLSDetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    static WLSDetectorConstruction* Instance() {return _fgInstance;}

    double GetScintillatorHalfThickness() {return _scintillatorHalfThickness;}
    double GetScintillatorHalfWidth()     {return _scintillatorHalfWidth;}
    double GetScintillatorHalfLength()    {return _scintillatorHalfLength;}

    void UpdateGeometry();
 
    void SetBarLength             (double barLength) {_barLength=barLength;}
    void SetBarWidth              (double barWidth) {_barWidth=barWidth;}
    void SetBarThickness          (double barThickness) {_barThickness=barThickness;}
    void SetFiberSeparation       (double fiberSeparation) {_fiberSeparation=fiberSeparation;}
    void SetHoleRadiusX           (double holeRadius) {_holeRadiusX=holeRadius;}
    void SetHoleRadiusY           (double holeRadius) {_holeRadiusY=holeRadius;}
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
    double GetHoleRadiusX()      {return _holeRadiusX;}
    double GetHoleRadiusY()      {return _holeRadiusY;}
    double GetCoatingThickness() {return _coatingThickness;}
    double GetFiberRadius()      {return _fiberRadius;}
    double GetClad1Radius()      {return _clad1Radius;}
    double GetClad2Radius()      {return _clad2Radius;}
    double GetSipmLength()       {return _sipmLength;}
    double GetSipmRadius()       {return _sipmRadius;}
 
    G4Material*         FindMaterial(G4String);
    G4VPhysicalVolume*  GetScintillatorVolume() {return _physiScintillator;}

    std::vector<double> GetXBins() {return _xbins;}
    std::vector<double> GetYBins() {return _ybins;}
    std::vector<double> GetZBins() {return _zbins;}
    std::vector<double> GetBetaBins()  {return _betabins;}
    std::vector<double> GetThetaBins() {return _thetabins;}
    std::vector<double> GetPhiBins()   {return _phibins;}
    std::vector<double> GetRBins()     {return _rbins;}

  private:

    static WLSDetectorConstruction*  _fgInstance;  

    int  _lengthOption;
    bool _checkOverlaps;

    WLSMaterials* _materials;

    G4VPhysicalVolume *_physiWorld;
    G4VPhysicalVolume *_physiScintillator;
 
    double _worldSizeX;
    double _worldSizeY;
    double _worldSizeZ;

    double _mppcPolish;
    double _mppcReflectivity;
    double _mirrorPolish;
    double _mirrorReflectivity;
    double _extrusionPolish;
    double _fiberGuideBarPolish;
    double _fiberGuideBarReflectivity;
    double _holePolish;

    double _barLength, _barWidth, _barThickness;
    double _fiberSeparation;
    double _coatingThickness;
    double _holeRadiusX, _holeRadiusY;
    double _fiberRadius, _clad1Radius, _clad2Radius;
    double _sipmRadius, _sipmLength;

    double _scintillatorHalfThickness;
    double _scintillatorHalfWidth;
    double _scintillatorHalfLength;

    std::vector<double> _xbins, _ybins, _zbins;
    std::vector<double> _betabins, _thetabins, _phibins, _rbins;

    void UpdateGeometryParameters();
};

#endif
