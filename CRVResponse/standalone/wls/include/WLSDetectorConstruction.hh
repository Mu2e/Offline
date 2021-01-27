#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

//#define FIBERTEST

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

    WLSDetectorConstruction(double lengthOption, int reflectorOption);
    ~WLSDetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    static WLSDetectorConstruction* Instance() {return _fgInstance;}

    void UpdateGeometry();
 
    G4Material*         FindMaterial(G4String);
    G4VPhysicalVolume*  GetScintillatorVolume() {return _physiScintillator;}

    double GetScintillatorHalfThickness() {return _barThickness/2.0-_coatingThickness;}
    double GetScintillatorHalfWidth()     {return _barWidth/2.0-_coatingThickness;}
    double GetScintillatorHalfLength()    {return _barLength/2.0;}
    double GetScintillatorCornerRadius()  {return _extrusionCornerRadius-_coatingThickness;}
    double GetFiberSeparation()  {return _fiberSeparation;}
    double GetHoleRadiusX()      {return _holeRadiusX;}
    double GetHoleRadiusY()      {return _holeRadiusY;}
    double GetClad2Radius()      {return _clad2Radius;}

    int    GetReflectorOption()  {return _reflectorOption;}

    std::vector<double> GetXBins() {return _xbins;}
    std::vector<double> GetYBins() {return _ybins;}
    std::vector<double> GetZBins() {return _zbins;}
    std::vector<double> GetBetaBins()  {return _betabins;}
    std::vector<double> GetThetaBins() {return _thetabins;}
    std::vector<double> GetPhiBins()   {return _phibins;}
    std::vector<double> GetRBins()     {return _rbins;}

  private:

    static WLSDetectorConstruction*  _fgInstance;  

    WLSMaterials* _materials;

    G4VPhysicalVolume *_physiWorld;
    G4VPhysicalVolume *_physiScintillator;
 
    double _worldSizeX;
    double _worldSizeY;
    double _worldSizeZ;

    double _mppcReflectivity;
    double _blacktapeReflectivity;
    double _reflectorReflectivity;
    double _fiberGuideBarReflectivity;

    double _barLength, _barWidth, _barThickness;
    double _fiberSeparation;
    double _coatingThickness;
    double _extrusionCornerRadius;
    double _holeRadiusX, _holeRadiusY;
    double _fiberRadius, _clad1Radius, _clad2Radius;
    double _sipmWidth, _sipmLength, _sipmWindowLength;
    double _airGap;
    double _fiberGuideBarLength;
  
    int    _reflectorOption;
    bool   _reflectorAtPositiveSide;
    bool   _reflectorAtNegativeSide;

    std::vector<double> _xbins, _ybins, _zbins;
    std::vector<double> _betabins, _thetabins, _phibins, _rbins;

    void UpdateGeometryParameters();
};

#endif
