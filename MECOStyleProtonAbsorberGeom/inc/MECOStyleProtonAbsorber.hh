#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh

//
// Class to represent the system of MECO Style Proton Absorber
//
//
// Original author MyeongJae Lee
//
// Coordinates are given in the detector coordinate
// system in mm.
//

// Includes from C++
#include <vector>
#include <memory>

// Includes from Mu2e
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorberPart.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ProtonAbsorberId {
    public:
      enum  enum_type {
        unknown = -1,
        pabs1 = 0,
        pabs2 = 1,
        opabs1 = 2,
        opabs2 = 3};
      explicit ProtonAbsorberId (enum_type id):
        _id(id)
      {}
      enum_type id() const { return _id; }
      ProtonAbsorberId():
        _id(unknown)
      {}
    private:
      enum_type _id;
  };

  class InnerProtonAbsSupport {

    friend class MECOStyleProtonAbsorberMaker;

  public:

    InnerProtonAbsSupport( std::size_t nSets, std::size_t nWiresPerSet )
      : _nSets ( nSets ) , _nWiresPerSet ( nWiresPerSet ) {}

    void setWireAngleOffset(double offset) { _wireAngleOffset = offset; }
    const Tube& getWire( std::size_t iSet, std::size_t iWire ) const { return _supportWireMap.at(iSet).at(iWire); }
    std::size_t nSets() const { return _nSets; }
    std::size_t nWiresPerSet() const { return _nWiresPerSet; }
    double      wireAngleOffset() const { return _wireAngleOffset; }
    const Tube& getEndRing( std::size_t iRing ) const { return _endRingMap.at(iRing); }
    std::size_t nEndRings() const { return _nEndRings; }

  private:
    std::size_t _nSets;
    std::size_t _nWiresPerSet;
    double      _wireAngleOffset = 0;
    std::vector<std::vector<Tube>> _supportWireMap;
    std::size_t _nEndRings = 0;
    std::vector<Tube> _endRingMap;
  };

  class MECOStyleProtonAbsorber : virtual public Detector{

  friend class MECOStyleProtonAbsorberMaker;

  public:
    MECOStyleProtonAbsorber() ;

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    MECOStyleProtonAbsorberPart const& part ( unsigned int n ) const { return _parts.at(n); }
    double virtualDetectorHalfLength()  const { return _vdHL; }
    std::string fillMaterial()           const { return _materialName; }
    double distanceFromTargetEnd() const { return _distfromtargetend; }
    double halfLength () const { return _halflength; }
    double thickness() const { return _thickness;}
    bool isAvailable (int id) const ;

    // outer PA
    std::string outerPAfillMaterial()    const {return _oPAmaterialName; }
    double outerPAzcenter ()             const { return _oPAzcenter; }
    double outerPAhalflength ()          const { return _oPAhalflength; }
    double outerPAthickness ()           const { return _oPAthickness; }
    double slotWidth()                   const { return _oPAslotWidth; }
    double slotLength()                  const { return _oPAslotLength; }
    double slotOffset()                  const { return _oPAslotOffset; }
    int    outerPAVersion()              const { return _oPAVersion; }

    // Outer PA supports
    int oPAnSupports()                             const { return _oPAnSupports; }
    std::vector<double> oPAsupportInnerRadii()     const { return _oPASupportIR; }
    std::vector<double> oPAsupportOuterRadii()     const { return _oPASupportOR; }
    std::vector<double> oPAsupportHalflength()     const { return _oPASupportHL; }
    std::vector<double> oPAsupportZMidpoint()      const { return _oPASupportZM; }
    std::vector<double> oPAsupportExtra()          const { return _oPASupportHE; }
    std::vector<double> oPAsupportXRad()           const { return _oPASupportXR; }
    std::vector<double> oPAsupportDPhiX()          const { return _oPASupportPH; }
    std::vector<std::string> oPAsupportMaterials() const { return _oPASupportMaterials; }
    //stopping target supports
    int                 nOPASupportSlats()         const { return _nOPASupportSlats;}
    std::vector<int>    slatTypes()                const { return _oPASlatTypes; }
    int                 nOPASupportSlatTypes()     const { return _nOPASupportSlatTypes;}
    std::vector<double> slatHeights()              const { return _oPASlatHeights; }
    std::vector<double> slatWidths()               const { return _oPASlatWidths; }
    std::vector<double> slatLengths()              const { return _oPASlatLengths; }
    std::vector<double> slatAngles()               const { return _oPASlatAngles; }
    std::vector<double> slatSideThicknesses()      const { return _oPASlatSideThicknesses; }
    std::vector<double> slatTopThicknesses()       const { return _oPASlatTopThicknesses; }
    std::vector<double> slatFillParameter1()       const { return _oPASlatFillParameter1; }
    std::vector<double> slatFillParameter2()       const { return _oPASlatFillParameter2; }
    std::vector<double> slatFillParameter3()       const { return _oPASlatFillParameter3; }
    std::vector<double> slatFillParameter4()       const { return _oPASlatFillParameter4; }
    std::vector<std::string> slatMaterials()       const { return _oPASlatMaterials;}
    std::vector<std::string> slatFillMaterials()   const { return _oPASlatFillMaterials;}
    //cross supports between support rings
    int                 nCrossSupports         () const { return _nCrossSupports         ;}
    std::vector<double> crossSupportThicknesses() const { return _crossSupportThicknesses;}
    std::vector<double> crossSupportWidth      () const { return _crossSupportWidth      ;}
    std::vector<int>    crossSupportOneIndex   () const { return _crossSupportOneIndex   ;}
    std::vector<int>    crossSupportTwoIndex   () const { return _crossSupportTwoIndex   ;}
    std::vector<double> crossSupportPhis       () const { return _crossSupportPhis       ;}
    std::vector<double> crossSupportHeights    () const { return _crossSupportHeights    ;}
    std::vector<double> crossSupportRadii      () const { return _crossSupportRadii      ;}
    std::string         crossSupportMaterial   () const { return _crossSupportMaterial   ;}

    // support structure for inner PA
    bool buildSupports() const { return _buildSupports; }
    const InnerProtonAbsSupport* getIPAsupport() const { return _ipaSupport.get(); }

    // Pion Degrader
    bool        degraderBuild()    const { return _degraderBuild; }
    double      degraderRotation() const { return _degraderRot; }
    double      degraderZ0()       const { return _degraderZ0; }
    std::string degraderFilterMaterial() const { return _degraderFiltMaterial;}
    std::string degraderFrameMaterial()  const { return _degraderFramMaterial;}
    std::string degraderCountwtMaterial() const {return _degraderCowtMaterial;}
    std::string degraderRodMaterial() const { return  _degraderRodMaterial; }
    std::string degraderSupportMaterial() const {return _degraderSuptMaterial;}
    std::vector<double> degraderFrameDims() const {return  _degraderFrameDims;}
    std::vector<double> degraderFilterDims() const {return  _degraderFilterDims;}
    std::vector<double> degraderCounterwtDims() const { return  _degraderCounterDims; }
    std::vector<double> degraderRodDims() const { return _degraderRodDims;}
    std::vector<double> degraderPivotPos() const { return _degraderPivotPos;}
    std::vector<double> degraderSupportArmDims() const { return _degraderSupportArmDims;}
    std::vector<double> degraderSupportPlateDims() const { return _degraderSupportPlateDims;}

  protected:

    std::vector<MECOStyleProtonAbsorberPart> _parts;
    // some variables that affects both parts
    double _vdHL;        // Virtual Detector half length
    std::string _materialName;  // Proton Absorber material
    double _distfromtargetend;  //distance from the target end to the start of proton absorber
    double _halflength;
    double _thickness;
    bool _pabs1flag, _pabs2flag;

    //outer PA
    std::string _oPAmaterialName;
    double _oPAzcenter;
    double _oPAhalflength;
    double _oPAthickness;
    bool _oPA1flag, _oPA2flag;
    double _oPAslotWidth;  // width of slots in OPA for ST support wires
    double _oPAslotLength;  // length of slots in OPA for ST support wires
    double _oPAslotOffset;  // offset of slots in OPA relative to DS2 part center
    int    _oPAVersion;

    int                      _oPAnSupports;  // How many supports there are
    std::vector<double>      _oPASupportIR;  // Inner radii
    std::vector<double>      _oPASupportOR;  // Outer radii
    std::vector<double>      _oPASupportHL;  // half lengths
    std::vector<double>      _oPASupportZM;  // mid-point Z values
    std::vector<double>      _oPASupportHE;  // Has Extra on bottom (only for upstm)
    std::vector<double>      _oPASupportXR;  // Extra radial amount
    std::vector<double>      _oPASupportPH;  // dPhi of extra bit
    std::vector<std::string> _oPASupportMaterials;  // material for each support

    int                      _nOPASupportSlats;
    std::vector<int>         _oPASlatTypes;
    int                      _nOPASupportSlatTypes;
    std::vector<double>      _oPASlatHeights;
    std::vector<double>      _oPASlatWidths;
    std::vector<double>      _oPASlatLengths;
    std::vector<double>      _oPASlatAngles;
    std::vector<double>      _oPASlatSideThicknesses;
    std::vector<double>      _oPASlatTopThicknesses;
    std::vector<std::string> _oPASlatMaterials;
    //parameters for adding fill to ST support slats
    std::vector<double>      _oPASlatFillParameter1;
    std::vector<double>      _oPASlatFillParameter2;
    std::vector<double>      _oPASlatFillParameter3;
    std::vector<double>      _oPASlatFillParameter4;
    std::vector<std::string> _oPASlatFillMaterials;

    int                 _nCrossSupports         ;
    std::vector<double> _crossSupportThicknesses;
    std::vector<double> _crossSupportWidth      ;
    std::vector<int>    _crossSupportOneIndex   ;
    std::vector<int>    _crossSupportTwoIndex   ;
    std::vector<double> _crossSupportPhis       ;
    std::vector<double> _crossSupportHeights    ;
    std::vector<double> _crossSupportRadii      ;
    std::string         _crossSupportMaterial   ;


    // support structure for inner PA
    bool _buildSupports;
    std::unique_ptr<InnerProtonAbsSupport> _ipaSupport;

    // Info for Pion degrader
    bool                 _degraderBuild;
    double               _degraderRot;
    double               _degraderZ0;
    std::string          _degraderFiltMaterial;
    std::string          _degraderFramMaterial;
    std::string          _degraderCowtMaterial;
    std::string          _degraderRodMaterial;
    std::string          _degraderSuptMaterial;
    std::vector<double>  _degraderFrameDims;
    std::vector<double>  _degraderFilterDims;
    std::vector<double>  _degraderCounterDims;
    std::vector<double>  _degraderRodDims;
    std::vector<double>  _degraderPivotPos;
    std::vector<double>  _degraderSupportArmDims;
    std::vector<double>  _degraderSupportPlateDims;
  };
}
#endif
