#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh

//
// Class to represent the system of MECO Style Proton Absorber
//
// $Id: MECOStyleProtonAbsorber.hh,v 1.6 2014/02/28 21:11:19 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/28 21:11:19 $
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
#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eInterfaces/inc/Detector.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorberPart.hh"

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
    ~InnerProtonAbsSupport(){}

    const Tube& getWire( std::size_t iSet, std::size_t iWire ) const { return _supportWireMap.at(iSet).at(iWire); }
    std::size_t nSets() const { return _nSets; }
    std::size_t nWiresPerSet() const { return _nWiresPerSet; }
    const Tube& getEndRing( std::size_t iRing ) const { return _endRingMap.at(iRing); }
    std::size_t nEndRings() const { return _nEndRings; }

  private:
    std::size_t _nSets;
    std::size_t _nWiresPerSet;
    std::vector<std::vector<Tube>> _supportWireMap;
    std::size_t _nEndRings;
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
    std::string outerPAfillMaterial() const {return _oPAmaterialName; }
    double outerPAzcenter () const { return _oPAzcenter; }
    double outerPAhalflength () const { return _oPAhalflength; }
    double outerPAthickness () const { return _oPAthickness; }

    // support structure for inner PA
    bool buildSupports() const { return _buildSupports; }
    const InnerProtonAbsSupport* getIPAsupport() const { return _ipaSupport.get(); }

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

    // support structure for inner PA
    bool _buildSupports;
    std::unique_ptr<InnerProtonAbsSupport> _ipaSupport;

  };
}
#endif 
