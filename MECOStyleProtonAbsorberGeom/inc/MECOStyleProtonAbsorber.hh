#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorber_hh

//
// Class to represent the system of MECO Style Proton Absorber
//
// $Id: MECOStyleProtonAbsorber.hh,v 1.4 2013/06/19 03:41:01 mjlee Exp $
// $Author: mjlee $
// $Date: 2013/06/19 03:41:01 $
//
// Original author MyeongJae Lee
//
// Coordinates are given in the detector coordinate
// system in mm.
//

// Includes from C++
#include <vector>

// Includes from Mu2e
#include "Mu2eInterfaces/inc/Detector.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorberPart.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

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
    bool _oPAflag;
  };
}
#endif 
