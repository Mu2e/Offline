#ifndef TargetGeom_Target_hh
#define TargetGeom_Target_hh

//
// Class to represent the system of target foils.
// For now these are just disks perpendicular to the z axis.
//
// $Id: Target.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate 
// system in cm.
//

// Includes from C++
#include <vector>

// Includes from Mu2e
#include "TargetGeom/inc/TargetFoil.hh"
#include "GeometryService/inc/Detector.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class Target : public Detector{

  friend class TargetMaker;

  public:
    Target();
    ~Target(){;};

    virtual std::string name() const { return "Target";}
    
    int nFoils() const { return _foils.size(); }

    TargetFoil const& foil( unsigned int n ) const { return _foils.at(n); }

    double cylinderRadius() const {return _radius;};
    double cylinderLength() const {return _zLen;};
    double cylinderCenter() const {return _z0;};

    std::string const fillMaterial() const {return _fillMaterial;};
    
  protected:

    // All dimensions in mm.
    std::vector<TargetFoil> _foils;

    // an enclosing cylinder.  This is defined to be on the z axis, in
    // detector coordinates
    double _radius;
    double _z0; // center in Z;
    double _zLen; // Length

    std::string _fillMaterial;
    
    
};

}
#endif /* TargetGeom_Target_hh */
