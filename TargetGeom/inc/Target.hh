#ifndef TargetGeom_Target_hh
#define TargetGeom_Target_hh

//
// Class to represent the system of target foils.
// For now these are just disks perpendicular to the z axis.
//
// $Id: Target.hh,v 1.8 2012/02/24 20:55:48 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 20:55:48 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate
// system in cm.
//

// Includes from C++
#include <vector>

// Includes from Mu2e
#include "Mu2eInterfaces/inc/Detector.hh"
#include "TargetGeom/inc/TargetFoil.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class Target : virtual public Detector{

  friend class TargetMaker;

  public:
    Target();

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int nFoils() const { return _foils.size(); }

    TargetFoil const& foil( unsigned int n ) const { return _foils.at(n); }

    double cylinderRadius() const {return _radius;}
    double cylinderLength() const {return _zLen;}
    double cylinderCenter() const {return _z0;}

    std::string const fillMaterial() const {return _fillMaterial;}

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
