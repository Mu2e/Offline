#ifndef TargetGeom_Target_hh
#define TargetGeom_Target_hh

//
// Class to represent the system of target foils.
// For now these are just disks perpendicular to the z axis.
//
// $Id: Target.hh,v 1.11 2013/05/31 18:07:18 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:18 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate
// system in mm.
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
    Target() : _radius(), _zLen() {}

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int nFoils() const { return _foils.size(); }

    TargetFoil const& foil( unsigned int n ) const { return _foils.at(n); }

    double cylinderRadius() const {return _radius;}
    double cylinderLength() const {return _zLen;}
    const CLHEP::Hep3Vector& centerInMu2e() const {return _centerInMu2e;}

    std::string const fillMaterial() const {return _fillMaterial;}

  protected:

    // All dimensions in mm.
    std::vector<TargetFoil> _foils;

    // an enclosing cylinder.  This is defined to be on the z axis, in
    // detector coordinates
    double _radius;
    double _zLen; // Length
    CLHEP::Hep3Vector _centerInMu2e;

    std::string _fillMaterial;

  };
}
#endif /* TargetGeom_Target_hh */
