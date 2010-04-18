#ifndef SUPPORT_HH
#define SUPPORT_HH

//
// Describe the properites of a support for the TTracker.
//
//
//  $Id: Support.hh,v 1.1 2010/04/18 00:37:16 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/04/18 00:37:16 $
//
//  Original author Rob Kutschke
//

#include <string>

#include "TrackerGeom/inc/TubsParams.hh"

namespace mu2e {

  struct Support{

    Support():
      innerRadius(0.),
      outerRadius(0.),
      halfThickness(0.),
      materialName(){}

    Support( double inRad, double outRad, double halfThick, const std::string& name):
      innerRadius(inRad),
      outerRadius(outRad),
      halfThickness(halfThick),
      materialName(name){}

    // Accept compiler supplied destructor, copy c'tor and assignment
    // operator.

    // Shape parameters of the support:
    // Inner and outer radii and half of thickness.
    double innerRadius;
    double outerRadius;
    double halfThickness;

    // Name of the G4material that makes up the support.
    std::string materialName;

    // Return the information formatted to make a G4Tubs.
    TubsParams getTubsParams() const{
      return TubsParams(innerRadius, outerRadius, halfThickness);
    }
    
  };

}

#endif
