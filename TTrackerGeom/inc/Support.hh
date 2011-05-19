#ifndef TTrackerGeom_Support_hh
#define TTrackerGeom_Support_hh

//
// Describe the properites of a support for the TTracker.
//
//
//  $Id: Support.hh,v 1.4 2011/05/19 22:23:06 wb Exp $
//  $Author: wb $
//  $Date: 2011/05/19 22:23:06 $
//
//  Original author Rob Kutschke
//

#include <string>

#include "TrackerGeom/inc/TubsParams.hh"

namespace mu2e {

  class Support{

  public:

    Support():
      innerRadius_(0.),
      outerRadius_(0.),
      halfThickness_(0.),
      materialName_(){}

    Support( double inRad, double outRad, double halfThick, const std::string& name):
      innerRadius_(inRad),
      outerRadius_(outRad),
      halfThickness_(halfThick),
      materialName_(name){}

    // Accept compiler supplied destructor, copy c'tor and assignment
    // operator.

    // Return the information formatted to make a G4Tubs.
    TubsParams getTubsParams() const{
      return TubsParams(innerRadius_, outerRadius_, halfThickness_);
    }

    double innerRadius() const { return innerRadius_; }
    double outerRadius() const { return outerRadius_; }
    double halfThickness() const { return halfThickness_; }
    std::string materialName() const { return materialName_; }

  private:

    // Shape parameters of the support:
    // Inner and outer radii and half of thickness.
    double innerRadius_;
    double outerRadius_;
    double halfThickness_;

    // Name of the G4material that makes up the support.
    std::string materialName_;

  };

}

#endif /* TTrackerGeom_Support_hh */
