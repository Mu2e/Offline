#ifndef TrackerGeom_Support_hh
#define TrackerGeom_Support_hh

//
// Describe the properites of a support for the Tracker.
//
//
//  $Id: Support.hh,v 1.5 2012/03/30 15:13:35 gandr Exp $
//  $Author: gandr $
//  $Date: 2012/03/30 15:13:35 $
//
//  Original author Rob Kutschke
//

#include <string>

#include "GeomPrimitives/inc/TubsParams.hh"

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

#endif /* TrackerGeom_Support_hh */
