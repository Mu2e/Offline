//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.cc,v 1.1 2011/05/10 16:44:02 vrusu Exp $ 
// $Author: vrusu $
// $Date: 2011/05/10 16:44:02 $
// 

// Mu2e include files
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  TrackerCalibrations::TrackerCalibrations( SimpleConfig const& config ){
    
    // Here we should eventually interface to some database

    
  }
  
  const double TrackerCalibrations::TimeDivisionResolution(StrawIndex strawIndex, double znorm) const {
    double resopar0 = 64.2;
    double resopar1 = 60.7;

    double reso  = resopar0 + resopar1 * (znorm - 0.5) * (znorm - 0.5); //resolution in mm
    return reso;

  }

  const double TrackerCalibrations::SignalVelocity(StrawIndex strawIndex) const {
    double distvsdeltat = 231.;
    return distvsdeltat; //mm/ns
  }

  const double TrackerCalibrations::TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const{ 
    return 0.5 * SignalVelocity(strawIndex) * deltaT; 
  }
}
