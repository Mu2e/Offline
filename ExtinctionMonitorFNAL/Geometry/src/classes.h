// $Id: classes.h,v 1.2 2012/08/28 05:02:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/28 05:02:31 $
//
// Original author Andrei Gaponenko
//

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtMonFNALBuilding>;
template class art::Wrapper<mu2e::ExtMonFNAL::ExtMon>;
