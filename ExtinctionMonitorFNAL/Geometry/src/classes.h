// $Id: classes.h,v 1.1 2012/08/24 15:06:59 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/24 15:06:59 $
//
// Original author Andrei Gaponenko
//

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtMonFNALBuilding>;
template class art::Wrapper<mu2e::ExtMonFNAL::ExtMon>;
