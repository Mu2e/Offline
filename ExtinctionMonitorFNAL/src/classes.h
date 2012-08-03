// $Id: classes.h,v 1.3 2012/08/03 00:32:10 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/03 00:32:10 $
//
// Original author Andrei Gaponenko
//

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnet.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALSensorStack.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtMonFNALBuilding>;
template class art::Wrapper<mu2e::ExtMonFNAL::ExtMon>;
