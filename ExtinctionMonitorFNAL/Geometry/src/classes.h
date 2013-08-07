// $Id: classes.h,v 1.4 2013/08/07 20:42:11 wieschie Exp $
// $Author: wieschie $
// $Date: 2013/08/07 20:42:11 $
//
// Original author Andrei Gaponenko
//

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlane.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtMonFNALBuilding>;
template class art::Wrapper<mu2e::ExtMonFNAL::ExtMon>;
