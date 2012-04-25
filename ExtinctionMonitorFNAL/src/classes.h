// $Id: classes.h,v 1.2 2012/04/25 18:47:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/25 18:47:32 $
//
// Original author Andrei Gaponenko
//

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtMonFNALBuilding>;
template class art::Wrapper<mu2e::ExtMonFNAL::ExtMon>;
