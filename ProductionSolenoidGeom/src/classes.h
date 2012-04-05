// $Id: classes.h,v 1.1 2012/04/05 18:43:39 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/05 18:43:39 $
//
// Original author Andrei Gaponenko
//

#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::PSShield>;
template class art::Wrapper<mu2e::PSEnclosure>;
