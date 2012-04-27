// $Id: classes.h,v 1.2 2012/04/27 05:37:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/27 05:37:32 $
//
// Original author Andrei Gaponenko
//

#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::PSShield>;
template class art::Wrapper<mu2e::PSEnclosure>;

template class std::vector<mu2e::PSShield::Groove>;
