// $Id: classes.h,v 1.3 2012/06/06 19:29:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/06/06 19:29:31 $
//
// Original author Andrei Gaponenko
//

#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::PSVacuum>;
template class art::Wrapper<mu2e::PSShield>;
template class art::Wrapper<mu2e::PSEnclosure>;

template class std::vector<mu2e::PSShield::Groove>;
