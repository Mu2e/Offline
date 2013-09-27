// $Id: classes.h,v 1.1 2013/09/27 13:07:05 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/09/27 13:07:05 $
//
// Original author Andrei Gaponenko
//

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1a.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1b.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::ExtNeutShieldUpstream1a>;
template class art::Wrapper<mu2e::ExtNeutShieldUpstream1b>;
