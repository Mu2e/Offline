// $Id: classes.h,v 1.2 2012/04/25 18:19:14 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/25 18:19:14 $
//
// Original author Andrei Gaponenko
//

#include "Mu2eBuildingGeom/inc/BuildingBasics.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "art/Persistency/Common/Wrapper.h"

template class art::Wrapper<mu2e::BuildingBasics>;
template class art::Wrapper<mu2e::Mu2eBuilding>;
