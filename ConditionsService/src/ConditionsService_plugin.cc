//
// Plugin for the Conditions Service.
//
// $Id: ConditionsService_plugin.cc,v 1.2 2011/05/17 15:35:59 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:35:59 $
//
// Original author Rob Kutschke
//

#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "ConditionsService/inc/ConditionsService.hh"

using mu2e::ConditionsService;
DEFINE_ART_SERVICE(ConditionsService);
