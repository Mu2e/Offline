//
// Plugin for the Random Number Service.
//
// $Id: RandomNumberService_plugin.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "RandomNumberService/inc/RandomNumberService.hh"

using mu2e::RandomNumberService;

DEFINE_FWK_SERVICE(RandomNumberService);
