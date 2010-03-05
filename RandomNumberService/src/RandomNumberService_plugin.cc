//
// Plugin for the Random Number Service.
//
// $Id: RandomNumberService_plugin.cc,v 1.2 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/05 16:07:38 $
//
// Original author Rob Kutschke
//
// Notes
// 1) This is a special service; it has the usual service interface plus it has
//    additional interface, RandomNumberGenerator, that is used by the framework and IOModules.
//    Therefore this service must be known by the name of the base class, not the name
//    of the concrete class.
//
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "RandomNumberService/inc/RandomNumberService.hh"

using mu2e::RandomNumberService;

typedef edm::serviceregistry::AllArgsMaker<edm::RandomNumberGenerator,RandomNumberService> RandomMaker;
DEFINE_FWK_SERVICE_MAKER(RandomNumberService, RandomMaker);
