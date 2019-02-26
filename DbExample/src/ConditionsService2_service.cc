//
//
//
#include <iostream>
#include <typeinfo>
#include "DbService/inc/DbService.hh"
#include "DbExample/inc/ConditionsService2.hh"

#include "DbExample/inc/DetData1Cache.hh"
#include "DbExample/inc/DetData2Cache.hh"
#include "DbExample/inc/DetData3Cache.hh"

using namespace std;

namespace mu2e {

  ConditionsService2::ConditionsService2(Parameters const& sTable,
                                       art::ActivityRegistry& iRegistry) :
    _config(sTable()),_simpleConfig(sTable().fileName()),
    _verbose(sTable().verbose())
  {

    // create this here to force DbService to be active before Conditions2
    art::ServiceHandle<DbService> d;

    auto dd1 = std::make_shared<mu2e::DetData1Cache>();
    _caches[dd1->name()] = dd1;
    auto dd2 = std::make_shared<mu2e::DetData2Cache>();
    _caches[dd2->name()] = dd2;
    auto dd3 = std::make_shared<mu2e::DetData3Cache>();
    _caches[dd3->name()] = dd3;

    //iRegistry.sPostBeginJob.watch(this, &ConditionsService2::postBeginJob);
  }

  //void ConditionsService2::postBeginJob() {
  //}


  //ConditionsService2::~ConditionsService2(){
  //}

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::ConditionsService2);
