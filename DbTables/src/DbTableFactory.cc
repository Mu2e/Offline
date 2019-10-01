#include "cetlib_except/exception.h"

#include "DbTables/inc/DbTableFactory.hh"
#include "DbTables/inc/TstCalib1.hh"
#include "DbTables/inc/TstCalib2.hh"
#include "DbTables/inc/TstCalib3.hh"
#include "DbTables/inc/TrkDelayPanel.hh"
#include "DbTables/inc/TrkPreampRStraw.hh"
#include "DbTables/inc/TrkPreampStraw.hh"
#include "DbTables/inc/TrkThresholdRStraw.hh"
#include "DbTables/inc/TrkAlignTracker.hh"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "DbTables/inc/TrkAlignPanel.hh"

mu2e::DbTable::ptr_t mu2e::DbTableFactory::newTable(std::string const& name) {
  if (name=="TstCalib1") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib1());
  } else if (name=="TstCalib2") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib2());
  } else if (name=="TstCalib3") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib3());
  } else if (name=="TrkDelayPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkDelayPanel());
  } else if (name=="TrkPreampRStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPreampRStraw());
  } else if (name=="TrkPreampStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPreampStraw());
  } else if (name=="TrkThresholdRStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkThresholdRStraw());
  } else if (name=="TrkAlignTracker") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignTracker());
  } else if (name=="TrkAlignPlane") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPlane());
  } else if (name=="TrkAlignPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPanel());
  } else {
    throw cet::exception("DBFILE_BAD_TABLE_NAME") 
      << "DbTableFactory::newTable call with bad table name: "+name+"\n";
    
  }
}

