#include "cetlib_except/exception.h"

#include "Offline/DbTables/inc/DbTableFactory.hh"
#include "Offline/DbTables/inc/TstCalib1.hh"
#include "Offline/DbTables/inc/TstCalib2.hh"
#include "Offline/DbTables/inc/TstCalib3.hh"
#include "Offline/DbTables/inc/TrkDelayPanel.hh"
#include "Offline/DbTables/inc/TrkDelayRStraw.hh"
#include "Offline/DbTables/inc/TrkPreampStraw.hh"
#include "Offline/DbTables/inc/TrkAlignElement.hh"
#include "Offline/DbTables/inc/TrkAlignStraw.hh"
#include "Offline/DbTables/inc/TrkElementStatus.hh"
#include "Offline/DbTables/inc/AnaTrkQualDb.hh"
#include "Offline/DbTables/inc/CalRoIDMapDIRACToOffline.hh"
#include "Offline/DbTables/inc/CalRoIDMapOfflineToDIRAC.hh"

#include "Offline/DbTables/inc/SimEfficiencies.hh"

mu2e::DbTable::ptr_t mu2e::DbTableFactory::newTable(std::string const& name) {
  if (name=="TstCalib1") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib1());
  } else if (name=="TstCalib2") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib2());
  } else if (name=="TstCalib3") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib3());
  } else if (name=="TrkDelayPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkDelayPanel());
  } else if (name=="TrkDelayRStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkDelayRStraw());
  } else if (name=="TrkPreampStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPreampStraw());
  } else if (name=="TrkAlignTracker") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignTracker());
  } else if (name=="TrkAlignPlane") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPlane());
  } else if (name=="TrkAlignPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPanel());
  } else if (name=="TrkAlignStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignStraw());
  } else if (name=="TrkPlaneStatus") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPlaneStatus());
  } else if (name=="TrkPanelStatus") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPanelStatus());
  } else if (name=="TrkStrawStatusLong") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkStrawStatusLong());
  } else if (name=="TrkStrawStatusShort") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkStrawStatusShort());
  } else if (name=="AnaTrkQualDb") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::AnaTrkQualDb());
  } else if (name=="SimEfficiencies") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::SimEfficiencies());
  } else if (name=="CalRoIDMapDIRACToOffline") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalRoIDMapDIRACToOffline());
  } else if (name=="CalRoIDMapOfflineToDIRAC") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalRoIDMapOfflineToDIRAC());    
  } else {
    throw cet::exception("DBFILE_BAD_TABLE_NAME")
      << "DbTableFactory::newTable call with bad table name: "+name+"\n";

  }
}
