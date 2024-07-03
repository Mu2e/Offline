#include "Offline/DbTables/inc/DbTableFactory.hh"
#include "Offline/DbTables/inc/AnaTrkQualDb.hh"
#include "Offline/DbTables/inc/CRVBadChan.hh"
#include "Offline/DbTables/inc/CRVPhoton.hh"
#include "Offline/DbTables/inc/CRVSiPM.hh"
#include "Offline/DbTables/inc/CRVTime.hh"
#include "Offline/DbTables/inc/STMEnergyPar.hh"
#include "Offline/DbTables/inc/STMPedestals.hh"
#include "Offline/DbTables/inc/SimEfficiencies.hh"
#include "Offline/DbTables/inc/SimEfficiencies2.hh"

#include "Offline/DbTables/inc/CalSourceEnergyCalib.hh"
#include "Offline/DbTables/inc/CalCosmicEnergyCalib.hh"
#include "Offline/DbTables/inc/CalCosmicTimeCalib.hh"
#include "Offline/DbTables/inc/CalLaserEnergyCalib.hh"
#include "Offline/DbTables/inc/CalCosmicTimeCalib.hh"
#include "Offline/DbTables/inc/CalLaserTimeCalib.hh"
#include "Offline/DbTables/inc/CalLaserEnergyCalib.hh"
#include "Offline/DbTables/inc/CalEnergyCalib.hh"
#include "Offline/DbTables/inc/CalTimeCalib.hh"
#include "Offline/DbTables/inc/CalCosmicT0Align.hh"

#include "Offline/DbTables/inc/TrkAlignElement.hh"
#include "Offline/DbTables/inc/TrkAlignElementSim.hh"
#include "Offline/DbTables/inc/TrkAlignStraw.hh"
#include "Offline/DbTables/inc/TrkAlignStrawSim.hh"
#include "Offline/DbTables/inc/TrkDelayPanel.hh"
#include "Offline/DbTables/inc/TrkDelayRStraw.hh"
#include "Offline/DbTables/inc/TrkElementStatus.hh"
#include "Offline/DbTables/inc/TrkPreampStraw.hh"
#include "Offline/DbTables/inc/TstCalib1.hh"
#include "Offline/DbTables/inc/TstCalib2.hh"
#include "Offline/DbTables/inc/TstCalib3.hh"
#include "Offline/DbTables/inc/TstAdhoc1.hh"
#include "cetlib_except/exception.h"


mu2e::DbTable::ptr_t mu2e::DbTableFactory::newTable(std::string const& name) {
  if (name == "TstCalib1") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib1());
  } else if (name == "TstCalib2") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib2());
  } else if (name == "TstCalib3") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstCalib3());
  } else if (name == "TstAdhoc1") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TstAdhoc1());
  } else if (name == "TrkDelayPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkDelayPanel());
  } else if (name == "TrkDelayRStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkDelayRStraw());
  } else if (name == "TrkPreampStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPreampStraw());
  } else if (name == "TrkAlignTracker") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignTracker());
  } else if (name == "TrkAlignPlane") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPlane());
  } else if (name == "TrkAlignPanel") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPanel());
  } else if (name == "TrkAlignStraw") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignStraw());
  } else if (name == "TrkAlignTrackerSim") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignTrackerSim());
  } else if (name == "TrkAlignPlaneSim") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPlaneSim());
  } else if (name == "TrkAlignPanelSim") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignPanelSim());
  } else if (name == "TrkAlignStrawSim") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkAlignStrawSim());
  } else if (name == "TrkPlaneStatus") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPlaneStatus());
  } else if (name == "TrkPanelStatus") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkPanelStatus());
  } else if (name == "TrkStrawStatusLong") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkStrawStatusLong());
  } else if (name == "TrkStrawStatusShort") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::TrkStrawStatusShort());
  } else if (name == "AnaTrkQualDb") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::AnaTrkQualDb());
  } else if (name == "SimEfficiencies") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::SimEfficiencies());
  } else if (name == "SimEfficiencies2") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::SimEfficiencies2());
  } else if (name == "STMEnergyPar") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::STMEnergyPar());
  } else if (name == "STMPedestals") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::STMPedestals());
  } else if (name == "CRVBadChan") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CRVBadChan());
  } else if (name == "CRVPhoton") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CRVPhoton());
  } else if (name == "CRVSiPM") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CRVSiPM());
  } else if (name == "CRVTime") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CRVTime());
  }  else if (name=="CalSourceEnergyCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalSourceEnergyCalib());
  }  else if (name=="CalCosmicTimeCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalCosmicTimeCalib());
  } else if (name=="CalCosmicEnergyCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalCosmicEnergyCalib());
  } else if (name=="CalLaserEnergyCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalLaserEnergyCalib());
  } else if (name=="CalLaserTimeCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalLaserTimeCalib());
  } else if (name=="CalEnergyCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalEnergyCalib());
  } else if (name=="CalTimeCalib") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalTimeCalib());
  } else if (name=="CalCosmicT0Align") {
    return std::shared_ptr<mu2e::DbTable>(new mu2e::CalCosmicT0Align());
  }else {
    throw cet::exception("DBFILE_BAD_TABLE_NAME")
        << "DbTableFactory::newTable call with bad table name: " + name + "\n";
  }
}
