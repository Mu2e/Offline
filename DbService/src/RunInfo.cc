#include "Offline/DbService/inc/RunInfo.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"
#include "cetlib_except/exception.h"

using namespace mu2e;

//**************************************************

RunInfo::RunInfo(std::string csv, std::string conditions,
                 std::string transitions) {
  //  run_number,run_type,condition_id,artdaq_partition,host_name,configuration_name,configuration_version,context_name,context_version,trigger_table_name,trigger_table_version,online_software_version,commit_time,shifter_note
  // 1,4,7,14,mu2edaq13.fnal.gov,dcsConfig,53,dcsContext,52,None,None,None,2024-05-06 16:17:21.973417-05:00,None
  _csv = csv;
  StringVec sv = splitString(csv);
  if (sv.size() != 14) {
    throw cet::exception("RUNTOOL_BAD_INFO_CSV")
        << " RunInfo: failed to construct from " << sv.size() << "fields in\n"
        << csv << "\n";
  }

  _run_number = std::stoi(sv[0]);
  _run_type = std::stoi(sv[1]);
  _condition_id = std::stoi(sv[2]);
  _artdaq_partition = std::stoi(sv[3]);
  _host_name = sv[4];
  _configuration_name = sv[5];
  _configuration_version = std::stoi(sv[6]);
  _context_name = sv[7];
  _context_version = std::stoi(sv[8]);
  _trigger_table_name = sv[9];
  if(sv[10]=="None") {
    _trigger_table_version = 0;
  } else {
    _trigger_table_version = std::stoi(sv[10]);
  }
  _online_software_version = sv[11];
  _commit_time = sv[12];
  _shifter_note = sv[13];

  _conditions = conditions;
  _transitions = transitions;
}
