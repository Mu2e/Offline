#include "Offline/DbService/inc/RunTool.hh"

#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include "Offline/GeneralUtilities/inc/splitString.hh"

#include "cetlib_except/exception.h"

#include <algorithm>
#include <iomanip>

using namespace mu2e;

//**************************************************

RunTool::RunTool() {
  DbIdList idList;  // read file of db connection details
  _reader.setDbId(idList.getDbId("run_info"));
  _reader.setVerbose(0);
}

//**************************************************
std::map<int, std::string> RunTool::flags(FlagType ftype) {
  std::map<int, std::string> flags;
  // run query url, answer returned in csv
  int rc;
  std::string tcsv, select, table, order;
  StringVec where;

  select = "";
  if (ftype == FlagType::run) {
    table = "online.run_type";
  } else if (ftype == FlagType::transition) {
    table = "online.run_transition_type";
  }
  // where.push_back("name:eq:" + name);
  rc = _reader.query(tcsv, select, table, where, order);
  if (rc != 0 || tcsv.empty()) {
    throw cet::exception("RUNTOOL_BAD_FLAG_QUERY")
        << " RunTool::get failed to retrieve flags for " << table << "\n";
  }

  StringVec rows = splitString(tcsv, "\n", "\"", "\\", true, false);
  for (auto const& csv : rows) {
    StringVec sv = splitString(csv);
    // Table columns are: id, description, name
    // We want the "name" column (index 2) for transitions, not description
    // (index 1)
    if (sv.size() >= 3) {
      flags.insert({std::stoi(sv[0]), sv[2]});  // Use name column
    } else if (sv.size() >= 2) {
      flags.insert({std::stoi(sv[0]),
                    sv[1]});  // Fallback to description if name not available
    }
  }
  return flags;
}

//**************************************************
RunInfo::RunVec RunTool::listRuns(const RunSelect& runsel, bool configs,
                                  bool transitions, bool subruns) {
  RunInfo::RunVec runv;
  std::string tcsv, select, table, order;
  StringVec where;

  // fetch the list of all runs from the run table
  select = "run_number,comment,create_time,run_type_id";
  table = "online.run";
  order = "-run_number";
  int rc = _reader.query(tcsv, select, table, where, order);
  if (rc != 0 || tcsv.empty()) {
    throw cet::exception("RUNTOOL_BAD_RUN_QUERY")
        << " RunTool::listRuns failed to retrieve runs \n";
  }

  StringVec rows = DbUtil::splitCsvLines(tcsv);

  // apply the selection and fetch configs, transitions, and subruns if
  // requested

  // setup run cuts
  std::string rr = runsel.run();
  int rmin = 0;
  int rmax = 999999;
  StringVec rlims = splitString(rr, "-");
  if (rlims.size() > 0 && rlims[0].size() > 0) {
    rmin = std::stoi(rlims[0]);
  }
  if (rlims.size() > 1 && rlims[1].size() > 0) {
    rmax = std::stoi(rlims[1]);
  }
  if (rlims.size() == 1) {
    rmax = rmin;
  }

  // setup time cuts
  std::string tt = runsel.time();
  std::time_t tmin = 0;
  std::time_t tmax = 10000000000;
  StringVec tlims = splitString(tt, "/");
  if (tlims.size() > 0 && tlims[0].size() > 0) {
    rc = TimeUtility::parseTimeTZ(tlims[0], tmin);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << tlims[0] << " \n";
    }
  }
  if (tlims.size() > 1 && tlims[1].size() > 0) {
    rc = TimeUtility::parseTimeTZ(tlims[1], tmax);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << tlims[1] << " \n";
    }
  }

  // the list of acceptable run types
  std::string ty = runsel.type();
  StringVec tys = splitString(ty, ",");
  std::vector<int> itypes;
  for (auto t : tys) itypes.emplace_back(std::stoi(t));

  std::time_t dtime = 0;
  if (runsel.days() > 0) {
    dtime = std::time(0) - runsel.days() * 24 * 60 * 60;
  }

  // loop over all runs, starting with highest run number
  for (auto const& csv : rows) {
    // apply cuts on this run
    bool pass = true;

    // Parse the run table row: run_number,comment,create_time,run_type_id
    StringVec sv = DbUtil::splitCsv(csv);

    int run = std::stoi(sv[0]);
    std::string comment = sv[1];
    std::string create_time = sv[2];
    int run_type_id = std::stoi(sv[3]);

    // limit on the range
    if (run < rmin || run > rmax) pass = false;
    // limit on the number of runs
    if (runsel.last() > 0 && runv.size() >= size_t(runsel.last())) pass = false;

    std::time_t rtime;
    rc = TimeUtility::parseTimeTZ(create_time, rtime);
    if (rc != 0) {
      throw cet::exception("RUNTOOL_BAD_TIME")
          << " RunTool::listRuns can't parse time" << create_time << " \n";
    }
    // limit on time range
    if (rtime < tmin || rtime > tmax) pass = false;
    // limit on daysage
    if (rtime < dtime) pass = false;
    // limits on run types
    if (itypes.size() > 0) {
      if (std::find(itypes.begin(), itypes.end(), run_type_id) == itypes.end())
        pass = false;
    }

    if (pass) {
      RunInfo runInfo;

      // Set the main run record
      RunRecord runRecord;
      runRecord.setRunNumber(run);
      runRecord.setComment(comment);
      runRecord.setCreateTime(create_time);
      runRecord.setRunTypeId(run_type_id);
      runInfo.setRunRecord(runRecord);

      // Fetch config records if requested
      if (configs) {
        tcsv.clear();
        select = "run_number,subsystem,settings,create_time,version";
        table = "online.config";
        where.clear();
        std::string ww = std::string("run_number:eq:") + std::to_string(run);
        where.emplace_back(ww);
        order = "subsystem";
        rc = _reader.query(tcsv, select, table, where, order);
        if (rc != 0) {
          throw cet::exception("RUNTOOL_BAD_CONFIG_QUERY")
              << " RunTool::listRuns failed to retrieve configs for run " << run
              << "\n";
        }

        if (!tcsv.empty()) {
          StringVec crows = splitString(tcsv, "\n", "", "", true, false);
          for (auto const& crow : crows) {
            StringVec cfields = DbUtil::splitCsv(crow, true);
            RunConfig config;
            config.setRunNumber(std::stoi(cfields[0]));
            config.setSubsystem(cfields[1]);
            std::string json = DbUtil::simplfyQeString(cfields[2]);
            config.setSettings(json);
            config.setCreateTime(cfields[3]);
            config.setVersion(std::stoi(cfields[4]));
            runInfo.addConfig(config);
          }
        }
      }  // end configs

      // Fetch transition records if requested
      if (transitions) {
        tcsv.clear();
        select = "run_number,type_id,cause_id,transition_time";
        table = "online.run_transition";
        where.clear();
        std::string ww = std::string("run_number:eq:") + std::to_string(run);
        where.emplace_back(ww);
        order = "transition_time";
        rc = _reader.query(tcsv, select, table, where, order);
        if (rc != 0) {
          throw cet::exception("RUNTOOL_BAD_TRANSITION_QUERY")
              << " RunTool::listRuns failed to retrieve transitions for run "
              << run << "\n";
        }

        if (!tcsv.empty()) {
          StringVec transRows =
              splitString(tcsv, "\n", "\"", "\\", true, false);
          for (auto const& transCsv : transRows) {
            StringVec transSv = splitString(transCsv);
            if (transSv.size() >= 3) {
              RunTransition trans;
              trans.setRunNumber(std::stoi(transSv[0]));
              trans.setTypeId(std::stoi(transSv[1]));
              trans.setCauseId((transSv.size() > 2 && transSv[2] != "None")
                                   ? std::stoi(transSv[2])
                                   : 0);
              if (transSv.size() > 3) trans.setTransitionTime(transSv[3]);
              runInfo.addTransition(trans);
            }
          }
        }
      }  // end transitions

      // Fetch subrun records if requested
      if (subruns) {
        tcsv.clear();
        select =
            "run_number,subrun_number,n_events,n_on_spill,n_off_spill,min_ewt,"
            "max_ewt,start_time_unix,stop_time_unix,event_mode_counts,created_"
            "at,n_null";
        table = "online.subrun";
        where.clear();
        std::string ww = std::string("run_number:eq:") + std::to_string(run);
        where.emplace_back(ww);
        order = "subrun_number";
        rc = _reader.query(tcsv, select, table, where, order);
        if (rc != 0) {
          throw cet::exception("RUNTOOL_BAD_SUBRUN_QUERY")
              << " RunTool::listRuns failed to retrieve subruns for run " << run
              << "\n";
        }

        if (!tcsv.empty()) {
          StringVec subrunRows =
              splitString(tcsv, "\n", "\"", "\\", true, false);
          for (auto const& subrunCsv : subrunRows) {
            StringVec subrunSv = splitString(subrunCsv);
            if (subrunSv.size() >= 2) {
              RunSubRun subrun;
              subrun.setRunNumber(std::stoi(subrunSv[0]));
              subrun.setSubrunNumber(std::stoi(subrunSv[1]));
              if (subrunSv.size() > 2 && subrunSv[2] != "None")
                subrun.setNEvents(std::stol(subrunSv[2]));
              if (subrunSv.size() > 3 && subrunSv[3] != "None")
                subrun.setNOnSpill(std::stol(subrunSv[3]));
              if (subrunSv.size() > 4 && subrunSv[4] != "None")
                subrun.setNOffSpill(std::stol(subrunSv[4]));
              if (subrunSv.size() > 5 && subrunSv[5] != "None")
                subrun.setMinEwt(std::stol(subrunSv[5]));
              if (subrunSv.size() > 6 && subrunSv[6] != "None")
                subrun.setMaxEwt(std::stol(subrunSv[6]));
              if (subrunSv.size() > 7 && subrunSv[7] != "None")
                subrun.setStartTimeUnix(std::stoi(subrunSv[7]));
              if (subrunSv.size() > 8 && subrunSv[8] != "None")
                subrun.setStopTimeUnix(std::stoi(subrunSv[8]));
              if (subrunSv.size() > 9) subrun.setEventModeCounts(subrunSv[9]);
              if (subrunSv.size() > 10) subrun.setCreatedAt(subrunSv[10]);
              if (subrunSv.size() > 11 && subrunSv[11] != "None")
                subrun.setNNull(std::stol(subrunSv[11]));
              runInfo.addSubrun(subrun);
            }
          }
        }
      }  // end subruns

      runv.emplace_back(runInfo);
    }  // end if pass
  }

  return runv;
}

//**************************************************
void RunTool::printRun(const RunInfo& info) {
  const RunRecord& runRecord = info.runRecord();
  std::string create_time = runRecord.createTime();
  if (create_time.size() > 19) {
    create_time = create_time.substr(0, 19);
  }

  // Print basic run info
  std::cout << runRecord.runNumber();
  std::cout << std::setw(5) << runRecord.runTypeId();
  std::cout << std::setw(23) << create_time;
  std::cout << "\n";

  // Print configs if present
  if (info.configs().size() > 0) {
    std::cout << "  configs (" << info.configs().size() << "):\n";
    for (const auto& config : info.configs()) {
      std::cout << std::setw(12) << config.subsystem() << std::setw(3)
                << config.version() << "   " << config.createTime() << "\n";
    }
  }

  // Print transitions if present with column headers
  if (info.transitions().size() > 0) {
    // Get transition type names
    std::map<int, std::string> transitionTypes = flags(FlagType::transition);

    std::cout << "  transitions (" << info.transitions().size() << "):\n";
    std::cout << "    " << std::setw(8) << "type_id" << std::setw(20)
              << "type_name"
              << "  transition_time\n";
    for (const auto& trans : info.transitions()) {
      std::string type_name = "";
      auto it = transitionTypes.find(trans.typeId());
      if (it != transitionTypes.end()) {
        type_name = it->second;
      }

      std::cout << "    " << std::setw(8) << trans.typeId() << std::setw(20)
                << type_name << "  " << trans.transitionTime() << "\n";
    }
  }

  // Print subruns if present with column headers
  if (info.subruns().size() > 0) {
    std::cout << "  subruns (" << info.subruns().size() << "):\n";
    std::cout << "    " << std::setw(8) << "subrun" << std::setw(12) << "events"
              << std::setw(15) << "min_ewt" << std::setw(15) << "max_ewt"
              << "  start_time                stop_time\n";

    // Save current TZ environment variable to restore later
    const char* old_tz = std::getenv("TZ");
    std::string saved_tz = old_tz ? old_tz : "";

    // Set timezone to Central US (America/Chicago) once before the loop
    setenv("TZ", "America/Chicago", 1);
    tzset();

    for (const auto& subrun : info.subruns()) {
      // Convert Unix timestamps to ISO8601 format in Central US time
      std::time_t start_t = subrun.startTimeUnix();
      std::time_t stop_t = subrun.stopTimeUnix();
      char start_buf[32], stop_buf[32];

      std::strftime(start_buf, sizeof(start_buf), "%Y-%m-%d %H:%M:%S",
                    std::localtime(&start_t));
      std::strftime(stop_buf, sizeof(stop_buf), "%Y-%m-%d %H:%M:%S",
                    std::localtime(&stop_t));

      std::cout << "    " << std::setw(8) << subrun.subrunNumber()
                << std::setw(12) << subrun.nEvents() << std::setw(15)
                << subrun.minEwt() << std::setw(15) << subrun.maxEwt() << "  "
                << start_buf << "  " << stop_buf << "\n";
    }

    // Restore original TZ environment variable
    if (!saved_tz.empty()) {
      setenv("TZ", saved_tz.c_str(), 1);
    } else {
      unsetenv("TZ");
    }
    tzset();
  }
}
